#include <MatFillerPCarithAvg.H>

MatFillerPCarithAvg::MatFillerPCarithAvg()
  : MatIDFiller()
{
}

MatFillerPCarithAvg::MatFillerPCarithAvg(const Array<Geometry>&  geomArray,
					 const Array<IntVect>&   refRatio,
					 const PArray<Material>& materials)
  : MatIDFiller(geomArray,refRatio,materials)
{
}

void 
MatFillerPCarithAvg::define(const Array<Geometry>& _geomArray,
			    const Array<IntVect>&  _refRatio,
			    const PArray<Material>& _materials)
{
  MatIDFiller::define(_geomArray,_refRatio,materials);
}

using std::cout;
using std::endl;
void writeMF(std::ostream& os, const MultiFab& mf) {
  const DistributionMapping& dm = mf.DistributionMap();
  for (int i=0; i<mf.size(); ++i) {
    if (ParallelDescriptor::MyProc()==dm[i]) {
      cout << mf[i] << endl;
    }
    ParallelDescriptor::Barrier();
  }
}

void 
MatFillerPCarithAvg::SetProperty(Real t, int level, MultiFab& mf, const std::string& pname, int dComp, int nGrow)
{
  bool ioproc = ParallelDescriptor::IOProcessor();
  BoxArray unfilled(mf.boxArray()); unfilled.grow(nGrow);
  
  if (level<ba_mixed.size() && ba_mixed[level].size()>0) {
    BoxArray baM = BoxLib::intersect(ba_mixed[level],unfilled); 
    if (baM.size()>0) {
      baM.removeOverlap();
      MultiFab mixed(baM,1,0);
      FillCoarseCells(t,level,mixed,pname,dComp);
      mf.copy(mixed);
      BoxList remaining;
      for (int i=0, N=unfilled.size(); i<N; ++i) {
	BoxArray bl = BoxLib::complementIn(unfilled[i],ba_mixed[level]);
	remaining.join(BoxList(bl));
      }
      remaining.simplify();
      unfilled = BoxArray(remaining);
    }
  }

  if (unfilled.ok()) {
    // Fill based on material ID at level
    BL_ASSERT(materialID[level].boxArray().isDisjoint());
    BoxList fillable = BoxLib::intersect((BoxList)materialID[level].boxArray(),(BoxList)unfilled);
    BoxList remaining;
    if (fillable.size()>0) {
      fillable.simplify();
      for (int j=0; j<unfilled.size(); ++j) {
	remaining.join(BoxLib::complementIn(unfilled[j],fillable));
      }
      remaining.simplify();
      BoxArray ba_fillable(fillable);
      MultiFab fillData(ba_fillable,1,0);
      MultiFab fillID(ba_fillable,1,0); fillID.setVal(1000);
      fillID.copy(materialID[level]);
      
      for (MFIter mfi(fillData); mfi.isValid(); ++mfi) {
	const Box& ovlp = mfi.validbox();
	const FArrayBox& idfab = fillID[mfi];
	FArrayBox& matfab = fillData[mfi];
	for (IntVect iv=ovlp.smallEnd(), End=ovlp.bigEnd(); iv<=End; ovlp.next(iv)) {
	  int matID = idfab(iv,0);
	  const Property* p = materials[matID].Prop(pname);
	  BL_ASSERT(p!=0);
	  p->eval(t,level,Box(iv,iv),matfab,dComp);
	}
      }
      mf.copy(fillData);
    }


    if (remaining.isNotEmpty()) {
      IntVect cumRatio = IntVect(D_DECL(1,1,1));
      FArrayBox cfab;
      for (int lev=level-1; lev>=0 && remaining.isNotEmpty(); --lev) {
        cumRatio *= RefRatio(lev);
        BoxArray bac = BoxArray(BoxList(remaining).coarsen(cumRatio));

	if (ioproc) {
	  cout << "Remaining coarsened (level,lev): " << level << ", " << lev << ": " << bac << endl;
	}

        IntVect refm = IntVect(cumRatio)-IntVect(D_DECL(1,1,1));
        Box rbox(IntVect(D_DECL(0,0,0)),refm);

        const BoxArray& bacd = materialID[lev].boxArray();
	BoxArray bac_interp = BoxLib::intersect(bacd,bac); bac_interp.removeOverlap();

	if (ioproc) {
	  cout << "Remaining coarsened that can be filled for level " << level << " from lev " << lev << ": " << bac_interp << endl;
	}

	MultiFab mfc_interp(bac_interp,1,0);
	SetProperty(t,lev,mfc_interp,pname,0,0);
	
	if (ioproc) {
	  cout << "phi for level " << level << " at lev " << lev << endl;
	}
	writeMF(cout,mfc_interp);


	BoxArray baf_interp(bac_interp); baf_interp.refine(cumRatio);
	MultiFab mff_interp(baf_interp,1,0);
	for (MFIter mfi(mfc_interp); mfi.isValid(); ++mfi) {
	  const Box& cbox = mfi.validbox();
	  const FArrayBox& cfab = mfc_interp[mfi];
	  FArrayBox& ffab = mff_interp[mfi];
	  for (IntVect civ=cbox.smallEnd(), End=cbox.bigEnd(); civ<=End; cbox.next(civ)) {
	    IntVect fiv = civ * cumRatio;
	    Box setBox = rbox; setBox.shift(fiv);
	    ffab.setVal(cfab(civ,0),setBox,0,1);
	  }
	}
	mf.copy(mff_interp,0,dComp,1);

	if (ioproc) {
	  cout << "interpolated phi for level " << level << " at lev " << lev << endl;
	}
	writeMF(cout,mff_interp);

	BoxList blcd(bacd);
	BoxList blc_remaining;
        for (int i=0; i<bac.size(); ++i) {
	  blc_remaining.join(BoxLib::complementIn(bac[i],blcd));
	}
	blc_remaining.simplify();
	BoxArray baf_remaining(blc_remaining); baf_remaining.refine(cumRatio);  unfilled=baf_remaining;


      }
      if (remaining.isNotEmpty()) {
	remaining.simplify();
      }
    }
  }
}

void
arithmetic_average(const FArrayBox& fineFab,
		   int              sComp,
		   FArrayBox&       crseFab,
		   const Box&       crseBox,
		   int              dComp,
		   int              nComp,
		   const IntVect&   ref)
{
  IntVect refm = IntVect(ref)-IntVect(D_DECL(1,1,1));
  Box lbox(IntVect(D_DECL(0,0,0)),refm);
  
  Real Ninv = 1.0 / lbox.numPts();
  for (IntVect civ = crseBox.smallEnd(); civ<=crseBox.bigEnd(); crseBox.next(civ)) {
    IntVect fbase = ref*civ;
    Box lfbox(fbase,fbase+refm);
    for (int n=0; n<nComp; ++n) {
      crseFab(civ,dComp+n) = fineFab.sum(lfbox,sComp+n,1) * Ninv;
    }
  }
}

void
MatFillerPCarithAvg::FillCoarseCells(Real t, int level, MultiFab& mfc, const std::string& pname, int dComp)
{
  int finestLevel = nLevs-1;
  if (level<finestLevel) {
    const BoxArray& bac = mfc.boxArray();
    const IntVect& ref = RefRatio(level);
    BoxArray baf = BoxArray(bac).refine(ref);
    MultiFab mff(baf,1,0);
    SetProperty(t,level+1,mff,pname,0,0);
    for (MFIter mfi(mfc); mfi.isValid(); ++mfi) {
      const Box& cbox = mfi.validbox();
      arithmetic_average(mff[mfi],0,mfc[mfi],cbox,0,1,ref);
    }
  }
}
