#include <MatFillerPC.H>

MatFillerPC::MatFillerPC()
  : MatIDFiller()
{
}

MatFillerPC::MatFillerPC(const Array<Geometry>&  geomArray,
                         const Array<IntVect>&   refRatio,
                         const PArray<Material>& materials)
  : MatIDFiller(geomArray,refRatio,materials)
{
}

void 
MatFillerPC::define(const Array<Geometry>& _geomArray,
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
MatFillerPC::SetProperty(Real t, int level, MultiFab& mf, const std::string& pname, int dComp, int nGrow)
{
  BoxArray unfilled(mf.boxArray()); unfilled.grow(nGrow);
  
  BoxArray baM;
  if (level<ba_mixed.size() && ba_mixed[level].size()>0) {
    baM = BoxLib::intersect(ba_mixed[level],unfilled); 
    if (baM.size()>0) {
      baM.removeOverlap();
      MultiFab mixed(baM,1,0);
      FillCoarseCells(t,level,mixed,pname,dComp); // here this level is "coarse", will average finer data
      mf.copy(mixed);
      BoxList bl_remaining;
      for (int i=0, N=unfilled.size(); i<N; ++i) {
	BoxArray bl = BoxLib::complementIn(unfilled[i],ba_mixed[level]);
	bl_remaining.join(BoxList(bl));
      }
      bl_remaining.simplify();
      unfilled = BoxArray(bl_remaining);
    }
  }


  // Here, unfilled split into baM(mixed)[treated above] and unfilled(nonmixed) 
  // Of the unfilled ones, first fill the ones we can at this level, then go to coarser

  if (unfilled.ok() && level<materialID.size()) {
    BL_ASSERT(materialID[level].boxArray().isDisjoint());

    // Find portion of unfilled that are fillable by nonmixed cells in the materialID struct at this level
    BoxList bl_fillable = BoxLib::intersect((BoxList)materialID[level].boxArray(),(BoxList)unfilled);
    BL_ASSERT( BoxLib::intersect((BoxArray)bl_fillable,baM).size()==0 );

    if (bl_fillable.size()>0) {
      bl_fillable.simplify();
      BoxArray ba_fillable(bl_fillable);
      MultiFab fillData(ba_fillable,1,0);
      MultiFab fillID(ba_fillable,1,0); 
      fillID.setVal(-2000); // FIME: remove
      fillID.copy(materialID[level]); // guaranteed to be filled completely
      
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

    // Find cells not filled at this level
    BoxList bl_not_filled;
    for (int i=0; i<unfilled.size(); ++i) {
      bl_not_filled.join(BoxLib::complementIn(unfilled[i],bl_fillable));
    }
    bl_not_filled.simplify();
    BoxArray remaining = (BoxArray)bl_not_filled;

    // If anything left, fill by interpolation of coarser
    if (remaining.size()!=0) {
      IntVect cumRatio = IntVect(D_DECL(1,1,1));
      FArrayBox cfab;
      for (int lev=level-1; lev>=0 && remaining.size()!=0; --lev) {
        cumRatio *= RefRatio(lev);
        if (materialID.size()>lev) {
          BoxArray bac_remaining = BoxArray(BoxList(remaining).coarsen(cumRatio));
          IntVect refm = IntVect(cumRatio)-IntVect(D_DECL(1,1,1));
          Box rbox(IntVect(D_DECL(0,0,0)),refm);
          
          // Get ovlp, including mixed and nonmixed cells, then remove mixed ones
          const BoxArray& ba_mat = materialID[lev].boxArray();
          BoxArray bac_interp = BoxLib::intersect(bac_remaining,ba_mat);
          if (lev<ba_mixed.size() && ba_mixed[lev].size()>0) {
            BoxList bl_mixed_local = (BoxList)BoxLib::intersect(bac_interp,ba_mixed[lev]);
            BoxList blc_mat_valid;
            for (int i=0; i<bac_interp.size(); ++i) {
              blc_mat_valid.join(BoxLib::complementIn(bac_interp[i],bl_mixed_local));
            }
            blc_mat_valid.simplify();
            bac_interp = (BoxArray)blc_mat_valid;
          }
          
          MultiFab mfc_interp(bac_interp,1,0);
          SetProperty(t,lev,mfc_interp,pname,0,0);

          // Piecewise constant "interpolation" of coarse values
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
          // mff_interp may be larger than remaining because it is generated by coarsen+refine
          BoxArray baf_interp_valid = BoxLib::intersect(baf_interp,BoxArray(remaining));
          MultiFab mff_interpt(baf_interp_valid,1,0);
          mff_interpt.copy(mff_interp); // Get only the part we want here
          mf.copy(mff_interpt,0,dComp,1); // Then get the part we want into the final container
          
          BoxList not_just_filled;
          BoxList blf_interp_valid(baf_interp_valid);
          for (int i=0; i<remaining.size(); ++i) {
            BoxList this_not_filled = BoxLib::complementIn(remaining[i],blf_interp_valid);
            not_just_filled.join(this_not_filled);
          }
          not_just_filled.simplify();
          remaining = BoxArray(not_just_filled);
        }
      }
    }
  }
}

void
MatFillerPC::FillCoarseCells(Real t, int level, MultiFab& mfc, const std::string& pname, int dComp)
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
      Average(mff[mfi],0,mfc[mfi],cbox,0,1,ref);
    }
  }
}

void
MatFillerPCarithAvg::Average(const FArrayBox& fineFab,
                             int              sComp,
                             FArrayBox&       crseFab,
                             const Box&       crseBox,
                             int              dComp,
                             int              nComp,
                             const IntVect&   ref) const
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

