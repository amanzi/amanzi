#include <MatFiller.H>
#include <TagBox.H>
#include <Cluster.H>
#include <BC_TYPES.H>
#include <MatFiller_F.H>

static Real grid_eff = 1;
static int tags_buffer = 1;
static int max_grid_size = 16;


MatFiller::MatFiller(const Array<Geometry>&  _geomArray,
                     const Array<IntVect>&   _refRatio,
                     const PArray<Material>& _materials)
{
  initialized = false;
  define(_geomArray,_refRatio,_materials);
  VerifyProperties();
}

void
MatFiller::define(const Array<Geometry>& _geomArray,
		  const Array<IntVect>&  _refRatio,
		  const PArray<Material>& _materials)
{
  if (initialized) {
    geomArray.clear();
    refRatio.clear();
    nLevs = 0;
    materials.clear();
    matIdx.clear();
    matNames.clear();
    initialized = false;
  }

  geomArray = _geomArray;
  refRatio = _refRatio;
  nLevs = geomArray.size();
  BL_ASSERT(nLevs == refRatio.size()+1);
  
  int cnt = 0;
  materials.resize(_materials.size(),PArrayManage);
  matNames.resize(materials.size());
  for (int i=0; i<materials.size(); ++i) {
    materials.set(i,new Material(_materials[i]));
    matNames[cnt] = materials[i].Name();
    matIdx[materials[i].Name()] = cnt++;
  }
  Initialize();
}

void
MatFiller::VerifyProperties()
{
  // Get the number of components for each property, ensure that all materials will return the
  //  the same number of components
  int nMat = materials.size();
  if (nMat>0) {
    const Array<std::string> property_names = materials[0].PropertyNames();
    for (int j=0; j<property_names.size(); ++j) {
      int nComp = -1;
      Property::CoarsenRule cRule = Property::INVALID_CR;
      Property::RefineRule rRule = Property::INVALID_RR;
      const std::string& pname = property_names[j];
      for (int i=0; i<materials.size(); ++i) {
        const Property* p = materials[i].Prop(pname);
        BL_ASSERT(p!=0);
        int this_nComp = p->nComp();
        Property::CoarsenRule this_cRule = p->coarsenRule();
        Property::RefineRule this_rRule = p->refineRule();
        if (i==0) {
          nComp = this_nComp;
          cRule = this_cRule;
          rRule = this_rRule;
        }
        else {
          BL_ASSERT(nComp == this_nComp);
          BL_ASSERT(cRule == this_cRule);
          BL_ASSERT(rRule == this_rRule);
        }
      }
      property_nComps[pname] = nComp;
      property_cRules[pname] = cRule;
      property_rRules[pname] = rRule;
    }
  }
}

void 
MatFiller::SetMaterialID(int level, MultiFab& mf, int nGrow)
{
  BoxArray unfilled(mf.boxArray());
  if (unfilled.size()==0) {
    return;
  }

  unfilled.grow(nGrow);
  MultiFab tmf(unfilled,1,0);

  BL_ASSERT(level<NumLevels());
  const Geometry& geom = geomArray[level];
  const Real* dx = geom.CellSize();
  FArrayBox tfab;
  for (MFIter mfi(tmf); mfi.isValid(); ++mfi) {
    FArrayBox& tfab = tmf[mfi];
    for (int j=0; j<materials.size(); ++j) {
      int matID = matIdx[materials[j].Name()];
      materials[j].setVal(tfab,matID,0,dx);
    }
  }
  
  if (level<ba_mixed.size() && ba_mixed[level].size()>0) {
    BoxArray mixed = BoxLib::intersect(ba_mixed[level], unfilled);
    if (mixed.size()>0) {
      MultiFab mmf(mixed,1,0);
      mmf.setVal(-1); // Something invalid
      tmf.copy(mmf);
    }
  }

  for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
    mf[mfi].copy(tmf[mfi],0,0,1);
  }
}

const BoxArray& 
MatFiller::Mixed(int level) const 
{
  static BoxArray DUMMY;
  return (level >= ba_mixed.size() ? DUMMY : ba_mixed[level]);
}

void
MatFiller::Initialize()
{
  ba_mixed = FindMixedCells();

  int finestLevel = ba_mixed.size();
  materialID.resize(finestLevel+1,PArrayManage);
  int nGrow = 0;
  materialID.set(0,new MultiFab(BoxArray(geomArray[0].Domain()), 1, nGrow));  
  for (int lev=0; lev<ba_mixed.size(); ++lev) {
    BoxList fbl(ba_mixed[lev]); fbl.refine(RefRatio(lev));
    fbl.simplify(); fbl.maxSize(max_grid_size);
    BoxArray fba(fbl);
    BL_ASSERT(fba.isDisjoint());
    materialID.set(lev+1,new MultiFab(fba, 1, nGrow));
  }
  for (int lev=0; lev<materialID.size(); ++lev) {
    SetMaterialID(lev,materialID[lev],0);
  }
  initialized = true;
}

Array<BoxArray>
MatFiller::FindMixedCells()
{
  Array<BoxArray> ba_array(nLevs-1);
  PArray<TagBox> tbar;
  if (nLevs>1) {
    tbar.resize(nLevs-1, PArrayManage);
    for (int lev=nLevs-2; lev>=0; --lev) {
      Box gbox = BoxLib::grow(geomArray[lev].Domain(),tags_buffer);
      tbar.set(lev, new TagBox(gbox)); // Workaround for bug in buffer code
    }
  }

  Array<IntVect> cr(nLevs);
  cr[nLevs-1] = IntVect(D_DECL(1,1,1));
  for (int lev=nLevs-2; lev>=0; --lev) {
    cr[lev] = cr[lev+1] * RefRatio(lev);
  }

  Box cbox;
  FArrayBox finestFab, coarseFab;
  const Real* dxFinest = geomArray[nLevs-1].CellSize();
  const Box& box = geomArray[0].Domain();
  for (IntVect civ=box.smallEnd(), CEnd=box.bigEnd(); civ<=CEnd; box.next(civ)) {
    Box coarsestBox(civ,civ);
    Box finestBox = Box(coarsestBox).refine(cr[0]);
    finestFab.resize(finestBox,1); finestFab.setVal(-1);
    for (int j=0; j<materials.size(); ++j) {
      int matID = matIdx[materials[j].Name()];
      materials[j].setVal(finestFab,matID,0,dxFinest);
    }

    for (int lev=nLevs-2; lev>=0; --lev) {
      Box thisBox = Box(finestBox).coarsen(cr[lev]);
      IntVect rm1 = cr[lev] - IntVect::TheUnitVector();
      TagBox& tb = tbar[lev];

      for (IntVect thisIv=thisBox.smallEnd(), TEnd=thisBox.bigEnd(); thisIv<=TEnd; thisBox.next(thisIv)) {      
        IntVect thisFineIv = thisIv * cr[lev];
        Box thisFineBox(thisFineIv,thisFineIv+rm1);      
        if (finestFab.min(thisFineBox,0) != finestFab.max(thisFineBox,0)) {
          tb.setVal(TagBox::SET,Box(thisIv,thisIv),0,1);
        }
      }
    }
  }


  for (int lev=nLevs-2; lev>=0; --lev) {
    if (lev < nLevs-2) {
      TagBox* fromFine = tbar[lev+1].coarsen(RefRatio(lev));
      fromFine->buffer(tags_buffer,tags_buffer);
      tbar[lev].merge(*fromFine);
      delete fromFine;
    }

    int num_tags = tbar[lev].numTags();
    if (num_tags>0) {
      Array<IntVect> tags(num_tags);
      long len = tbar[lev].collate(tags.dataPtr(), 0);
      ClusterList clist(tags.dataPtr(), len);
      clist.chop(grid_eff);
      BoxList bl = clist.boxList(); bl.simplify();
      ba_array[lev] = BoxLib::intersect(BoxArray(bl),geomArray[lev].Domain());
      ba_array[lev].maxSize(max_grid_size);
    }
  }

  return ba_array;
}

void 
MatFiller::SetProperty(Real               t,
                       int                level,
                       MultiFab&          mf,
                       const std::string& pname,
                       int                dComp,
                       int                nGrow,
                       void*              ctx)
{
  if (geomArray[level].isAnyPeriodic()) {
    BoxLib::Abort("Periodic not yet supported");
  }

  BoxArray unfilled(mf.boxArray()); 
  if (unfilled.size()==0) {
    return;
  }
  unfilled.grow(nGrow);

  // Make a handy data structure
  std::vector<const Property*> props(materials.size());
  for (int i=0; i<materials.size(); ++i) {
    const Property* p = materials[i].Prop(pname);
    BL_ASSERT(p!=0);
    props[i] = p;
  }

  std::map<std::string,int>::const_iterator it=property_nComps.find(pname);
  BL_ASSERT(it!=property_nComps.end());
  int nComp = it->second;
  BL_ASSERT(mf.nComp() >= dComp + nComp);
  MultiFab tmf(unfilled,nComp,0);
  
  BoxArray baM;
  if (level<ba_mixed.size() && ba_mixed[level].size()>0) {
    baM = BoxLib::intersect(ba_mixed[level],unfilled); 
    if (baM.size()>0) {
      baM.removeOverlap();
      MultiFab mixed(baM,nComp,0);
      FillCoarseCells(t,level,mixed,pname,0,nComp,ctx);
      tmf.copy(mixed);
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

    // Find portion of unfilled that are fillable by non-mixed cells in the materialID struct at this level
    BoxList bl_fillable = BoxLib::intersect((BoxList)materialID[level].boxArray(),(BoxList)unfilled);
    BL_ASSERT( BoxLib::intersect((BoxArray)bl_fillable,baM).size()==0 );

    if (bl_fillable.size()>0) {
      bl_fillable.simplify();
      BoxArray ba_fillable(bl_fillable);
      MultiFab fillData(ba_fillable,nComp,0);
      MultiFab fillID(ba_fillable,1,0); 
      fillID.copy(materialID[level]); // guaranteed to be filled completely
      for (MFIter mfi(fillData); mfi.isValid(); ++mfi) {
	const Box& ovlp = mfi.validbox();
	const FArrayBox& idfab = fillID[mfi];
	FArrayBox& matfab = fillData[mfi];
	for (IntVect iv=ovlp.smallEnd(), End=ovlp.bigEnd(); iv<=End; ovlp.next(iv)) {
	  props[idfab(iv,0)]->eval(t,level,Box(iv,iv),matfab,0,ctx);
	}
      }
      tmf.copy(fillData);
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

          if (bac_interp.size()>0) {
            MultiFab mfc_interp(bac_interp,nComp,0);
            SetProperty(t,lev,mfc_interp,pname,0,0,ctx);
            
            BoxArray baf_interp(bac_interp); baf_interp.refine(cumRatio);
            MultiFab mff_interp(baf_interp,nComp,0);
            for (MFIter mfi(mfc_interp); mfi.isValid(); ++mfi) {
              const Box& cbox = mfi.validbox();
              const FArrayBox& cfab = mfc_interp[mfi];
              FArrayBox& ffab = mff_interp[mfi];
              RefineData(cfab,0,ffab,cbox,0,nComp,cumRatio,pname);
            }
            // mff_interp may be larger than remaining because it is generated by coarsen+refine
            BoxArray baf_interp_valid = BoxLib::intersect(baf_interp,BoxArray(remaining));
            MultiFab mff_interpt(baf_interp_valid,nComp,0);
            mff_interpt.copy(mff_interp); // Get only the part we want here
            tmf.copy(mff_interpt,0,0,nComp); // Then get the part we want 
            
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

  for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
    mf[mfi].copy(tmf[mfi],0,dComp,nComp);
  }

  FillCellsOutsideDomain(t,level,mf,pname,dComp,nComp,nGrow);
}

void
MatFiller::FillCellsOutsideDomain(Real               t,
                                  int                level,
                                  MultiFab&          mf,
                                  const std::string& pname,
                                  int                dComp,
                                  int                nComp,
                                  int                nGrow)
{
  const Array<int> bc(2*BL_SPACEDIM,FOEXTRAP);
  const Geometry& geom = geomArray[level];
  const Box& domain = geom.Domain();
  const Real* dx = geom.CellSize();
  for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
    FArrayBox& fab = mf[mfi];
    for (int n=0; n<nComp; ++n) {
      FORT_FILCC(fab.dataPtr(dComp+n), ARLIM(fab.loVect()), ARLIM(fab.hiVect()),
                 domain.loVect(), domain.hiVect(),dx,dx,bc.dataPtr());
    }
  }
}

void
MatFiller::FillCoarseCells(Real               t,
                           int                level,
                           MultiFab&          mfc,
                           const std::string& pname,
                           int                dComp,
                           int                nComp,
                           void*              ctx)
{
  int finestLevel = nLevs-1;
  if (level<finestLevel) {
    const BoxArray& bac = mfc.boxArray();
    const IntVect& ref = RefRatio(level);
    BoxArray baf = BoxArray(bac).refine(ref);
    MultiFab mff(baf,nComp,0);
    SetProperty(t,level+1,mff,pname,0,0,ctx);
    for (MFIter mfi(mfc); mfi.isValid(); ++mfi) {
      const Box& cbox = mfi.validbox();
      CoarsenData(mff[mfi],0,mfc[mfi],cbox,dComp,nComp,ref,pname);
    }
  }
}

Property::CoarsenRule
MatFiller::coarsenRule(const std::string& pname) const
{
  std::map<std::string,Property::CoarsenRule>::const_iterator it = property_cRules.find(pname);
  return ( it == property_cRules.end() ? Property::INVALID_CR : it->second );
}

Property::RefineRule
MatFiller::refineRule(const std::string& pname) const
{
  std::map<std::string,Property::RefineRule>::const_iterator it = property_rRules.find(pname);
  return ( it == property_rRules.end() ? Property::INVALID_RR : it->second );
}

int
MatFiller::nComp(const std::string& property_name) const
{
  std::map<std::string,int>::const_iterator it = property_nComps.find(property_name);
  return ( it == property_nComps.end() ? -1 : it->second );
}

void
MatFiller::CoarsenData(const FArrayBox&   fineFab,
                       int                sComp,
                       FArrayBox&         crseFab,
                       const Box&         crseBox,
                       int                dComp,
                       int                nComp,
                       const IntVect&     ref,
                       const std::string& pname) const
{
  IntVect refm = IntVect(ref)-IntVect(D_DECL(1,1,1));
  Box lbox(IntVect(D_DECL(0,0,0)),refm);

  // Shift to +ve indices so that coarsening stuff works correctly
  IntVect cshift;
  for (int d=0; d<BL_SPACEDIM; ++d) {
    cshift[d] = std::max(0,-crseBox.smallEnd()[d]);
  }
  IntVect fshift = IntVect(cshift) * ref;
  Box cbox = crseFab.box(); cbox.shift(cshift);
  Box fbox = fineFab.box(); fbox.shift(fshift);
  Box cwbox = crseBox;      cwbox.shift(cshift);
  const int* rvect = ref.getVect();
  BL_ASSERT(fbox.contains(Box(cbox).refine(ref)));

  Property::CoarsenRule rule = coarsenRule(pname);
  if (rule==Property::Arithmetic) {
    FORT_CRSNARITH(fineFab.dataPtr(sComp),ARLIM(fbox.loVect()),ARLIM(fbox.hiVect()),
                   crseFab.dataPtr(dComp),ARLIM(cbox.loVect()),ARLIM(cbox.hiVect()),
                   cwbox.loVect(),cwbox.hiVect(),rvect,&nComp);
  }
  else if (rule == Property::ComponentHarmonic) {
    FORT_CRSNHARM(fineFab.dataPtr(sComp),ARLIM(fbox.loVect()),ARLIM(fbox.hiVect()),
                  crseFab.dataPtr(dComp),ARLIM(cbox.loVect()),ARLIM(cbox.hiVect()),
                  cwbox.loVect(),cwbox.hiVect(),rvect,&nComp);
  }
  else {
    BoxLib::Abort("No appropriate coarsening implementation");
  }
}

void
MatFiller::RefineData(const FArrayBox&   crseFab,
                      int                sComp,
                      FArrayBox&         fineFab,
                      const Box&         crseBox,
                      int                dComp,
                      int                nComp,
                      const IntVect&     ref,
                      const std::string& pname) const
{
  IntVect refm = IntVect(ref)-IntVect(D_DECL(1,1,1));
  Box lbox(IntVect(D_DECL(0,0,0)),refm);

  // Shift to +ve indices so that coarsening stuff works correctly
  IntVect cshift;
  for (int d=0; d<BL_SPACEDIM; ++d) {
    cshift[d] = std::max(0,-crseBox.smallEnd()[d]);
  }
  IntVect fshift = IntVect(cshift) * ref;
  Box cbox = crseFab.box(); cbox.shift(cshift);
  Box fbox = fineFab.box(); fbox.shift(fshift);
  Box cwbox = crseBox;      cwbox.shift(cshift);
  const int* rvect = ref.getVect();
  BL_ASSERT(fbox.contains(Box(cbox).refine(ref)));

  Property::RefineRule rule = refineRule(pname);
  if (rule==Property::PiecewiseConstant) {
    FORT_REFINEPC(crseFab.dataPtr(sComp),ARLIM(cbox.loVect()),ARLIM(cbox.hiVect()),
                  fineFab.dataPtr(dComp),ARLIM(fbox.loVect()),ARLIM(fbox.hiVect()),
                  cwbox.loVect(),cwbox.hiVect(),rvect,&nComp);
  }
  else {
    BoxLib::Abort("No appropriate refining implementation");
  }
}
