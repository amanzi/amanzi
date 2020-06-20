#include <MatFiller.H>
#include <TagBox.H>
#include <Cluster.H>
#include <BC_TYPES.H>
#include <MatFiller_F.H>

static Real grid_eff = 1;
static int tags_buffer = 1;
static int max_grid_size = 32;

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

const Geometry&
MatFiller::Geom(int level) const
{
  if (! Initialized() ) {
    BoxLib::Abort("MatFiller not initialized");
  }
  return geomArray[level];
}

const IntVect&
MatFiller::RefRatio(int crse_level) const
{
  return refRatio[crse_level];
}

iMultiFab&
MatFiller::MaterialID(int level)
{
  if (! Initialized() ) {
    BoxLib::Abort("MatFiller not initialized");
  }
  return materialID[level];
}

const std::map<std::string,int>&
MatFiller::MatIdx() const
{
  if (! Initialized() ) {
    BoxLib::Abort("MatFiller not initialized");
  }
  return matIdx;
}

int
MatFiller::NumLevels() const
{
  if (! Initialized() ) {
    BoxLib::Abort("MatFiller not initialized");
  }
  return nLevs;
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
MatFiller::SetMaterialID(int level, iMultiFab& mf, int nGrow, bool ignore_mixed) const
{
  BoxArray unfilled(mf.boxArray());
  if (unfilled.size()==0) {
    return;
  }

  unfilled.grow(nGrow);
  iMultiFab tmf(unfilled,1,0);

  const Geometry& geom = geomArray[level];
  const Real* dx = geom.CellSize();
  IArrayBox tfab;
  for (MFIter mfi(tmf); mfi.isValid(); ++mfi) {
    IArrayBox& tfab = tmf[mfi];
    for (int j=0; j<materials.size(); ++j) {
      int matID = -1;
      std::map<std::string,int>::const_iterator it = matIdx.find(materials[j].Name());
      if (it!=matIdx.end()) {
        matID = it->second;
      }
      BL_ASSERT(matID >= 0);
      materials[j].setVal(tfab,matID,0,dx);
    }
  }

  if (!ignore_mixed && level<ba_mixed.size() && ba_mixed[level].size()>0) {
    BoxArray mixed = BoxLib::intersect(ba_mixed[level], unfilled);
    if (mixed.size()>0) {
      iMultiFab mmf(mixed,1,0);
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
  BoxArray cba = BoxArray(geomArray[0].Domain());
  materialID.set(0,new iMultiFab(cba, 1, nGrow));  
  for (int lev=0; lev<ba_mixed.size(); ++lev) {
    BoxList fbl(ba_mixed[lev]); fbl.refine(RefRatio(lev));
    fbl.simplify(); fbl.maxSize(max_grid_size);
    BoxArray fba(fbl);
    BL_ASSERT(fba.isDisjoint());
    if (fba.size()!=0) { 
      materialID.set(lev+1,new iMultiFab(fba, 1, nGrow));
    }
    else {
      materialID.set(lev+1,new iMultiFab);
    }
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

  Array<IntVect> cr(nLevs);
  cr[nLevs-1] = IntVect(D_DECL(1,1,1));
  for (int lev=nLevs-2; lev>=0; --lev) {
    cr[lev] = cr[lev+1] * RefRatio(lev);
  }

  int nMat = materials.size();
  Array<BoxArray> material_bounds(nMat);
  for (int j=0; j<nMat; ++j) {
    UnionRegion union_region("test","all",materials[j].Regions());
    material_bounds[j] = union_region.approximate_bounds(Array<Real>(geomArray[nLevs-1].ProbLo(),BL_SPACEDIM),
                                                         Array<Real>(geomArray[nLevs-1].CellSize(),BL_SPACEDIM));
  }

  BoxArray cba(geomArray[0].Domain());
  //int init_mat_granularity = std::max(1.,(Real)geomArray[0].Domain().numPts()/(10*ParallelDescriptor::NProcs()));
  int init_mat_granularity = 1;
  cba.maxSize(init_mat_granularity);
  int cba_size = cba.size();
  MultiFab junk(cba,1,0,Fab_noallocate);
  const DistributionMapping& dm = junk.DistributionMap();
  int my_proc = ParallelDescriptor::MyProc();
  bool is_ioproc = ParallelDescriptor::IOProcessor();
  Array<int> my_grids;
  for (int i=0; i<cba_size ; ++i) {
    if (dm[i] == my_proc) {
      my_grids.push_back(i);
    }
  }

  // Find all boxes in cba that possibly contain more than 1 material
  // If only 1 material, store material idx 
  std::vector< std::pair<int,Box> > isects;
  bool first_only = true;
  Array<int> grid_to_mat(cba_size,-1); // A convenient way to share results over procs later

  for (int i=0; i<ParallelDescriptor::NProcs(); ++i) {
    if (my_proc==i) {
      for (int j=0; j<my_grids.size(); ++j) {
        const Box& cbox = cba[my_grids[j]];
        Box fbox = Box(cbox).refine(cr[0]);

        Array<int> hits;
        for (int k=0; k<nMat; ++k) {
          bool bounds_provided = material_bounds[k].size() > 0;
          if (bounds_provided) {
            material_bounds[k].intersections(fbox,isects,first_only,0);
          }
          if (!bounds_provided || isects.size()>0) {
            hits.push_back(k);
          }
        }
        int nhits = hits.size();
        if (nhits == 0) {
          BoxLib::Abort("Material not defined for part of domain");
        }
        grid_to_mat[my_grids[j]] = (int)( hits.size() > 1  ?  -1  :  hits[0]);
      }
    }
  }
  ParallelDescriptor::ReduceIntSum(grid_to_mat.dataPtr(),cba_size); // Share results

  // Build reduced cba, containing only the interesting boxes
  BoxList cbl_red;
  for (int i=0; i<grid_to_mat.size(); ++i) {
    if (grid_to_mat[i] < 0) {
      cbl_red.push_back(cba[i]);
    }
  }
  BoxArray cba_red(cbl_red);
  MultiFab junk_red(cba_red,1,0,Fab_noallocate);
  const DistributionMapping& dm_red = junk_red.DistributionMap();
  int num_red = cba_red.size();

  if (num_red) {
    PArray<TagBoxArray> tbar;
    if (nLevs>1) {
      tbar.resize(nLevs-1, PArrayManage);
      for (int lev=nLevs-2; lev>=0; --lev) {
	BoxArray gba = BoxArray(cba_red).refine(cr[0]).coarsen(cr[lev]);
	tbar.set(lev, new TagBoxArray(gba,tags_buffer)); 
      }
    }

    IArrayBox finestFab, coarseFab;
    const Real* dxFinest = geomArray[nLevs-1].CellSize();
    for (int i=0; i<num_red; ++i) {
      if (dm_red[i] == my_proc) {
	const Box& box = cba_red[i];
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
	    TagBox& tb = tbar[lev][i];
          
	    for (IntVect thisIv=thisBox.smallEnd(), TEnd=thisBox.bigEnd(); thisIv<=TEnd; thisBox.next(thisIv)) {      
	      IntVect thisFineIv = thisIv * cr[lev];
	      Box thisFineBox(thisFineIv,thisFineIv+rm1);      
	      if (finestFab.min(thisFineBox,0) != finestFab.max(thisFineBox,0)) {
		tb.setVal(TagBox::SET,Box(thisIv,thisIv),0,1);
	      }
	    }
	  }
	}
      }
    }

    for (int lev=nLevs-2; lev>=0; --lev) {
      if (lev < nLevs-2) {
	tbar[lev+1].coarsen(RefRatio(lev));
	for (int i=0; i<num_red; ++i) {
	  if (dm_red[i] == my_proc) {
	    tbar[lev][i].merge(tbar[lev+1][i]);
	  }
	}
      }
      tbar[lev].buffer(tags_buffer);

      std::vector<IntVect> tags;
      tbar[lev].collate(tags);
      long int num_tags = tags.size();
      if (num_tags>0) {
	ClusterList clist(&(tags[0]), num_tags);
	clist.chop(grid_eff);
	BoxList bl = clist.boxList(); bl.simplify();
	ba_array[lev] = BoxLib::intersect(BoxArray(bl),geomArray[lev].Domain());
	ba_array[lev].maxSize(max_grid_size);
      }
    }
  }
  return ba_array;
}

bool
MatFiller::CanDerive(const std::string& property_name) const
{
  return property_nComps.find(property_name) != property_nComps.end();
}

int
MatFiller::NComp(const std::string& property_name) const
{
  std::map<std::string,int>::const_iterator it=property_nComps.find(property_name);
  return (it==property_nComps.end() ? 0 : it->second);
}

bool 
MatFiller::SetProperty(Real               t,
                       int                level,
                       MultiFab&          mf,
                       const std::string& pname,
                       int                dComp,
                       int                nGrow,
                       void*              ctx,
                       bool               ignore_mixed) const
{
  if (geomArray[level].isAnyPeriodic()) {
    BoxLib::Abort("Periodic not yet supported");
  }

  BoxArray unfilled(mf.boxArray()); 
  if (unfilled.size()==0) {
    return true;
  }
  unfilled.grow(nGrow);

  // Find number of comps for this property, but also use this to decide
  // if we know about this property at all.
  std::map<std::string,int>::const_iterator it=property_nComps.find(pname);
  if (it==property_nComps.end()) {
    return false;
  }
  int nComp = it->second;

  BL_ASSERT(mf.nComp() >= dComp + nComp);
  MultiFab tmf(unfilled,nComp,0);

  Array<Real> tmpv;
  Array<Array<Real> > values(nComp, Array<Real>(materials.size()));
  for (int i=0; i<materials.size(); ++i) {
    const Property* p = materials[i].Prop(pname);
    if (p==0) {
      return false;
    }
    p->Evaluate(t,tmpv);
    for (int n=0; n<nComp; ++n) {
      values[n][i] = tmpv[n];
    }
  }

  BoxArray baM;
  if (level<ba_mixed.size() && ba_mixed[level].size()>0) {
    baM = BoxLib::intersect(ba_mixed[level],unfilled); 
    if (baM.size()>0) {
      baM.removeOverlap();
      MultiFab mixed(baM,nComp,0);
      if (ignore_mixed) {
        iMultiFab id(baM,nComp,0);
        SetMaterialID(level,id,0,ignore_mixed);
        for (MFIter mfi(id); mfi.isValid(); ++mfi) {
          const IArrayBox& idfab = id[mfi];
          FArrayBox& mfab = mixed[mfi];
          const Box& bx = mfi.validbox();
          for (int n=0; n<nComp; ++n) {
            FORT_FILLP (mfab.dataPtr(n), ARLIM(mfab.loVect()),  ARLIM(mfab.hiVect()),
                        idfab.dataPtr(), ARLIM(idfab.loVect()), ARLIM(idfab.hiVect()),
                        bx.loVect(), bx.hiVect(), values[n].dataPtr());
          }
        }
      } else {
        bool ret = FillCoarseCells(t,level,mixed,pname,0,nComp,ctx);
        if (!ret) return false;
      }
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
      BoxArray ba_fillable(bl_fillable); ba_fillable.removeOverlap();
      MultiFab fillData(ba_fillable,nComp,0);
      iMultiFab fillID(ba_fillable,1,0);
      fillID.copy(materialID[level]); // guaranteed to be filled completely
      for (MFIter mfi(fillData); mfi.isValid(); ++mfi) {
	const Box& ovlp = mfi.validbox();
	const IArrayBox& idfab = fillID[mfi];
	FArrayBox& matfab = fillData[mfi];
	for (int n=0; n<nComp; ++n) {
	  FORT_FILLP (matfab.dataPtr(n), ARLIM(matfab.loVect()), ARLIM(matfab.hiVect()),
		      idfab.dataPtr(),   ARLIM(idfab.loVect()),  ARLIM(idfab.hiVect()),
		      ovlp.loVect(), ovlp.hiVect(), values[n].dataPtr());
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
            bool ret = SetProperty(t,lev,mfc_interp,pname,0,0,ctx);
            if (!ret) return false;
            
            BoxArray baf_interp(bac_interp); baf_interp.refine(cumRatio);
            MultiFab mff_interp(baf_interp,nComp,0);
            for (MFIter mfi(mfc_interp); mfi.isValid(); ++mfi) {
              const Box& cbox = mfi.validbox();
              const FArrayBox& cfab = mfc_interp[mfi];
              FArrayBox& ffab = mff_interp[mfi];
              RefineData(cfab,0,ffab,cbox,0,nComp,cumRatio,refineRule(pname));
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

  FillCellsOutsideDomain(t,level,mf,dComp,nComp,geomArray[level]);
  if (nGrow>0) {
    mf.FillBoundary(dComp,nComp);
  }

  return true;
}

void
MatFiller::FillCellsOutsideDomain(Real               t,
                                  int                level,
                                  MultiFab&          mf,
                                  int                dComp,
                                  int                nComp,
				  const Geometry&    geom)
{
  const Array<int> bc(2*BL_SPACEDIM,FOEXTRAP);
  const Box& domain = geom.Domain();
  const Real* dx = geom.CellSize();
  const Real* plo = geom.ProbLo();
  for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
    FArrayBox& fab = mf[mfi];

    Box bx = fab.box() & domain;
    if (bx.ok()) {
      for (int n=0; n<nComp; ++n) {
	FORT_FILCC(fab.dataPtr(dComp+n), ARLIM(fab.loVect()), ARLIM(fab.hiVect()),
		   domain.loVect(), domain.hiVect(),dx,plo,bc.dataPtr());
      }
    }
    else {
      std::cout << "no part of box inside domain: " << fab.box() << std::endl;
      std::cout << "domain: " << domain << std::endl;
      BoxLib::Abort();
    }
  }
}

bool
MatFiller::FillCoarseCells(Real               t,
                           int                level,
                           MultiFab&          mfc,
                           const std::string& pname,
                           int                dComp,
                           int                nComp,
                           void*              ctx) const
{
  int finestLevel = nLevs-1;
  if (level<finestLevel) {
    const BoxArray& bac = mfc.boxArray();
    const IntVect& ref = RefRatio(level);
    BoxArray baf = BoxArray(bac).refine(ref);
    MultiFab mff(baf,nComp,0);
    bool ret = SetProperty(t,level+1,mff,pname,0,0,ctx);
    if (!ret) return false;
    for (MFIter mfi(mfc); mfi.isValid(); ++mfi) {
      const Box& cbox = mfi.validbox();
      CoarsenData(mff[mfi],0,mfc[mfi],cbox,dComp,nComp,ref,coarsenRule(pname));
    }
  }
  return true;
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
                       const Property::CoarsenRule& rule)
{
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
                      const Property::RefineRule& rule)
{
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

  if (rule==Property::PiecewiseConstant) {
    FORT_REFINEPC(crseFab.dataPtr(sComp),ARLIM(cbox.loVect()),ARLIM(cbox.hiVect()),
                  fineFab.dataPtr(dComp),ARLIM(fbox.loVect()),ARLIM(fbox.hiVect()),
                  cwbox.loVect(),cwbox.hiVect(),rvect,&nComp);
  }
  else {
    BoxLib::Abort("No appropriate refining implementation");
  }
}
