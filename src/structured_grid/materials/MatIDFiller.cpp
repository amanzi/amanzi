#include <MatIDFiller.H>
#include <TagBox.H>
#include <Cluster.H>

static Real grid_eff = 1;
static int tags_buffer = 1;
static int max_grid_size = 16;

MatIDFiller::MatIDFiller()
{
  initialized = false;
}

MatIDFiller::MatIDFiller(const Array<Geometry>& _geomArray,
                         const Array<IntVect>&  _refRatio,
                         const PArray<Material>& _materials)
{
  initialized = false;
  define(_geomArray,_refRatio,_materials);
}

void
MatIDFiller::define(const Array<Geometry>& _geomArray,
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
MatIDFiller::SetMaterialID(int level, MultiFab& mf, int nGrow)
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
MatIDFiller::Mixed(int level) const 
{
  static BoxArray DUMMY;
  return (level >= ba_mixed.size() ? DUMMY : ba_mixed[level]);
}

void
MatIDFiller::Initialize()
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
MatIDFiller::FindMixedCells()
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

