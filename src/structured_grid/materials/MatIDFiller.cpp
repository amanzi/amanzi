#include <MatIDFiller.H>
#include <TagBox.H>
#include <Cluster.H>

static Real grid_eff = 1;
static int tags_buffer = 1;


MatIDFiller::Material::Material(const std::string&    _name, 
                                const PArray<Region>& _regions)
  : name(_name)
{
  regions.resize(_regions.size(),PArrayNoManage);
  for (int i=0; i<regions.size(); ++i) {
    regions.set(i,(Region*)&(_regions[i]));
  }
}

MatIDFiller::Material::Material(const MatIDFiller::Material& rhs) 
{
  name = rhs.name;
  if (rhs.regions.size()>0) {
    regions.resize(rhs.regions.size(),PArrayNoManage);
    for (int i=0; i<regions.size(); ++i) {
      regions.set(i,(Region*)&(rhs.regions[i]));
    }
  }
}


MatIDFiller::Material::~Material()
{}

void 
MatIDFiller::Material::setVal(FArrayBox& fab,Real val,int comp, const Real* dx) const
{
  int ng = 0;
  for (int k=0; k<regions.size(); ++k) {
    regions[k].setVal(fab,val,comp,dx,ng);
  }
}

MatIDFiller::MatIDFiller()
{
  initialized = false;
}

MatIDFiller::MatIDFiller(const Array<Geometry>& _geomArray,
                         const Array<IntVect>&  _refRatio,
                         const PArray<MatIDFiller::Material>& _materials)
{
  initialized = false;
  define(_geomArray,_refRatio,_materials);
}

void
MatIDFiller::define(const Array<Geometry>& _geomArray,
                    const Array<IntVect>&  _refRatio,
                    const PArray<MatIDFiller::Material>& _materials)
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
    materials.set(i,new MatIDFiller::Material(_materials[i]));
    matNames[cnt] = materials[i].Name();
    matIdx[materials[i].Name()] = cnt++;
  }
  Initialize();
}


void
MatIDFiller::SetMaterialID(int level, const Box& box, FArrayBox& fab)
{
  BoxArray unfilled(box);

  if (level<ba_mixed.size() && ba_mixed[level].size()>0) {
    std::vector< std::pair<int,Box> > isects = ba_mixed[level].intersections(box);
    for (int i=0; i<isects.size(); ++i) {
      fab.setVal(-1,isects[i].second,0,1); // Something invalid
    }
    unfilled = BoxLib::complementIn(box,ba_mixed[level]);
  }

  if (unfilled.ok()) {
    FArrayBox tfab;
    const Geometry& geom = Geom(level);
    const Real* dx = geom.CellSize();
    for (int i=0; i<unfilled.size(); ++i) {
      const Box& bx = unfilled[i];
      tfab.resize(bx,1); tfab.setVal(-1);
      for (int j=0; j<materials.size(); ++j) {
        int matID = matIdx[materials[j].Name()];
	materials[j].setVal(tfab,matID,0,dx);
      }
      fab.copy(tfab,bx,0,bx,0,1);
    }
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
  materialID.set(0,new MultiFab(BoxArray(Geom(0).Domain()), 1, 0));  
  for (int lev=0; lev<ba_mixed.size(); ++lev) {
    BoxArray fba = BoxArray(ba_mixed[lev]).refine(RefRatio(lev));
    fba.removeOverlap();
    materialID.set(lev+1,new MultiFab(fba, 1, 0));
  }
  for (int lev=0; lev<materialID.size(); ++lev) {
    for (MFIter mfi(materialID[lev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      FArrayBox& fab = materialID[lev][mfi];
      SetMaterialID(lev,box,fab);
    }
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
      Box gbox = BoxLib::grow(Geom(lev).Domain(),tags_buffer);
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
  const Real* dxFinest = Geom(nLevs-1).CellSize();
  const Box& box = Geom(0).Domain();
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

    Array<IntVect> tags(tbar[lev].numTags()); 
    long len = tbar[lev].collate(tags.dataPtr(), 0);
    ClusterList clist(tags.dataPtr(), len);
    clist.chop(grid_eff);
    BoxList bl = clist.boxList(); bl.simplify();
    ba_array[lev] = BoxLib::intersect(BoxArray(bl),Geom(lev).Domain());
  }

  return ba_array;
}

