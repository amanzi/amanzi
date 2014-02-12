#include <iostream>
#include <fstream>
#include <iomanip>
using std::cout;
using std::endl;
#include <ParmParse.H>

#include <Material.H>
#include <MatFiller.H>
#include <RegionManager.H>

static std::map<std::string, Region*> regions;

Array<const Region*>
RegionPtrArray(const Array<std::string>& region_names)
{
  Array<const Region*> ret(region_names.size());

  for (int i=0; i<region_names.size(); ++i)
  {
    const std::string& name = region_names[i];
    std::map<std::string, Region*>::const_iterator it = regions.find(name);
    if (it != regions.end()) {
      ret[i] = it->second;
    }
    else {
      std::string m = "Named region not found: " + name;
      BoxLib::Error(m.c_str());
    }
  }
  return ret;
}

PArray<Material>
SetMaterials(RegionManager& rm)
{
  PArray<Material> materials(5, PArrayManage);
  Array<std::string> region_names;
  Array<const Region*> regionset;
  materials.set(0,new Material("Soil",rm.RegionPtrArray(region_names)));
  region_names.clear();
  region_names.push_back("TankConcFloor");
  region_names.push_back("TankConcRoof1");
  region_names.push_back("TankConcRoof2");
  region_names.push_back("TankConcWall");
  materials.set(1,new Material("TankConc",rm.RegionPtrArray(region_names)));
  region_names.clear();
  region_names.push_back("TankFFfloor");
  region_names.push_back("TankFFwall");
  region_names.push_back("TankWaste");
  materials.set(2,new Material("TankFF",rm.RegionPtrArray(region_names)));
  region_names.clear();
  region_names.push_back("TankGrout");
  materials.set(3,new Material("TankGrout",rm.RegionPtrArray(region_names)));
  region_names.clear();
  region_names.push_back("TankLinerFloor");
  region_names.push_back("TankLinerRoof");
  region_names.push_back("TankLinerWall");
  materials.set(4,new Material("TankLiner",rm.RegionPtrArray(region_names)));
  region_names.clear();
  return materials;
}

void
DestroyRegions()
{
  for (std::map<std::string, Region*>::iterator it=regions.begin(), End=regions.end(); it!=End; ++it) {
    delete it->second;
  }
}

void WritePlotfile(const std::string         &pfversion,
                   const PArray<MultiFab>    &data,
                   const Real                 time,
                   const Real                *probLo,
                   const Real                *probHi,
                   const Array<int>          &refRatio,
                   const Array<Box>          &probDomain,
                   const Array<Array<Real> > &dxLevel,
                   const int                  coordSys,
                   const std::string         &oFile,
                   const Array<std::string>  &names,
                   const bool                 verbose,
		   const bool                 isCartGrid,
		   const Real                *vfeps,
		   const int                 *levelSteps);

int
main (int   argc,
      char* argv[])
{
  BoxLib::Initialize(argc,argv);

  ParmParse pp;

  int nLevs = 3; pp.query("nLevs",nLevs);
  BL_ASSERT(nLevs>0);
  Array<int> n_cells;
  pp.getarr("n_cells",n_cells,0,BL_SPACEDIM);

  Array<int> rRatio(nLevs-1,4);
  if (nLevs>1) {
    pp.getarr("refine_ratio",rRatio,0,nLevs-1);
  }

  Array<IntVect> refRatio(rRatio.size());
  for (int lev=0; lev<rRatio.size(); ++lev) {
    refRatio[lev] = rRatio[lev] * IntVect::TheUnitVector();
  }

  Array<Geometry> geomArray(nLevs);
  for (int lev=0; lev<nLevs; ++lev) {
    Box domain;
    if (lev==0) {
      IntVect be;
      for (int d=0; d<BL_SPACEDIM; ++d) {
        be[d] = n_cells[d] - 1;
      }
      domain = Box(IntVect(D_DECL(0,0,0)),be);
    }
    else {
      domain = Box(geomArray[lev-1].Domain()).refine(refRatio[lev-1]);
    }
    geomArray[lev] = Geometry(domain);
  }

  Region::domlo.resize(BL_SPACEDIM);
  Region::domhi.resize(BL_SPACEDIM);
  for (int d=0; d<BL_SPACEDIM; ++d) {
    Region::domlo[d] = Geometry::ProbLo()[d];
    Region::domhi[d] = Geometry::ProbHi()[d];
  }

  RegionManager rm;
  PArray<Material> materials = SetMaterials(rm);

  MatFiller matFiller(geomArray,refRatio,materials);

  bool fail = false;

  const std::map<std::string,int>& mat_map = matFiller.MatIdx();

  Array<int> bins(mat_map.size(),0);

  for (int lev=0; lev<nLevs; ++lev) {
    const MultiFab& mf = matFiller.MaterialID(lev);
    if (lev<nLevs) {
      const BoxArray& ba_mixed = matFiller.Mixed(lev);
      if (ba_mixed.size()>0) {
        Real maxVal = -1;
        for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
          const FArrayBox& fab = mf[mfi];
          std::vector< std::pair<int,Box> > isects = ba_mixed.intersections(mfi.validbox());
          for (int ii = 0, N = isects.size(); ii < N; ii++)
          {
            maxVal = std::max(maxVal,fab.max(isects[ii].second,0));
          }
        }        
        ParallelDescriptor::ReduceRealMax(maxVal);
        fail = (maxVal>-1);
      }
    }

    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
      const Box& vbox = mfi.validbox();
      const FArrayBox& fab = mf[mfi];
      for (IntVect iv=vbox.smallEnd(), BIG=vbox.bigEnd(); iv<=BIG; vbox.next(iv)) {
        int val = (int) fab(iv,0);
        if (val>0) {
          bins[val]++;
        }
      }    
    }
  }

  ParallelDescriptor::ReduceIntSum(bins.dataPtr(),bins.size());

  fail &= bins[0] == 0;
  fail &= bins[1] == 610;
  fail &= bins[2] == 10;
  fail &= bins[3] == 2146;
  fail &= bins[4] == 0;

  DestroyRegions();
  BoxLib::Finalize();
  if (fail) {
    BoxLib::Abort();
  }
  return 0;
}
