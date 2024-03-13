/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <iostream>
#include <fstream>
#include <iomanip>

// closing DSO objects
#include "VerboseObject_objs.hh"

using std::cout;
using std::endl;
#include <ParmParse.H>

#include <Material.H>
#include <MatFiller.H>
#include <Region.H>

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

void SetRegions()
{
  Array<Real> lo(2), hi(2);
  std::string r_name, r_purpose;
  r_name = "SoilLower";
  r_purpose = "all";
  lo[0] = 0;
  lo[1] = 0;
  hi[0] = 40;
  hi[1] = 10;
  regions["SoilLower"] = new BoxRegion(r_name,r_purpose,lo,hi);

  r_name = "SoilRight";
  r_purpose = "all";
  lo[0] = 12;
  lo[1] = 10;
  hi[0] = 40;
  hi[1] = 18;
  regions["SoilRight"] = new BoxRegion(r_name,r_purpose,lo,hi);

  r_name = "SoilUpper";
  r_purpose = "all";
  lo[0] = 0;
  lo[1] = 18;
  hi[0] = 40;
  hi[1] = 24;
  regions["SoilUpper"] = new BoxRegion(r_name,r_purpose,lo,hi);

  r_name = "TankConcFloor";
  r_purpose = "all";
  lo[0] = 0.5;
  lo[1] = 10;
  hi[0] = 12;
  hi[1] = 10.3;
  regions["TankConcFloor"] = new BoxRegion(r_name,r_purpose,lo,hi);

  r_name = "TankConcRoof1";
  r_purpose = "all";
  lo[0] = 0;
  lo[1] = 17.7;
  hi[0] = 11.67;
  hi[1] = 18;
  regions["TankConcRoof1"] = new BoxRegion(r_name,r_purpose,lo,hi);

  r_name = "TankConcRoof2";
  r_purpose = "all";
  lo[0] = 11.69;
  lo[1] = 17.7;
  hi[0] = 12;
  hi[1] = 18;
  regions["TankConcRoof2"] = new BoxRegion(r_name,r_purpose,lo,hi);

  r_name = "TankConcWall";
  r_purpose = "all";
  lo[0] = 11.7;
  lo[1] = 10.3;
  hi[0] = 12;
  hi[1] = 17.7;
  regions["TankConcWall"] = new BoxRegion(r_name,r_purpose,lo,hi);

  r_name = "TankFFfloor";
  r_purpose = "all";
  lo[0] = 0;
  lo[1] = 10;
  hi[0] = 0.5;
  hi[1] = 10.31;
  regions["TankFFfloor"] = new BoxRegion(r_name,r_purpose,lo,hi);

  r_name = "TankFFwall";
  r_purpose = "all";
  lo[0] = 11.67;
  lo[1] = 10.33;
  hi[0] = 11.69;
  hi[1] = 18;
  regions["TankFFwall"] = new BoxRegion(r_name,r_purpose,lo,hi);

  r_name = "TankGrout";
  r_purpose = "all";
  lo[0] = 0;
  lo[1] = 10.33;
  hi[0] = 11.67;
  hi[1] = 17.69;
  regions["TankGrout"] = new BoxRegion(r_name,r_purpose,lo,hi);

  r_name = "TankLinerFloor";
  r_purpose = "all";
  lo[0] = 0.5;
  lo[1] = 10.3;
  hi[0] = 11.7;
  hi[1] = 10.31;
  regions["TankLinerFloor"] = new BoxRegion(r_name,r_purpose,lo,hi);

  r_name = "TankLinerRoof";
  r_purpose = "all";
  lo[0] = 0;
  lo[1] = 17.69;
  hi[0] = 11.67;
  hi[1] = 17.7;
  regions["TankLinerRoof"] = new BoxRegion(r_name,r_purpose,lo,hi);

  r_name = "TankLinerWall";
  r_purpose = "all";
  lo[0] = 11.69;
  lo[1] = 10.31;
  hi[0] = 11.7;
  hi[1] = 17.7;
  regions["TankLinerWall"] = new BoxRegion(r_name,r_purpose,lo,hi);

  r_name = "TankWaste";
  r_purpose = "all";
  lo[0] = 0;
  lo[1] = 10.31;
  hi[0] = 11.69;
  hi[1] = 10.33;
  regions["TankWaste"] = new BoxRegion(r_name,r_purpose,lo,hi);
}

PArray<Material>
SetMaterials()
{
  PArray<Material> materials(5, PArrayManage);
  Array<std::string> region_names;
  Array<const Region*> regionset;
  region_names.push_back("SoilLower");
  region_names.push_back("SoilRight");
  region_names.push_back("SoilUpper");
  regionset = RegionPtrArray(region_names);
  materials.set(0,new Material("Soil",regionset));
  region_names.clear();
  region_names.push_back("TankConcFloor");
  region_names.push_back("TankConcRoof1");
  region_names.push_back("TankConcRoof2");
  region_names.push_back("TankConcWall");
  regionset = RegionPtrArray(region_names);
  materials.set(1,new Material("TankConc",regionset));
  region_names.clear();
  region_names.push_back("TankFFfloor");
  region_names.push_back("TankFFwall");
  region_names.push_back("TankWaste");
  regionset = RegionPtrArray(region_names);
  materials.set(2,new Material("TankFF",regionset));
  region_names.clear();
  region_names.push_back("TankGrout");
  regionset = RegionPtrArray(region_names);
  materials.set(3,new Material("TankGrout",regionset));
  region_names.clear();
  region_names.push_back("TankLinerFloor");
  region_names.push_back("TankLinerRoof");
  region_names.push_back("TankLinerWall");
  regionset = RegionPtrArray(region_names);
  materials.set(4,new Material("TankLiner",regionset));
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
  BoxLib::Initialize(argc,argv,false);
  ParmParse pp;

  int nLevs = 3;
  Array<int> n_cells(BL_SPACEDIM);
  n_cells[0] = 40; n_cells[1] = 24;

  Array<int> rRatio(nLevs-1,4);
  Array<IntVect> refRatio(nLevs-1);
  for (int lev=0; lev<nLevs-1; ++lev) {
    refRatio[lev] = rRatio[lev] * IntVect::TheUnitVector();
  }

  int coord = 0;
  Array<int> is_per(BL_SPACEDIM,0);
  Array<Real> prob_lo(BL_SPACEDIM,0);
  Array<Real> prob_hi(BL_SPACEDIM);
  D_EXPR(prob_hi[0]=40, prob_hi[1]=24, prob_hi[2]=10);
  const RealBox rb(prob_lo.dataPtr(),prob_hi.dataPtr());

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
    geomArray[lev] = Geometry(domain,&rb,coord,is_per.dataPtr());
  }

  Region::geometry_eps = 1.e-6;
  Region::domlo.resize(BL_SPACEDIM);
  Region::domhi.resize(BL_SPACEDIM);
  for (int d=0; d<BL_SPACEDIM; ++d) {
    Region::domlo[d] = Geometry::ProbLo()[d];
    Region::domhi[d] = Geometry::ProbHi()[d];
  }

  SetRegions();
  PArray<Material> materials = SetMaterials();

  MatFiller matFiller(geomArray,refRatio,materials);

  bool fail = false;


  const std::map<std::string,int>& mat_map = matFiller.MatIdx();

  Array<int> bins(mat_map.size(),0);

  for (int lev=0; lev<nLevs; ++lev) {
    const iMultiFab& mf = matFiller.MaterialID(lev);
    if (lev<nLevs) {
      const BoxArray& ba_mixed = matFiller.Mixed(lev);
      if (ba_mixed.size()>0) {
        int maxVal = -1;
        for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
          const IArrayBox& fab = mf[mfi];
          std::vector< std::pair<int,Box> > isects = ba_mixed.intersections(mfi.validbox());
          for (int ii = 0, N = isects.size(); ii < N; ii++)
          {
            maxVal = std::max(maxVal,fab.max(isects[ii].second,0));
          }
        }
        ParallelDescriptor::ReduceIntMax(maxVal);
        fail = (maxVal>-1);
      }
    }

    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
      const Box& vbox = mfi.validbox();
      const IArrayBox& fab = mf[mfi];
      for (IntVect iv=vbox.smallEnd(), BIG=vbox.bigEnd(); iv<=BIG; vbox.next(iv)) {
        int val = fab(iv,0);
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
