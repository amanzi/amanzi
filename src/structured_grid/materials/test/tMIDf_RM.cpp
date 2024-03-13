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
#include <RegionManager.H>

PArray<Material>
SetMaterials(RegionManager& rm)
{
  ParmParse pp;
  int nmat = pp.countval("materials");
  PArray<Material> materials;
  if (nmat) {
    materials.resize(nmat, PArrayManage);
    Array<std::string> material_names, region_names;
    pp.getarr("materials",material_names,0,nmat);
    for (int i=0; i<nmat; ++i) {
      const std::string prefix("materials." + material_names[i]);
      ParmParse ppr(prefix.c_str());
      ppr.getarr("regions",region_names,0,ppr.countval("regions"));
      materials.set(i,new Material(material_names[i],rm.RegionPtrArray(region_names)));
    }
  }
  return materials;
}

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
  MatFiller matFiller(geomArray,refRatio,SetMaterials(rm));

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

  BoxLib::Finalize();
  if (fail) {
    BoxLib::Abort();
  }
  return 0;
}
