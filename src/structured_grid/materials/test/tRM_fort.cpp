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
#include <VisMF.H>

#include <RockManager.H>

#include <tRM_fort_F.H>

static void
GradFill (MultiFab&          mf,
          const Array<Real>& grad,
          const Geometry&    geom)
{
  for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
    FArrayBox& fab = mf[mfi];
    const Box& box = mfi.validbox();
    for (IntVect iv=box.smallEnd(), End=box.bigEnd(); iv<=End; box.next(iv)) {
      Real val = 0;
      for (int d=0; d<BL_SPACEDIM; ++d) {
        val += grad[d]*(iv[d]+0.5)*geom.CellSize()[d];
      }
      fab(iv,0) = val;
    }
  }
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

  int nGrowRM = 3;
  RegionManager rm;
  RockManager rockManager(&rm);
  rockManager.FinalizeBuild(geomArray,refRatio,nGrowRM,false);

  Real time = 0.2;
  int level = 0;
  int nGrow = 0;
  bool ignore_mixed = true;
  BoxArray ba(geomArray[level].Domain());
  int maxSize=ba[0].length(0);  pp.query("maxSize",maxSize);
  ba.maxSize(maxSize);

  iMultiFab matID(ba,1,0);
  rockManager.GetMaterialID(level,matID,nGrow,ignore_mixed);

  MultiFab pc(ba,1,nGrow);
  MultiFab sat(ba,1,nGrow);

  MultiFab dsdp(ba,1,nGrow);
  MultiFab pc1(ba,1,nGrow);
  MultiFab kr(ba,1,nGrow);
  Array<Real> gradp(BL_SPACEDIM,0);
  gradp[1] = 101325. / (geomArray[level].ProbHi()[1]-geomArray[level].ProbLo()[1]);
  GradFill(pc,gradp,geomArray[level]);
  for (MFIter mfi(pc); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();
    rockManager.InverseCapillaryPressure(pc[mfi].dataPtr(),matID[mfi].dataPtr(),time,sat[mfi].dataPtr(),box.numPts());
    rockManager.CapillaryPressure(sat[mfi].dataPtr(),matID[mfi].dataPtr(),time,pc1[mfi].dataPtr(),box.numPts());
    rockManager.DInverseCapillaryPressure(sat[mfi].dataPtr(),matID[mfi].dataPtr(),time,dsdp[mfi].dataPtr(),box.numPts());
    rockManager.RelativePermeability(sat[mfi].dataPtr(),matID[mfi].dataPtr(),time,kr[mfi].dataPtr(),box.numPts());
  }

  MultiFab pc2(ba,1,nGrow);
  MultiFab sat2(ba,1,nGrow);
  MultiFab kr2(ba,1,nGrow);
  MultiFab dsdp2(ba,1,nGrow);
  int rmID = rockManager.ID();
  for (MFIter mfi(pc); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();
    const FArrayBox& s = sat[mfi];
    const IArrayBox& m = matID[mfi];
    FArrayBox& p = pc2[mfi];
    FArrayBox& s2 = sat2[mfi];
    FArrayBox& k2 = kr2[mfi];
    FArrayBox& ds = dsdp2[mfi];

    FORT_TEST_PCAP (&rmID,&time,box.loVect(),box.hiVect(),
		    s.dataPtr(), ARLIM(s.loVect()), ARLIM(s.hiVect()),
		    m.dataPtr(), ARLIM(m.loVect()), ARLIM(m.hiVect()),
		    p.dataPtr(), ARLIM(p.loVect()), ARLIM(p.hiVect()));

    FORT_TEST_INVPCAP (&rmID,&time,box.loVect(),box.hiVect(),
		       p.dataPtr(),  ARLIM(p.loVect()),  ARLIM(p.hiVect()),
		       m.dataPtr(),  ARLIM(m.loVect()),  ARLIM(m.hiVect()),
		       s2.dataPtr(), ARLIM(s2.loVect()), ARLIM(s2.hiVect()));

    FORT_TEST_RELPERM (&rmID,&time,box.loVect(),box.hiVect(),
		       s.dataPtr(),  ARLIM(s.loVect()),  ARLIM(s.hiVect()),
		       m.dataPtr(),  ARLIM(m.loVect()),  ARLIM(m.hiVect()),
		       k2.dataPtr(), ARLIM(k2.loVect()), ARLIM(k2.hiVect()));

    FORT_TEST_DSDPCAP (&rmID,&time,box.loVect(),box.hiVect(),
		       s.dataPtr(),  ARLIM(s.loVect()),  ARLIM(s.hiVect()),
		       m.dataPtr(),  ARLIM(m.loVect()),  ARLIM(m.hiVect()),
		       ds.dataPtr(), ARLIM(ds.loVect()), ARLIM(ds.hiVect()));
  }

  MultiFab::Subtract(pc2,pc1,0,0,1,nGrow);
  Real diff_pc = pc2.norm0(0);

  MultiFab::Subtract(sat2,sat,0,0,1,nGrow);
  Real diff_sat = sat2.norm0(0);

  MultiFab::Subtract(kr2,kr,0,0,1,nGrow);
  Real diff_kr = kr2.norm0(0);

  MultiFab::Subtract(dsdp2,dsdp,0,0,1,nGrow);
  Real diff_dsdp = dsdp2.norm0(0);

  Real max_diff = 1.e-12;
  bool fail = ( (   diff_pc > max_diff) ||
                (  diff_sat > max_diff) ||
                (   diff_kr > max_diff) ||
                ( diff_dsdp > max_diff) );
  if (fail) {
    BoxLib::Abort();
  }

  return 0;
}
