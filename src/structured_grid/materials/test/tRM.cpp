#include <iostream>
#include <fstream>
#include <iomanip>
using std::cout;
using std::endl;
#include <ParmParse.H>
#include <VisMF.H>

#include <RockManager.H>

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

  RegionManager rm;
  RockManager rockManager(&rm,geomArray,refRatio);

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

  Array<Real> grads(BL_SPACEDIM,0);
  grads[1] = 0.5 / (geomArray[level].ProbHi()[1]-geomArray[level].ProbLo()[1]);
  GradFill(sat,grads,geomArray[level]);
  sat.plus(0.5,0);
  for (MFIter mfi(pc); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();
    rockManager.CapillaryPressure(sat[mfi].dataPtr(),matID[mfi].dataPtr(),time,pc[mfi].dataPtr(),box.numPts());
  }

  MultiFab pc1(ba,1,nGrow);
  MultiFab dsdp(ba,1,nGrow);
  MultiFab relperm(ba,1,nGrow);
  MultiFab c(ba,1,nGrow);
  Array<Real> gradp(BL_SPACEDIM,0);
  gradp[1] = 101325. / (geomArray[level].ProbHi()[1]-geomArray[level].ProbLo()[1]);
  GradFill(pc,gradp,geomArray[level]);
  for (MFIter mfi(pc); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();
    rockManager.InverseCapillaryPressure(pc[mfi].dataPtr(),matID[mfi].dataPtr(),time,sat[mfi].dataPtr(),box.numPts());
    rockManager.CapillaryPressure(sat[mfi].dataPtr(),matID[mfi].dataPtr(),time,pc1[mfi].dataPtr(),box.numPts());
    rockManager.DInverseCapillaryPressure(sat[mfi].dataPtr(),matID[mfi].dataPtr(),time,dsdp[mfi].dataPtr(),box.numPts());
    rockManager.RelativePermeability(sat[mfi].dataPtr(),matID[mfi].dataPtr(),time,relperm[mfi].dataPtr(),box.numPts());

    c[mfi].setVal(0);
    for (IntVect iv(box.smallEnd()), End=box.bigEnd(); iv<=End; box.next(iv)) {
      if (iv[1] < End[1] && iv[1] > box.smallEnd()[1]) {
        c[mfi](iv,0) = 0.5*(sat[mfi](iv+IntVect(0,1),0) - sat[mfi](iv-IntVect(0,1),0))/(pc[mfi](iv+IntVect(0,1),0) - pc[mfi](iv-IntVect(0,1),0));
      }
    }

  }

#if 1
  MultiFab pc2(ba,1,nGrow);
  MultiFab::Copy(pc2,pc,0,0,1,nGrow);
  MultiFab::Subtract(pc2,pc1,0,0,1,nGrow);

  VisMF::Write(sat,"s");
  VisMF::Write(pc,"p");
  VisMF::Write(pc1,"p1");
  VisMF::Write(pc2,"dp");
  VisMF::Write(dsdp,"dsdp");
  VisMF::Write(c,"c");
  VisMF::Write(relperm,"k");
#endif
  return 0;
}
