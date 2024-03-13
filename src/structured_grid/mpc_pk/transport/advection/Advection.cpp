/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <fstream>

#include <Advection.H>
#include <Advection_F.H>

int
Advection::nGrowHyp()
{
  int nGrowHyp;
  BDS_SET_NGROWHYP(&nGrowHyp);
  return nGrowHyp;
}

int
Advection::nGrowForce()
{
  int nGrowForce;
  BDS_SET_NGROWFORCE(&nGrowForce);
  return nGrowForce;
}

int
Advection::nWork()
{
  int nWork;
  BDS_SET_NWORK(&nWork);
  return nWork;
}

void
Advection::FluxDivergence(const FArrayBox& CoCC, FArrayBox& CnCC,       int Ccomp,    int Cng,
                          const FArrayBox& SoCC, const FArrayBox& SnCC, int Scomp,    int Sng,
                          const PArray<FArrayBox>& UEC,                 int Ucomp,    int Ung,
                          FArrayBox& FDivCC,                            int FDivcomp, int FDivng,
                          PArray<FArrayBox>& FluxEC,                    int Flcomp,   int Flng,
                          const FArrayBox& FcCC,                        int Fccomp,   int Fcng,
                          const FArrayBox& FsCC,                        int Fscomp,   int Fsng,
                          const FArrayBox& VolCC, const PArray<FArrayBox>& AEC, const FArrayBox& PhiCC,
                          const Box& vbox, const Real* dx, Real dt, const int* bc, int nc,
                          bool is_conservative)
{
  int nGrowH = nGrowHyp();
  int nGrowF = nGrowForce();

  BL_ASSERT(Cng >= nGrowH);
  BL_ASSERT(Sng >= nGrowF);
  BL_ASSERT(Ung >= 0);
  BL_ASSERT(FDivng >= 0);
  BL_ASSERT(Flng >= 0);
  BL_ASSERT(Fcng >= nGrowF);
  BL_ASSERT(Fsng >= nGrowF);

  BL_ASSERT(CoCC.nComp()>=Ccomp+nc);
  BL_ASSERT(CnCC.nComp()>=Ccomp+nc);
  BL_ASSERT(SoCC.nComp()>=Scomp+1);
  BL_ASSERT(SnCC.nComp()>=Scomp+1);
  BL_ASSERT(FDivCC.nComp()>=FDivcomp+nc);
  BL_ASSERT(FcCC.nComp()>=Fccomp+nc);
  BL_ASSERT(FsCC.nComp()>=Fscomp+1);

  const Real* co_dat  = CoCC.dataPtr(Ccomp);
  Real*       cex_dat = FluxEC[0].dataPtr(Flcomp);
  Real*       cey_dat = FluxEC[1].dataPtr(Flcomp);
  const Real* u_dat   = UEC[0].dataPtr(Ucomp);
  const Real* v_dat   = UEC[1].dataPtr(Ucomp);
  const Real* ax_dat  = AEC[0].dataPtr();
  const Real* ay_dat  = AEC[1].dataPtr();
  const Real* vol_dat = VolCC.dataPtr();
  Real*       aofs_dat = FDivCC.dataPtr(FDivcomp);

  Box Slbox = BoxLib::grow(vbox,nGrowF);

  int NWORK = nWork();
  FArrayBox slope(Slbox,NWORK);
  Real* slope_dat = slope.dataPtr();

  FArrayBox capInv(Slbox,1);
  capInv.copy(PhiCC);
  capInv.mult(SoCC);
  capInv.invert(1);
  const Real* ci_dat = capInv.dataPtr();

  FArrayBox F(Slbox,nc); // F = (Fc - c*Fs) / (phi.s) evaluated at t_old
  for (int i=0; i<nc; ++i) {
    F.copy(FsCC,Fscomp,i,1);
    F.mult(CoCC,Ccomp+i,i,1);
    F.mult(-1,i,1);
    F.plus(FcCC,Fccomp+i,i,1);
    F.mult(capInv,0,i,1);
  }
  const Real* F_dat = F.dataPtr();
  int is_cons = is_conservative ? 1 : 0;

#if BL_SPACEDIM == 3
  Real*       cez_dat = FluxEC[2].dataPtr(Flcomp);
  const Real* w_dat   = UEC[2].dataPtr(Ucomp);
  const Real* az_dat  = AEC[2].dataPtr();
#endif

  BDS_EDGE_STATES(co_dat, &Cng, ci_dat, &nGrowF, D_DECL(cex_dat, cey_dat, cez_dat), &Flng,
                  D_DECL(u_dat, v_dat, w_dat), &Ung, F_dat, &Fcng,
                  slope_dat, &nGrowF, &NWORK, dx, &dt, &nc, &is_cons, vbox.loVect(), vbox.hiVect(), bc);

  ADV_FDIV(aofs_dat, &FDivng, D_DECL(cex_dat, cey_dat, cez_dat), &Flng,
           D_DECL(u_dat, v_dat, w_dat), &Ung, D_DECL(ax_dat, ay_dat, az_dat), vol_dat, &nc, vbox.loVect(), vbox.hiVect());

}

void
Advection::AdvUpdate(const FArrayBox& CoCC, FArrayBox& CnCC,       int Ccomp,    int Cng,
                     const FArrayBox& SoCC, const FArrayBox& SnCC, int Scomp,    int Sng,
                     FArrayBox& FDivCC,                            int FDivcomp, int FDivng,
                     PArray<FArrayBox>& FluxEC,                    int Flcomp,   int Flng,
                     const FArrayBox& PhiCC,                                     int Phing,
                     const Box& vbox, Real dt, int nc)
{
  BL_ASSERT(CoCC.nComp()>=Ccomp+nc);
  BL_ASSERT(CnCC.nComp()>=Ccomp+nc);
  BL_ASSERT(SoCC.nComp()>=Scomp+1);
  BL_ASSERT(SnCC.nComp()>=Scomp+1);
  BL_ASSERT(FDivCC.nComp()>=FDivcomp+nc);

  Real*       cn_dat    = CnCC.dataPtr(Ccomp);
  const Real* co_dat    = CoCC.dataPtr(Ccomp);
  const Real* so_dat    = SoCC.dataPtr(Ccomp);
  const Real* sn_dat    = SnCC.dataPtr(Ccomp);
  const Real* xflux_dat = FluxEC[0].dataPtr(Flcomp);
  const Real* yflux_dat = FluxEC[1].dataPtr(Flcomp);
  const Real* phi_dat   = PhiCC.dataPtr();
  Real*       aofs_dat  = FDivCC.dataPtr(FDivcomp);

  ADV_UPDATE(co_dat, cn_dat, &Cng, so_dat, sn_dat, &Sng, aofs_dat, &FDivng, phi_dat, &Phing,
             &dt, &nc, vbox.loVect(), vbox.hiVect());
}

void
Advection::EstimateMaxEigenvalues(const FArrayBox&        SatCC, int Scomp, int Satng,
                                  D_DECL(const FArrayBox& Ux,
                                         const FArrayBox& Uy,
                                         const FArrayBox& Uz),     int Ucomp, int Ung,
                                  const FArrayBox&        PhiCC,            int Phing,
                                  const Box&              vbox,
                                  Real*                   eigmax)
{
  BL_ASSERT(SatCC.nComp()>=Scomp);

  const Real* s_dat    = SatCC.dataPtr(Scomp);
  const Real* umac_dat = Ux.dataPtr(Ucomp);
  const Real* vmac_dat = Uy.dataPtr(Ucomp);
  const Real* phi_dat  = PhiCC.dataPtr();

#if BL_SPACEDIM == 3
  const Real* wmac_dat = Uz.dataPtr(Ucomp);
#endif

  BDS_EST_EIGEN(s_dat,&Satng,D_DECL(umac_dat,vmac_dat,wmac_dat),&Ung,phi_dat,&Phing,vbox.loVect(),vbox.hiVect(),eigmax);
}
