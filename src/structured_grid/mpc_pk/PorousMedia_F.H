/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#ifndef _PorousMedia_F_H_
#define _PorousMedia_F_H_

#ifdef BL_LANG_FORT
#    define FORT_FILCC             filcc
#    define FORT_AVGDOWN           avgdown
#    define FORT_SUMMASS           summass
#    define FORT_CEN2EDG           cen2edg
#    define FORT_EDGE_INTERP       edge_interp
#    define FORT_PC_EDGE_INTERP    pc_edge_interp
#    define FORT_MK_CPL            mkcpl
#    define FORT_MK_INV_CPL        mkinvcpl
#    define FORT_MK_LAMBDA         mklambda
#    define FORT_MK_DLAMBDA        mkdlambda
#    define FORT_RICHARD_AJAC      richard_ajac
#    define FORT_RICHARD_NJAC      richard_njac
#    define FORT_RICHARD_NJAC2     richard_njac2
#    define FORT_RICHARD_ALPHA     richard_alpha
#    define FORT_RICHARD_VELBC     richard_velbc
#    define FORT_AVG_UMAC          avgumac
#    define FORT_CPY_UMAC          cpyumac
#    define FORT_LINEAR_PRESSURE   linearpres
#    define FORT_HYD_HEAD          hydhead
#    define FORT_STEADYSTATE       steadystate
#    define FORT_MAX_TRACDIFF_DT   maxtracdiff
#    define FORT_TENSORDIFF        tensordiff
#    define FORT_VISCEXTRAP        viscextrap
#    define FORT_INITKEDGE         initkedge
#else
#  if defined(BL_FORT_USE_UPPERCASE)
#    define FORT_FILCC             FILCC
#    define FORT_AVGDOWN           AVGDOWN
#    define FORT_SUMMASS           SUMMASS
#    define FORT_CEN2EDG           CEN2EDG
#    define FORT_EDGE_INTERP       EDGE_INTERP
#    define FORT_PC_EDGE_INTERP    PC_EDGE_INTERP
#    define FORT_MK_CPL            MKCPL
#    define FORT_MK_INV_CPL        MKINVCPL
#    define FORT_MK_LAMBDA         MKLAMBDA
#    define FORT_MK_DLAMBDA        MKDLAMBDA
#    define FORT_RICHARD_AJAC      RICHARD_AJAC
#    define FORT_RICHARD_NJAC      RICHARD_NJAC
#    define FORT_RICHARD_NJAC2     RICHARD_NJAC2
#    define FORT_RICHARD_ALPHA     RICHARD_ALPHA
#    define FORT_RICHARD_VELBC     RICHARD_VELBC
#    define FORT_AVG_UMAC          AVGUMAC
#    define FORT_CPY_UMAC          CPYUMAC
#    define FORT_LINEAR_PRESSURE   LINEARPRES
#    define FORT_HYD_HEAD          HYDHEAD
#    define FORT_STEADYSTATE       STEADYSTATE
#    define FORT_PATCH_GHOST       PATCH_GHOST
#    define FORT_MAX_TRACDIFF_DT   MAXTRACDIFF
#    define FORT_TENSORDIFF        TENSORDIFF
#    define FORT_VISCEXTRAP        VISCEXTRAP
#    define FORT_INITKEDGE         INITKEDGE
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_FILCC             filcc
#    define FORT_AVGDOWN           avgdown
#    define FORT_SUMMASS           summass
#    define FORT_CEN2EDG           cen2edg
#    define FORT_EDGE_INTERP       edge_interp
#    define FORT_PC_EDGE_INTERP    pc_edge_interp
#    define FORT_MK_CPL            mkcpl
#    define FORT_MK_INV_CPL        mkinvcpl
#    define FORT_MK_LAMBDA         mklambda
#    define FORT_MK_DLAMBDA        mkdlambda
#    define FORT_RICHARD_AJAC      richard_ajac
#    define FORT_RICHARD_NJAC      richard_njac
#    define FORT_RICHARD_NJAC2     richard_njac2
#    define FORT_RICHARD_ALPHA     richard_alpha
#    define FORT_RICHARD_VELBC     richard_velbc
#    define FORT_AVG_UMAC          avgumac
#    define FORT_CPY_UMAC          cpyumac
#    define FORT_LINEAR_PRESSURE   linearpres
#    define FORT_HYD_HEAD          hydhead
#    define FORT_STEADYSTATE       steadystate
#    define FORT_MAX_TRACDIFF_DT   maxtracdiff
#    define FORT_TENSORDIFF        tensordiff
#    define FORT_VISCEXTRAP        viscextrap
#    define FORT_INITKEDGE         initkedge
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_FILCC             filcc_
#    define FORT_AVGDOWN           avgdown_
#    define FORT_SUMMASS           summass_
#    define FORT_CEN2EDG           cen2edg_
#    define FORT_EDGE_INTERP       edge_interp_
#    define FORT_PC_EDGE_INTERP    pc_edge_interp_
#    define FORT_MK_CPL            mkcpl_
#    define FORT_MK_INV_CPL        mkinvcpl_
#    define FORT_MK_LAMBDA         mklambda_
#    define FORT_MK_DLAMBDA        mkdlambda_
#    define FORT_RICHARD_AJAC      richard_ajac_
#    define FORT_RICHARD_NJAC      richard_njac_
#    define FORT_RICHARD_NJAC2     richard_njac2_
#    define FORT_RICHARD_ALPHA     richard_alpha_
#    define FORT_RICHARD_VELBC     richard_velbc_
#    define FORT_AVG_UMAC          avgumac_
#    define FORT_CPY_UMAC          cpyumac_
#    define FORT_LINEAR_PRESSURE   linearpres_
#    define FORT_HYD_HEAD          hydhead_
#    define FORT_STEADYSTATE       steadystate_
#    define FORT_MAX_TRACDIFF_DT   maxtracdiff_
#    define FORT_TENSORDIFF        tensordiff_
#    define FORT_VISCEXTRAP        viscextrap_
#    define FORT_INITKEDGE         initkedge_
#  endif

#include <ArrayLim.H>
#include "Prob_PM_F.H"

extern "C" {

    void FORT_FILCC (const Real * q, ARLIM_P(q_lo), ARLIM_P(q_hi),
                     const int * domlo, const int * domhi,
                     const Real * dx_crse, const Real * xlo,
                     const int * bc);

    void FORT_AVGDOWN (const Real* crse, ARLIM_P(clo),     ARLIM_P(chi),
                       const int* nvar,
                       const Real* fine, ARLIM_P(fine_lo), ARLIM_P(fine_hi),
                       const Real* cvol, ARLIM_P(cv_lo),   ARLIM_P(cv_hi),
                       const Real* fvol, ARLIM_P(fv_lo),   ARLIM_P(fv_hi),
		       const int* lo, const int* hi, const int* ratio);

    void FORT_SUMMASS (const Real* dat, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                       ARLIM_P(grid_lo), ARLIM_P(grid_hi),
                       const Real* delta, Real* mass);

    void FORT_CEN2EDG(const int* lo, const int* hi,
                      ARLIM_P(clo), ARLIM_P(chi), const Real* cfabdat,
                      ARLIM_P(elo), ARLIM_P(ehi), Real* efabdat,
                      const int* nc, const int* dir, const int* isharm);

    void FORT_EDGE_INTERP(const int* flo, const int* fhi,
                          const int* nc, const int* refRatio, const int* dir,
                          Real*       fine, ARLIM_P(flo), ARLIM_P(fhi));

    void FORT_PC_EDGE_INTERP(const int* lo, const int* hi,
			     const int* nc, const int* refRatio, const int* dir,
			     const Real* crse, ARLIM_P(clo), ARLIM_P(chi),
			     Real*       fine, ARLIM_P(flo), ARLIM_P(fhi));

    void FORT_MAX_TRACDIFF_DT(const int* lo, const int* hi,
                              const Real* psv,  ARLIM_P(psv),   ARLIM_P(psv),
                              const Real* psadx,ARLIM_P(psadx), ARLIM_P(psadx),
                              const Real* psady,ARLIM_P(psady), ARLIM_P(psady),
#if (BL_SPACEDIM==3)
                              const Real* psadz,ARLIM_P(psadz), ARLIM_P(psadz),
#endif
                              const Real* delta, Real* max_dt);
   void FORT_MK_CPL (const Real* pc,  ARLIM_P(pclo),  ARLIM_P(pchi),
		     const Real* s,   ARLIM_P(slo),   ARLIM_P(shi),
		     const Real* phi, ARLIM_P(philo), ARLIM_P(phihi),
		     const Real* kp,  ARLIM_P(kplo),  ARLIM_P(kphi),
		     const Real* cf,  ARLIM_P(cflo),  ARLIM_P(cfhi),
		     const int* ncf, const int* lo, const int* hi, const int* bc);

   void FORT_MK_INV_CPL (const Real* pc,  ARLIM_P(pclo),  ARLIM_P(pchi),
			 const Real* s,   ARLIM_P(slo),   ARLIM_P(shi),
			 const Real* phi, ARLIM_P(philo), ARLIM_P(phihi),
			 const Real* kp,  ARLIM_P(kplo),  ARLIM_P(kphi),
			 const Real* cf,  ARLIM_P(cflo),  ARLIM_P(cfhi), const int* ncf);

  void FORT_MK_LAMBDA (const Real* lbd, ARLIM_P(lbdlo), ARLIM_P(lbdhi),
		       const Real* s,   ARLIM_P(slo),   ARLIM_P(shi),
		       const Real* kr,  ARLIM_P(krlo),  ARLIM_P(krhi),
		       const int* nkr);

  void FORT_MK_DLAMBDA (const Real* lbd,  ARLIM_P(lbdlo),  ARLIM_P(lbdhi),
			const Real* s,   ARLIM_P(slo),   ARLIM_P(shi),
			const Real* kr,  ARLIM_P(krlo),  ARLIM_P(krhi),
			const int* nkr);

  void FORT_RICHARD_AJAC(const Real* S, ARLIM_P(slo), ARLIM_P(shi),

			 const Real* dfx, ARLIM_P(dfxlo), ARLIM_P(dfxhi),
			 const Real* dfy, ARLIM_P(dfylo), ARLIM_P(dfyhi),
#if (BL_SPACEDIM==3)
			 const Real* dfz, ARLIM_P(dfzlo), ARLIM_P(dfzhi),
#endif
			 const Real* ux, ARLIM_P(uxlo), ARLIM_P(uxhi),
			 const Real* uy, ARLIM_P(uylo), ARLIM_P(uyhi),
#if (BL_SPACEDIM==3)
			 const Real* uz, ARLIM_P(uzlo), ARLIM_P(uzhi),
#endif
			 const Real* kpx, ARLIM_P(kpxlo), ARLIM_P(kpxhi),
			 const Real* kpy, ARLIM_P(kpylo), ARLIM_P(kpyhi),
#if (BL_SPACEDIM==3)
			 const Real* kpz, ARLIM_P(kpzlo), ARLIM_P(kpzhi),
#endif
			 const Real* lbd, ARLIM_P(lbdlo), ARLIM_P(lbdhi),
			 const Real* pc, ARLIM_P(pclo), ARLIM_P(pchi),
			 const Real* P, ARLIM_P(plo), ARLIM_P(phi),
			 const Real* K, ARLIM_P(klo), ARLIM_P(khi),
			 const Real* kr,  ARLIM_P(krlo),  ARLIM_P(krhi),
			 const int* nkr,
			 const Real* cf,  ARLIM_P(cflo),  ARLIM_P(cfhi),
			 const int* ncf,
			 const int* lo, const int* hi,
			 const int* domlo, const int* domhi,
			 const Real* dx, const int* bc,
			 const int* inflow_bc_lo, const int* inflow_bc_hi,
			 const int* do_upwind);

  void FORT_RICHARD_NJAC(const Real* S, ARLIM_P(slo), ARLIM_P(shi),

			 const Real* dfx, ARLIM_P(dfxlo), ARLIM_P(dfxhi),
			 const Real* dfy, ARLIM_P(dfylo), ARLIM_P(dfyhi),
#if (BL_SPACEDIM==3)
			 const Real* dfz, ARLIM_P(dfzlo), ARLIM_P(dfzhi),
#endif
			 const Real* ux, ARLIM_P(uxlo), ARLIM_P(uxhi),
			 const Real* uy, ARLIM_P(uylo), ARLIM_P(uyhi),
#if (BL_SPACEDIM==3)
			 const Real* uz, ARLIM_P(uzlo), ARLIM_P(uzhi),
#endif
			 const Real* kpx, ARLIM_P(kpxlo), ARLIM_P(kpxhi),
			 const Real* kpy, ARLIM_P(kpylo), ARLIM_P(kpyhi),
#if (BL_SPACEDIM==3)
			 const Real* kpz, ARLIM_P(kpzlo), ARLIM_P(kpzhi),
#endif
			 const Real* lbd, ARLIM_P(lbdlo), ARLIM_P(lbdhi),
			 const Real* pc, ARLIM_P(pclo), ARLIM_P(pchi),
			 const Real* P, ARLIM_P(plo), ARLIM_P(phi),
			 const Real* K, ARLIM_P(klo), ARLIM_P(khi),
			 const Real* kr,  ARLIM_P(krlo),  ARLIM_P(krhi),
			 const int* nkr,
			 const Real* cf,  ARLIM_P(cflo),  ARLIM_P(cfhi),
			 const int* ncf,
			 const int* lo, const int* hi,
			 const int* domlo, const int* domhi,
			 const Real* dx, const int* bc,
			 const int* inflow_bc_lo, const int* inflow_bc_hi,
			 const Real* deps, const int* do_upwind);

  void FORT_RICHARD_NJAC2(const Real* dfx, ARLIM_P(dfxlo), ARLIM_P(dfxhi),
			  const Real* dfy, ARLIM_P(dfylo), ARLIM_P(dfyhi),
#if (BL_SPACEDIM==3)
			  const Real* dfz, ARLIM_P(dfzlo), ARLIM_P(dfzhi),
#endif
			  const Real* ux, ARLIM_P(uxlo), ARLIM_P(uxhi),
			  const Real* uy, ARLIM_P(uylo), ARLIM_P(uyhi),
#if (BL_SPACEDIM==3)
			  const Real* uz, ARLIM_P(uzlo), ARLIM_P(uzhi),
#endif
			  const Real* kpx, ARLIM_P(kpxlo), ARLIM_P(kpxhi),
			  const Real* kpy, ARLIM_P(kpylo), ARLIM_P(kpyhi),
#if (BL_SPACEDIM==3)
			  const Real* kpz, ARLIM_P(kpzlo), ARLIM_P(kpzhi),
#endif
			  const Real* lbd, ARLIM_P(lbdlo), ARLIM_P(lbdhi),
			  const Real* pc, ARLIM_P(pclo), ARLIM_P(pchi),
			  const Real* P, ARLIM_P(plo), ARLIM_P(phi),
			  const Real* K, ARLIM_P(klo), ARLIM_P(khi),
			  const Real* kr,  ARLIM_P(krlo),  ARLIM_P(krhi),
			  const int* nkr,
			  const Real* cf,  ARLIM_P(cflo),  ARLIM_P(cfhi),
			  const int* ncf,
			  const int* lo, const int* hi,
			  const int* domlo, const int* domhi,
			  const Real* dx, const int* bc,
			  const int* inflow_bc_lo, const int* inflow_bc_hi,
			  const Real* deps, const int* do_upwind);

    void FORT_RICHARD_ALPHA (const Real* alpha,  ARLIM_P(alo),  ARLIM_P(ahi),
			     const Real* s,   ARLIM_P(slo),   ARLIM_P(shi),
			     const Real* phi, ARLIM_P(philo), ARLIM_P(phihi),
			     const Real* kp,  ARLIM_P(kplo),  ARLIM_P(kphi),
			     const Real* cf,  ARLIM_P(cflo),  ARLIM_P(cfhi),
			     const int* ncf, const int* lo, const int* hi);

    void FORT_RICHARD_VELBC(const Real* rg, ARLIM_P(rglo), ARLIM_P(rghi),
			    const Real* ux, ARLIM_P(uxlo), ARLIM_P(uxhi),
			    const Real* uy, ARLIM_P(uylo), ARLIM_P(uyhi),
#if (BL_SPACEDIM==3)
			    const Real* uz, ARLIM_P(uzlo), ARLIM_P(uzhi),
#endif
			    const int*  lo,    const int*  hi,
			    const int*  domlo, const int*  domhi,
			    const Real* dx,
			    const int* inflow_bc_lo, const int* inflow_bc_hi,
			    const Real* dt);

  void FORT_AVG_UMAC (const Real* umac, ARLIM_P(umlo), ARLIM_P(umhi),
		      const Real* vmac, ARLIM_P(vmlo), ARLIM_P(vmhi),
#if (BL_SPACEDIM == 3)
		      const Real* wmac, ARLIM_P(wmlo), ARLIM_P(wmhi),
#endif
		      const Real* u   , ARLIM_P( ulo), ARLIM_P( uhi),
		      const int* lo, const int* hi, const int* upwind_flag);

  void FORT_CPY_UMAC (const Real* umac, ARLIM_P(umlo), ARLIM_P(umhi),
		      const Real* vmac, ARLIM_P(vmlo), ARLIM_P(vmhi),
#if (BL_SPACEDIM == 3)
		      const Real* wmac, ARLIM_P(wmlo), ARLIM_P(wmhi),
#endif
		      const Real* u   , ARLIM_P( ulo), ARLIM_P( uhi),
		      const int* lo, const int* hi, const int* ishift);

  void FORT_LINEAR_PRESSURE(const int* lo, const int* hi,
                            Real* p, ARLIM_P(p), ARLIM_P(p), const int* ncomps,
                            const Real* dx, const Real* problo, const Real* probhi, const Real* ref_val,
                            const Real* ref_loc, const Real* gradp);

  void FORT_HYD_HEAD(const int* lo, const int* hi,
                     const Real* pin_hout, ARLIM_P(plo), ARLIM_P(phi),
                     const Real* rhog, const int* gdir, const Real* dx, const Real* plo,
                     const int* ncomps);

  void FORT_STEADYSTATE(Real* rho, ARLIM_P(rlo), ARLIM_P(rhi),
			const Real* rhoval, const Real* muval, const int* ncomps,
			const Real* kappa, ARLIM_P(klo), ARLIM_P(khi),
			const int*  mid,   ARLIM_P(mlo), ARLIM_P(mhi),
			const int* rmID, const Real* time,
			const Real* inflow_vel, const int* nc, const Real* grav,
			const int* lo, const int* hi);

  void FORT_TENSORDIFF(const int* lo, const int* hi,
                       const Real* a, ARLIM_P(alo), ARLIM_P(ahi),
                       const Real* Dcc, ARLIM_P(Dcclo), ARLIM_P(Dcchi),

                       const Real* u, ARLIM_P(ulo), ARLIM_P(uhi),
                       const Real* bx, ARLIM_P(bxlo), ARLIM_P(bxhi),
                       const Real* b1x, ARLIM_P(b1xlo), ARLIM_P(b1xhi),

                       const Real* v, ARLIM_P(vlo), ARLIM_P(vhi),
                       const Real* by, ARLIM_P(bylo), ARLIM_P(byhi),
                       const Real* b1y, ARLIM_P(b1ylo), ARLIM_P(b1yhi),

#if BL_SPACEDIM>2
                       const Real* w, ARLIM_P(wlo), ARLIM_P(whi),
                       const Real* bz, ARLIM_P(bzlo), ARLIM_P(bzhi),
                       const Real* b1z, ARLIM_P(b1zlo), ARLIM_P(b1zhi),
#endif
                       const int* bc);

  void FORT_VISCEXTRAP (Real* visc, ARLIM_P(visclo), ARLIM_P(isclo),
                        const int* lo, const int* hi,
                        const int* comp);

  void FORT_INITKEDGE (const Real* cdat,ARLIM_P(c_lo), ARLIM_P(c_hi),
                       Real*       edat,ARLIM_P(e_lo), ARLIM_P(e_hi),
                       const int* lo, const int* hi, const int* dir);
}

#endif
#endif /*_PorousMedia_F_H_*/
