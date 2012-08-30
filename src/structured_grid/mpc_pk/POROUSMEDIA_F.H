#ifndef _POROUSMEDIA_F_H_
#define _POROUSMEDIA_F_H_

#ifdef BL_LANG_FORT
#    define FORT_FILCC             filcc
#    define FORT_AVGDOWN           avgdown
#    define FORT_SUMMASS           summass
#    define FORT_SUMMASSSQ         summasssq
#    define FORT_CEN2EDG           cen2edg
#    define FORT_EDGE_INTERP       edge_interp
#    define FORT_EDGE_AVGDOWN      edge_avgdown
#    define FORT_PC_CF_EDGE_INTERP pc_edge_interp
#    define FORT_CV_EDGE_INTERP    cv_edge_interp
#    define FORT_CVF_EDGE_INTERP   cvf_edge_interp
#    define FORT_REACH             reach
#    define FORT_SPECTEMPVISC      spectempvisc
#    define FORT_CALCEOS           fcalceos
#    define FORT_RECOMP_UPDATE     recomp_update
#    define FORT_CHEM_DUMMY        chem_dummy
#    define FORT_MK_CPL            mkcpl
#    define FORT_MK_INV_CPL        mkinvcpl
#    define FORT_SMOOTH_CPL        smoothcpl
#    define FORT_MK_LAMBDA         mklambda
#    define FORT_MK_DLAMBDA        mkdlambda
#    define FORT_GETDIFFUSE_CPL    getdiffuse_cpl
#    define FORT_GETDIFFUSE_CPL_dp getdiffuse_cpl_dp
#    define FORT_RICHARD_COEF      richard_coef
#    define FORT_RICHARD_AJAC      richard_ajac
#    define FORT_RICHARD_NJAC      richard_njac
#    define FORT_RICHARD_NJAC2     richard_njac2
#    define FORT_RICHARD_ALPHA     richard_alpha
#    define FORT_RICHARD_VELBC     richard_velbc
#    define FORT_DIV_UMAC          fdiv_umac
#    define FORT_UPHASE            fuphase
#    define FORT_UPHASE_P          fuphasep
#    define FORT_AVG_UMAC          avgumac
#    define FORT_CPY_UMAC          cpyumac
#    define FORT_HYDRO             hydro
#    define FORT_HYDRO_PRESSURE    hydro_pressure
#    define FORT_ROCKHOLD          rockhold
#    define FORT_STEADYSTATE       steadystate
#    define FORT_STEADYSTATE_FAB   steadystatefab
#    define FORT_PATCH_GHOST       patch_ghost
#else
#  if defined(BL_FORT_USE_UPPERCASE)
#    define FORT_FILCC             FILCC
#    define FORT_AVGDOWN           AVGDOWN
#    define FORT_SUMMASS           SUMMASS
#    define FORT_SUMMASSSQ         SUMMASSSQ
#    define FORT_CEN2EDG           CEN2EDG
#    define FORT_EDGE_INTERP       EDGE_INTERP
#    define FORT_EDGE_AVGDOWN      EDGE_AVGDOWN
#    define FORT_PC_CF_EDGE_INTERP PC_EDGE_INTERP
#    define FORT_CV_EDGE_INTERP    CV_EDGE_INTERP
#    define FORT_CVF_EDGE_INTERP   CVF_EDGE_INTERP
#    define FORT_REACH             REACH
#    define FORT_SPECTEMPVISC      SPECTEMPVISC
#    define FORT_CALCEOS           FCALCEOS
#    define FORT_RECOMP_UPDATE     RECOMP_UPDATE
#    define FORT_CHEM_DUMMY        CHEM_DUMMY
#    define FORT_MK_CPL            MKCPL
#    define FORT_MK_INV_CPL        MKINVCPL
#    define FORT_SMOOTH_CPL        smoothcpl
#    define FORT_MK_LAMBDA         MKLAMBDA
#    define FORT_MK_DLAMBDA        MKDLAMBDA
#    define FORT_GETDIFFUSE_CPL    GETDIFFUSE_CPL
#    define FORT_GETDIFFUSE_CPL_dp GETDIFFUSE_CPL_DP
#    define FORT_RICHARD_COEF      RICHARD_COEF
#    define FORT_RICHARD_AJAC      RICHARD_AJAC
#    define FORT_RICHARD_NJAC      RICHARD_NJAC
#    define FORT_RICHARD_NJAC2     RICHARD_NJAC2
#    define FORT_RICHARD_ALPHA     RICHARD_ALPHA
#    define FORT_RICHARD_VELBC     RICHARD_VELBC
#    define FORT_DIV_UMAC          FDIV_UMAC
#    define FORT_UPHASE            FUPHASE
#    define FORT_UPHASE_P          FUPHASEP
#    define FORT_AVG_UMAC          AVGUMAC
#    define FORT_CPY_UMAC          CPYUMAC
#    define FORT_HYDRO             HYDRO
#    define FORT_HYDRO_PRESSURE    HYDRO_PRESSURE
#    define FORT_ROCKHOLD          ROCKHOLD
#    define FORT_STEADYSTATE       STEADYSTATE
#    define FORT_STEADYSTATE_FAB   STEADYSTATEFAB
#    define FORT_PATCH_GHOST       PATCH_GHOST
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_FILCC             filcc
#    define FORT_AVGDOWN           avgdown
#    define FORT_SUMMASS           summass
#    define FORT_SUMMASSSQ         summasssq
#    define FORT_CEN2EDG           cen2edg
#    define FORT_EDGE_INTERP       edge_interp
#    define FORT_EDGE_AVGDOWN      edge_avgdown
#    define FORT_PC_CF_EDGE_INTERP pc_edge_interp
#    define FORT_CV_EDGE_INTERP    cv_edge_interp
#    define FORT_CVF_EDGE_INTERP   cvf_edge_interp
#    define FORT_REACH             reach
#    define FORT_SPECTEMPVISC      spectempvisc
#    define FORT_CALCEOS           fcalceos
#    define FORT_RECOMP_UPDATE     recomp_update
#    define FORT_CHEM_DUMMY        chem_dummy
#    define FORT_MK_CPL            mkcpl
#    define FORT_MK_INV_CPL        mkinvcpl
#    define FORT_SMOOTH_CPL        smoothcpl
#    define FORT_MK_LAMBDA         mklambda
#    define FORT_MK_DLAMBDA        mkdlambda
#    define FORT_GETDIFFUSE_CPL    getdiffuse_cpl
#    define FORT_GETDIFFUSE_CPL_dp getdiffuse_cpl_dp
#    define FORT_RICHARD_COEF      richard_coef
#    define FORT_RICHARD_AJAC      richard_ajac
#    define FORT_RICHARD_NJAC      richard_njac
#    define FORT_RICHARD_NJAC2     richard_njac2
#    define FORT_RICHARD_ALPHA     richard_alpha
#    define FORT_RICHARD_VELBC     richard_velbc
#    define FORT_DIV_UMAC          fdiv_umac
#    define FORT_UPHASE            fuphase
#    define FORT_UPHASE_P          fuphasep
#    define FORT_AVG_UMAC          avgumac
#    define FORT_CPY_UMAC          cpyumac
#    define FORT_HYDRO             hydro
#    define FORT_HYDRO_PRESSURE    hydro_pressure
#    define FORT_ROCKHOLD          rockhold
#    define FORT_STEADYSTATE       steadystate
#    define FORT_STEADYSTATE_FAB   steadystatefab
#    define FORT_PATCH_GHOST       patch_ghost
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_FILCC             filcc_
#    define FORT_AVGDOWN           avgdown_
#    define FORT_SUMMASS           summass_
#    define FORT_SUMMASSSQ         summasssq_
#    define FORT_CEN2EDG           cen2edg_
#    define FORT_EDGE_INTERP       edge_interp_
#    define FORT_EDGE_AVGDOWN      edge_avgdown_
#    define FORT_PC_CF_EDGE_INTERP pc_edge_interp_
#    define FORT_CV_EDGE_INTERP    cv_edge_interp_
#    define FORT_CVF_EDGE_INTERP   cvf_edge_interp_
#    define FORT_REACH             reach_
#    define FORT_SPECTEMPVISC      spectempvisc_
#    define FORT_CALCEOS           fcalceos_
#    define FORT_RECOMP_UPDATE     recomp_update_
#    define FORT_CHEM_DUMMY        chem_dummy_
#    define FORT_MK_CPL            mkcpl_
#    define FORT_MK_INV_CPL        mkinvcpl_
#    define FORT_SMOOTH_CPL        smoothcpl_
#    define FORT_MK_LAMBDA         mklambda_
#    define FORT_MK_DLAMBDA        mkdlambda_
#    define FORT_GETDIFFUSE_CPL    getdiffuse_cpl_
#    define FORT_GETDIFFUSE_CPL_dp getdiffuse_cpl_dp_
#    define FORT_RICHARD_COEF      richard_coef_
#    define FORT_RICHARD_AJAC      richard_ajac_
#    define FORT_RICHARD_NJAC      richard_njac_
#    define FORT_RICHARD_NJAC2     richard_njac2_
#    define FORT_RICHARD_ALPHA     richard_alpha_
#    define FORT_RICHARD_VELBC     richard_velbc_
#    define FORT_DIV_UMAC          fdiv_umac_
#    define FORT_UPHASE            fuphase_
#    define FORT_UPHASE_P          fuphasep_
#    define FORT_AVG_UMAC          avgumac_
#    define FORT_CPY_UMAC          cpyumac_
#    define FORT_HYDRO             hydro_
#    define FORT_HYDRO_PRESSURE    hydro_pressure_
#    define FORT_ROCKHOLD          rockhold_
#    define FORT_STEADYSTATE       steadystate_
#    define FORT_STEADYSTATE_FAB   steadystatefab_
#    define FORT_PATCH_GHOST       patch_ghost_
#  endif

#include <ArrayLim.H>
#include "PROB_PM_F.H"

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

    void FORT_SUMMASSSQ (const Real* dat, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                         ARLIM_P(grid_lo), ARLIM_P(grid_hi),
                         const Real* delta, Real* mass);

    void FORT_CEN2EDG(const int* lo, const int* hi,
                      ARLIM_P(clo), ARLIM_P(chi), const Real* cfabdat,
                      ARLIM_P(elo), ARLIM_P(ehi), Real* efabdat,
                      const int* nc, const int* dir, const int* isharm);

    void FORT_EDGE_INTERP(const int* flo, const int* fhi,
                          const int* nc, const int* refRatio, const int* dir,
                          Real*       fine, ARLIM_P(flo), ARLIM_P(fhi));

    void FORT_EDGE_AVGDOWN(const int* lo, const int* hi,
			   const int* nc, const int* refRatio, const int* dir,
			   const Real* crse, ARLIM_P(clo), ARLIM_P(chi),
			   Real*       fine, ARLIM_P(flo), ARLIM_P(fhi));
    
    void FORT_PC_CF_EDGE_INTERP(const int* lo, const int* hi,
                                const int* nc, const int* refRatio, const int* dir,
                                const Real* crse, ARLIM_P(clo), ARLIM_P(chi),
                                Real*       fine, ARLIM_P(flo), ARLIM_P(fhi));

    void FORT_CV_EDGE_INTERP(const int* lo, const int* hi,
			     const int* domlo, const int* domhi,
			     const int* nc, const int* refRatio, 
			     const int* dir, const Real* dx,
			     const Real* crse, ARLIM_P(clo), ARLIM_P(chi),
			     Real*       fine, ARLIM_P(flo), ARLIM_P(fhi));

    void FORT_CVF_EDGE_INTERP(const int* clo, const int* chi,
			      const int* nc, const int* refRatio,
			      const Real* cdx, const Real* fdx,
			      const Real* crsex, ARLIM_P(cxlo), ARLIM_P(cxhi),
			      const Real* crsey, ARLIM_P(cylo), ARLIM_P(cyhi),
			      const Real* finex, ARLIM_P(fxlo), ARLIM_P(fxhi),
			      const Real* finey, ARLIM_P(fylo), ARLIM_P(fyhi));
  
    void FORT_SPECTEMPVISC(const int* lo, const int* hi,
			   const Real* S, ARLIM_P(slo), ARLIM_P(shi),
			   const Real* B, ARLIM_P(blo), ARLIM_P(bhi),
			   const Real* P, ARLIM_P(plo), ARLIM_P(phi),
			   const Real* visc_coef, const int* ncomp, 
			   const int* diffuse);

    //
    // subroutines related to parabolic pressure solve
    //
    void FORT_CALCEOS(const Real* coef, ARLIM_P(coeflo), ARLIM_P(coefhi),
		      const Real* s,    ARLIM_P(slo),    ARLIM_P(shi),
		      const Real* rho,  ARLIM_P(rholo),  ARLIM_P(rhohi),
		      const Real* phi,  ARLIM_P(philo),  ARLIM_P(phihi),
		      const int* ncomp, const int* lo, const int* hi,
		      const Real* dt);

    void FORT_RECOMP_UPDATE(const int* lo, const int* hi,
			    const Real* updat, ARLIM_P(uplo), ARLIM_P(uphi),
			    const Real* phidat, ARLIM_P(philo), ARLIM_P(phihi),
			    const Real* xflux, ARLIM_P(fxlo), ARLIM_P(fxhi),
			    const Real* yflux, ARLIM_P(fylo), ARLIM_P(fyhi),
#if (BL_SPACEDIM==3)
			    const Real* zflux, ARLIM_P(fzlo), ARLIM_P(fzhi),
#endif
			    const Real* vol,   ARLIM_P(vlo),  ARLIM_P(vhi),
			    int* nc);

   void FORT_CHEM_DUMMY(const Real* s, ARLIM_P(slo), ARLIM_P(shi), 
			Real* dt, int* nc);

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
   
   void FORT_SMOOTH_CPL (const Real* pc,  ARLIM_P(pclo),  ARLIM_P(pchi),
			 const Real* cf,  ARLIM_P(cflo),  ARLIM_P(cfhi),
			 const int* ncf, const int* lo, const int* hi);  

  void FORT_MK_LAMBDA (const Real* lbd, ARLIM_P(lbdlo), ARLIM_P(lbdhi),
		       const Real* s,   ARLIM_P(slo),   ARLIM_P(shi),
		       const Real* kr,  ARLIM_P(krlo),  ARLIM_P(krhi),
		       const int* nkr);

  void FORT_MK_DLAMBDA (const Real* lbd,  ARLIM_P(lbdlo),  ARLIM_P(lbdhi),
			const Real* s,   ARLIM_P(slo),   ARLIM_P(shi),
			const Real* kr,  ARLIM_P(krlo),  ARLIM_P(krhi),
			const int* nkr);

   void FORT_GETDIFFUSE_CPL(const Real* lbd, ARLIM_P(lbdlo), ARLIM_P(lbdhi),
			    const Real* dfx, ARLIM_P(dfxlo), ARLIM_P(dfxhi),
			    const Real* dfy, ARLIM_P(dfylo), ARLIM_P(dfyhi),
#if (BL_SPACEDIM==3)
			    const Real* dfz, ARLIM_P(dfzlo), ARLIM_P(dfzhi),
#endif
			    const int* lo, const int* hi,
			    const int* domlo, const int* domhi,
			    const int* bc, const int* ncomp);

  void FORT_GETDIFFUSE_CPL_dp(const Real* S, const Real* h_cc,
			      ARLIM_P(slo), ARLIM_P(shi),
			      const Real* lbd, ARLIM_P(lbdlo), ARLIM_P(lbdhi),
			      const Real* dfx, ARLIM_P(dfxlo), ARLIM_P(dfxhi),
			      const Real* dfy, ARLIM_P(dfylo), ARLIM_P(dfyhi),
#if (BL_SPACEDIM==3)
			      const Real* dfz, ARLIM_P(dfzlo), ARLIM_P(dfzhi),
#endif
			      const Real* P, ARLIM_P(plo), ARLIM_P(phi), 
			      const Real* K, ARLIM_P(klo), ARLIM_P(khi),
			      const Real* cf,  ARLIM_P(cflo),  ARLIM_P(cfhi),
			      const int* ncf,
			      const int* lo, const int* hi,
			      const int* domlo, const int* domhi,
			      const int* bc, const int* ncomp);

   void FORT_RICHARD_COEF(const Real* lbd, ARLIM_P(lbdlo), ARLIM_P(lbdhi),
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
			  const int* lo, const int* hi,
			  const int* domlo, const int* domhi,const int* bc, 
			  const int* nc, const int* do_upwind);

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

    void FORT_DIV_UMAC(const Real* sdat, ARLIM_P(slo), ARLIM_P(shi),
		       const Real* xflux, ARLIM_P(fxlo), ARLIM_P(fxhi),
		       const Real* yflux, ARLIM_P(fylo), ARLIM_P(fyhi),
#if (BL_SPACEDIM==3)
		       const Real* zflux, ARLIM_P(fzlo), ARLIM_P(fzhi),
#endif
		       const int* lo, const int* hi, const Real* dx);
  
    void FORT_UPHASE (const Real* uphase, ARLIM_P(uplo), ARLIM_P(uphi),
		      const Real* vphase, ARLIM_P(vplo), ARLIM_P(vphi),
#if (BL_SPACEDIM == 3)
		      const Real* wphase, ARLIM_P(wplo), ARLIM_P(wphi),
#endif
		      const Real* umac, ARLIM_P(umlo), ARLIM_P(umhi),
		      const Real* vmac, ARLIM_P(vmlo), ARLIM_P(vmhi),
#if (BL_SPACEDIM == 3)
		      const Real* wmac, ARLIM_P(wmlo), ARLIM_P(wmhi),
#endif
		      const Real* pc, ARLIM_P(pclo), ARLIM_P(pchi),
		      const Real* lbd, ARLIM_P(lbdlo), ARLIM_P(lbdhi),
		      const Real* kpx, ARLIM_P(kpxlo), ARLIM_P(kpxhi),
		      const Real* kpy, ARLIM_P(kpylo), ARLIM_P(kpyhi),
#if (BL_SPACEDIM == 3)
		      const Real* kpz, ARLIM_P(kpzlo), ARLIM_P(kpzhi),
#endif
		      const int* bc, 
		      const int* dom_lo, const int* dom_hi, 
		      const Real* dx, const int* lo, const int* hi);   

    void FORT_UPHASE_P (const Real* uphase, ARLIM_P(uplo), ARLIM_P(uphi),
			const Real* vphase, ARLIM_P(vplo), ARLIM_P(vphi),
#if (BL_SPACEDIM == 3)
			const Real* wphase, ARLIM_P(wplo), ARLIM_P(wphi),
#endif
			const Real* lbd, ARLIM_P(lbdlo), ARLIM_P(lbdhi),
			const Real* pc, ARLIM_P(pclo), ARLIM_P(pchi),
			const Real* kpx, ARLIM_P(kpxlo), ARLIM_P(kpxhi),
			const Real* kpy, ARLIM_P(kpylo), ARLIM_P(kpyhi),
#if (BL_SPACEDIM == 3)
			const Real* kpz, ARLIM_P(kpzlo), ARLIM_P(kpzhi),
#endif
			const int* lo, const int* hi,
			const int* domlo, const int* domhi, 
			const Real* dx, const int* bc);   
  
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

  void FORT_HYDRO(Real* rho, ARLIM_P(rlo), ARLIM_P(rhi), 
		  const Real* rhoval, const int* ncomps,
		  const Real* cplcoef, ARLIM_P(clo), ARLIM_P(chi),
		  const int* ncplcoef,const Real* dx, 
		  const Real* wt_loc, const Real* grav);

  void FORT_HYDRO_PRESSURE(Real* p, ARLIM_P(p), ARLIM_P(p), 
			   const Real* rhoval, const int* ncomps,
			   const Real* dx, 
			   const Real* wt_loc, const Real* grav);

  void FORT_ROCKHOLD(Real* rho, ARLIM_P(rlo), ARLIM_P(rhi), 
		     const Real* rhoval, const int* ncomps,
		     const Real* depth, const Real* pressure, const int* nz,
		     const Real* cplcoef, ARLIM_P(clo), ARLIM_P(chi),
		     const int* ncplcoef,const Real* dx,  const Real* grav, const Real* probhi);

  void FORT_STEADYSTATE(Real* rho, ARLIM_P(rlo), ARLIM_P(rhi), 
			const Real* rhoval, const Real* muval, const int* ncomps,
			const Real* kappa, ARLIM_P(klo), ARLIM_P(khi),
			const Real* krcoef, ARLIM_P(clo), ARLIM_P(chi),
			const int* nkrcoef,const Real* dx, 
			const Real* inflow_vel, const int* nc, const Real* grav);  

  void FORT_STEADYSTATE_FAB(Real* rho, ARLIM_P(rlo), ARLIM_P(rhi), 
			    const Real* rhoval, const Real* muval, const int* ncomps,
			    const Real* kappa, ARLIM_P(klo), ARLIM_P(khi),
			    const Real* krcoef, ARLIM_P(clo), ARLIM_P(chi), const int* nkrcoef,
			    const Real* vel, ARLIM_P(vlo), ARLIM_P(vhi),
			    const Real* dx, const int* nc, const Real* grav);
  
  void FORT_PATCH_GHOST(const Real* p, ARLIM_P(p_lo), ARLIM_P(p_hi),
			const int* ncomp, int* face, 
			const int* domlo, const int* domhi);

}
    
#endif
#endif /*_POROUSMEDIA_F_H_*/
