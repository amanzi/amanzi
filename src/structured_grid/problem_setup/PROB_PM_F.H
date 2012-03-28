#ifndef _PROB_PM_F_H_
#define _PROB_PM_F_H_

#ifdef BL_LANG_FORT
#    define FORT_INITDATA    initdata
#    define FORT_INIT_TRACER inittracer
#    define FORT_RECINITDATA recinitdata
#    define FORT_INITPARAMS  initparams
#    define FORT_TCRPARAMS   tcrparams
#    define FORT_ONE_N_FILL  onenfill
#    define FORT_ALL_N_FILL  allnfill
#    define FORT_ALL_T_FILL  alltfill
#    define FORT_REC_N_FILL  recnfill
#    define FORT_ADVFILL     advfill
#    define FORT_XVELFILL    xvelfill
#    define FORT_YVELFILL    yvelfill
#    define FORT_ZVELFILL    zvelfill
#    define FORT_PRESFILL    presfill
#    define FORT_ENTHFILL    enthfill
#    define FORT_DIVUFILL    divufill
#    define FORT_ADVERROR    adverror
#    define FORT_VALGTERROR  valgterr
#    define FORT_DIFFGTERROR diffgterr
#    define FORT_VALLTERROR  vallterr
#    define FORT_KPERROR     kperror
#    define FORT_INIT_UMAC   initumac
#    define FORT_CHEM        chem
#    define FORT_USER_KAPPA  userkappa
#    define FORT_MK_FORCE    mkforce
#    define FORT_MK_TFORCE   mktforce
#    define FORT_FIND_INV_FLUX findinvflux
#else
#  if  defined(BL_FORT_USE_UPPERCASE)
#    define FORT_INITDATA    INITDATA
#    define FORT_INIT_TRACER INITTRACER
#    define FORT_RECINITDATA RECINITDATA
#    define FORT_INITPARAMS  INITPARAMS
#    define FORT_TCRPARAMS   TCRPARAMS
#    define FORT_MK_MACCOEF  MKMACCOEF
#    define FORT_MK_MACDIVU  MKMACDIVU
#    define FORT_ONE_N_FILL  ONENFILL
#    define FORT_ALL_N_FILL  ALLNFILL
#    define FORT_ALL_T_FILL  ALLTFILL
#    define FORT_REC_N_FILL  RECNFILL
#    define FORT_ADVFILL     ADVFILL
#    define FORT_XVELFILL    XVELFILL
#    define FORT_YVELFILL    YVELFILL
#    define FORT_ZVELFILL    ZVELFILL
#    define FORT_PRESFILL    PRESFILL
#    define FORT_ENTHFILL    ENTHFILL
#    define FORT_DIVUFILL    DIVUFILL
#    define FORT_ADVERROR    ADVERROR
#    define FORT_VALGTERROR  VALGTERR
#    define FORT_VALLTERROR  VALLTERR
#    define FORT_DIFFGTERROR DIFFGTERR
#    define FORT_KPERROR     KPERROR
#    define FORT_INIT_UMAC   INITUMAC
#    define FORT_CHEM        CHEM
#    define FORT_USER_KAPPA  USERKAPPA
#    define FORT_MK_FORCE    MKFORCE
#    define FORT_MK_TFORCE   MKTFORCE
#    define FORT_MAKEFORCE   MAKEFORCE
#    define FORT_FIND_INV_FLUX FINDINVFLUX
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_INITDATA    initdata
#    define FORT_INIT_TRACER inittracer
#    define FORT_RECINITDATA recinitdata
#    define FORT_INITPARAMS  initparams
#    define FORT_TCRPARAMS   tcrparams
#    define FORT_MK_MACCOEF  mkmaccoef
#    define FORT_MK_MACDIVU  mkmacdivu
#    define FORT_ONE_N_FILL  onenfill
#    define FORT_ALL_N_FILL  allnfill
#    define FORT_ALL_T_FILL  alltfill
#    define FORT_REC_N_FILL  recnfill
#    define FORT_ADVFILL     advfill
#    define FORT_XVELFILL    xvelfill
#    define FORT_YVELFILL    yvelfill
#    define FORT_ZVELFILL    zvelfill
#    define FORT_PRESFILL    presfill
#    define FORT_ENTHFILL    enthfill
#    define FORT_DIVUFILL    divufill
#    define FORT_ADVERROR    adverror
#    define FORT_VALGTERROR  valgterr
#    define FORT_VALLTERROR  vallterr
#    define FORT_DIFFGTERROR diffgterr
#    define FORT_KPERROR     kperror
#    define FORT_INIT_UMAC   initumac
#    define FORT_MK_CPL      mkcpl
#    define FORT_CHEM        chem
#    define FORT_USER_KAPPA  userkappa
#    define FORT_MK_FORCE    mkforce
#    define FORT_MK_TFORCE   mktforce
#    define FORT_FIND_INV_FLUX findinvflux
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_INITDATA    initdata_
#    define FORT_INIT_TRACER inittracer_
#    define FORT_RECINITDATA recinitdata_
#    define FORT_INITPARAMS  initparams_
#    define FORT_TCRPARAMS   tcrparams_
#    define FORT_MK_MACCOEF  mkmaccoef_
#    define FORT_MK_MACDIVU  mkmacdivu_
#    define FORT_ONE_N_FILL  onenfill_
#    define FORT_ALL_N_FILL  allnfill_
#    define FORT_ALL_T_FILL  alltfill_
#    define FORT_REC_N_FILL  recnfill_
#    define FORT_ADVFILL     advfill_
#    define FORT_XVELFILL    xvelfill_
#    define FORT_YVELFILL    yvelfill_
#    define FORT_ZVELFILL    zvelfill_
#    define FORT_PRESFILL    presfill_
#    define FORT_ENTHFILL    enthfill_
#    define FORT_DIVUFILL    divufill_
#    define FORT_ADVERROR    adverror_
#    define FORT_VALGTERROR  valgterr_
#    define FORT_VALLTERROR  vallterr_
#    define FORT_DIFFGTERROR diffgterr_
#    define FORT_KPERROR     kperror_
#    define FORT_INIT_UMAC   initumac_
#    define FORT_CHEM        chem_
#    define FORT_USER_KAPPA  userkappa_
#    define FORT_MK_FORCE    mkforce_
#    define FORT_MK_TFORCE   mktforce_
#    define FORT_FIND_INV_FLUX findinvflux_
#  endif

#include "ArrayLim.H"

  extern "C" {

     void FORT_INITDATA (const int* level, const Real* time,
                         const Real* scal, ARLIM_P(s_lo), ARLIM_P(s_hi),
			 const Real* rhoinit, const int* nscal,
			 const Real* dx);

     void FORT_INIT_TRACER (const int* level, const Real* time,
			    const Real* scal, ARLIM_P(s_lo), ARLIM_P(s_hi),
			    const Real* rhoinit, const int* nscal,
			    const int* ntracer, const Real* dx);

     void FORT_RECINITDATA (const int* level, const Real* time,
			    const int* lo, const int* hi,
			    const Real* scal, ARLIM_P(s_lo), ARLIM_P(s_hi),
			    const Real* dx,const Real* xlo, const Real* xhi);

     void FORT_INITPARAMS (int* ncomps, int* nphases, int* model_in, Real* rhoval, 
			   Real* muval, int* phaseidx,Real* gravity);
    
     void FORT_TCRPARAMS (int* ntracers);

     void FORT_COREREACT_PARAMS(int*nreact, Real* rhosec); 
    
     void FORT_ONE_N_FILL (Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
			   const int* dom_lo, const int* dom_hi,
			   const Real* dx, const Real* grd_lo,
			   const Real* time, const int* bc);
    
     void FORT_ALL_N_FILL (Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
			   const int* dom_lo, const int* dom_hi,
			   const Real* dx, const Real* grd_lo,
			   const Real* time, const int* bc);

     void FORT_ALL_T_FILL (Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
			   const int* dom_lo, const int* dom_hi,
			   const Real* dx, const Real* grd_lo,
			   const Real* time, const int* bc);

     void FORT_ADVFILL  (Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi), 
                         const int* dom_lo, const int* dom_hi,
                         const Real* dx, const Real* grd_lo,
                         const Real* time, const int* bc);

     void FORT_XVELFILL (Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                         const int* dom_lo, const int* dom_hi,
                         const Real* dx, const Real* grd_lo,
                         const Real* time, const int* bc);

     void FORT_YVELFILL (Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                         const int* dom_lo, const int* dom_hi,
                         const Real* dx, const Real* grd_lo,
                         const Real* time, const int* bc);

#if (BL_SPACEDIM == 3)
     void FORT_ZVELFILL (Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                         const int* dom_lo, const int* dom_hi,
                         const Real* dx, const Real* grd_lo,
                         const Real* time, const int* bc);
#endif
     void FORT_PRESFILL (Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                         const int* dom_lo, const int* dom_hi,
                         const Real* dx, const Real* grd_lo,
                         const Real* time, const int* bc);

     void FORT_ENTHFILL (Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi), 
                         const int* dom_lo, const int* dom_hi,
                         const Real* dx, const Real* grd_lo,
                         const Real* time, const int* bc);

     void FORT_DIVUFILL (Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                         const int* dom_lo, const int* dom_hi,
                         const Real* dx, const Real* grd_lo,
                         const Real* time, const int* bc);
 
     void FORT_ADVERROR  (int* tag, ARLIM_P(tlo), ARLIM_P(thi),
                          const int* tagval, const int* clearval,
                          Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                          const int* lo, const int* hi, const int* nvar,
                          const int* domain_lo, const int* domain_hi,
                          const Real* dx, const Real* xlo,
                          const Real* prob_lo, const Real* time,
                          const int* level);
    
     void FORT_VALGTERROR(int* tag, ARLIM_P(tlo), ARLIM_P(thi),
                          const int* tagval, const int* clearval,
                          const Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                          const Real* mask, ARLIM_P(mask_lo), ARLIM_P(mask_hi),
                          const int* lo, const int* hi, const int* nvar,
                          const int* domain_lo, const int* domain_hi,
                          const Real* dx, const Real* xlo,
                          const Real* prob_lo, const Real* time,
                          const int* level, const Real* value);
    
     void FORT_VALLTERROR(int* tag, ARLIM_P(tlo), ARLIM_P(thi),
                          const int* tagval, const int* clearval,
                          const Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                          const Real* mask, ARLIM_P(mask_lo), ARLIM_P(mask_hi),
                          const int* lo, const int* hi, const int* nvar,
                          const int* domain_lo, const int* domain_hi,
                          const Real* dx, const Real* xlo,
                          const Real* prob_lo, const Real* time,
                          const int* level, const Real* value);
    
     void FORT_DIFFGTERROR(int* tag, ARLIM_P(tlo), ARLIM_P(thi),
                           const int* tagval, const int* clearval,
                           const Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                           const Real* mask, ARLIM_P(mask_lo), ARLIM_P(mask_hi),
                           const int* lo, const int* hi, const int* nvar,
                           const int* domain_lo, const int* domain_hi,
                           const Real* dx, const Real* xlo,
                           const Real* prob_lo, const Real* time,
                           const int* level, const Real* value);
    
     void FORT_INIT_UMAC (const Real* umac, ARLIM_P(umlo), ARLIM_P(umhi),
                          const Real* vmac, ARLIM_P(vmlo), ARLIM_P(vmhi),
#if (BL_SPACEDIM == 3)
                          const Real* wmac, ARLIM_P(wmlo), ARLIM_P(wmhi),
#endif
			  const Real* pc, ARLIM_P(pclo), ARLIM_P(pchi),
			  const Real* lbd, ARLIM_P(lbdlo), ARLIM_P(lbdhi),
                          const Real* ldx, ARLIM_P(ldxlo), ARLIM_P(ldxhi),
                          const Real* ldy, ARLIM_P(ldylo), ARLIM_P(ldyhi),
#if (BL_SPACEDIM == 3)
                          const Real* ldz, ARLIM_P(ldzlo), ARLIM_P(ldzhi),
#endif
                          const Real* kpx, ARLIM_P(kpxlo), ARLIM_P(kpxhi),
                          const Real* kpy, ARLIM_P(kpylo), ARLIM_P(kpyhi),
#if (BL_SPACEDIM == 3)
                          const Real* kpz, ARLIM_P(kpzlo), ARLIM_P(kpzhi),
#endif
                          const Real* RhoG, ARLIM_P(rglo), ARLIM_P(rghi),
                          const Real* RhoDX, ARLIM_P(rdxlo), ARLIM_P(rdxhi),
                          const Real* RhoDY, ARLIM_P(rdylo), ARLIM_P(rdyhi),
#if (BL_SPACEDIM == 3)
                          const Real* RhoDZ, ARLIM_P(rdzlo), ARLIM_P(rdzhi),
#endif
                          const Real* n, ARLIM_P(n_lo), ARLIM_P(n_hi),
                          const Real* phi, ARLIM_P(p_lo), ARLIM_P(p_hi), 
                          const Real* const_diff_coef,
			  const int* bc, const int* pbc,
                          const int* dom_lo, const int* dom_hi, 
			  const Real* dx, const int* lo, const int* hi,
			  const Real* wt_lo, const Real* wt_hi,
			  const int* inflow_bc_lo, const int* inflow_bc_hi);

    void FORT_USER_KAPPA (const Real* kp_dat, ARLIM_P(kp_lo), ARLIM_P(kp_hi), 
			  const Real* kappaval, const int* nkpval, 
			  const Real* dkappa, const Real* dx, 
			  const Real* problo, const Real* probhi);

    void FORT_KPERROR (const int* tagdat,   ARLIM_P(tg_lo), ARLIM_P(tg_hi), 
		       const Real* kpdat,   ARLIM_P(kp_lo), ARLIM_P(kp_hi),
		       const int* set, const Real* val, const Real* dx,
		       const Real* problo, const Real* probhi,
		       const int* lo, const int* hi, 
		       const int* domlo, const int* domhi,const int* level);

    void FORT_MK_FORCE (const Real* time, 
			const Real* fdat, ARLIM_P(fd_lo), ARLIM_P(fd_hi),
			const Real* dx, const Real* problo, const Real* probhi,
			const int* scomp, const int * ncomp, 
			const int* do_rho_scale); 

    void FORT_MK_TFORCE (const Real* time, 
			 const Real* fdat, ARLIM_P(fd_lo), ARLIM_P(fd_hi),
			 const Real* dx, const Real* problo, const Real* probhi,
			 const int* scomp, const int * ncomp, 
			 const int* do_rho_scale);

     void FORT_FIND_INV_FLUX(Real* sol, const Real* vel, const int* nc,
			     const Real* vtot,const Real* gstar,const Real* muval,
			     const int* ncomps, const int* kr_type, const Real* kr_coef);

#ifdef COREREACT
    void FORT_CHEM (const Real* s, ARLIM_P(slo), ARLIM_P(shi), int* nc,
		    const int* lo, const int* hi, Real* dt, 
		    const Real* press, ARLIM_P(p_lo), ARLIM_P(p_hi),
		    const Real* kappa, ARLIM_P(k_lo), ARLIM_P(k_hi),
		    const Real* phi, ARLIM_P(ph_lo), ARLIM_P(ph_hi),
		    const Real* cr, ARLIM_P(c_lo), ARLIM_P(c_hi), 
		    int* nd, int* idx_cr, const Real* mweight);
#endif

  }
#endif
#endif /*_PROB_PM_F_H_*/
