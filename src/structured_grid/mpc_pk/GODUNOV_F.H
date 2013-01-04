#ifndef _GODUNOV_F_H_
#define _GODUNOV_F_H_

#ifdef BL_LANG_FORT
#    define FORT_TEST_UMAC          testumac
#    define FORT_ESTATE_FPU         estatefpu
#    define FORT_ESTATE_LIN         estatelin
#    define FORT_ESTATE_RMN         estatermn
#    define FORT_ESTATE_CPL         estatecpl
#    define FORT_ESTATE_PMR         estatepmr
#    define FORT_ESTATE_TRACER      estatetracer
#    define FORT_SYNC_RMN           syncrmn
#    define FORT_ADV_FORCING        advforcing
#    define FORT_ADV_RMN_FORCING    advrmnforcing
#    define FORT_SYNC_ADV_FORCING   syncadvforcing
#    define FORT_SYNC_RMN_FORCING   syncrmnforcing
#    define FORT_CONSSCALMINMAX     consscalminmax
#    define FORT_CONVSCALMINMAX     convscalminmax
#    define FORT_UPDATE_AOFS_TF     updateaofstf
#    define FORT_UPDATE_AOFS_TRACER updateaofstracer
#    define FORT_SUM_TF_DIVU_VISC   sumtfdivuvisc
#    define FORT_SUM_TF_DIVU        sumtfdivu
#    define FORT_SET_PARAMS         setparams
#    define FORT_GET_EIG            geteig
#    define FORT_EST_EIG            esteig
#    define FORT_EST_EIG_LIN        esteiglin
#    define FORT_EST_EIG_TRC        esteigtrc
#    define FORT_EST_EIG_CPL        esteigcpl
#else
#  if defined(BL_FORT_USE_UPPERCASE)
#    define FORT_TEST_UMAC          TESTUMAC
#    define FORT_ESTATE_FPU         ESTATEFPU
#    define FORT_ESTATE_LIN         estatelin
#    define FORT_ESTATE_RMN         ESTATERMN
#    define FORT_ESTATE_CPL         ESTATECPL
#    define FORT_ESTATE_PMR         ESTATEPMR
#    define FORT_ESTATE_TRACER      ESTATETRACER
#    define FORT_SYNC_RMN           SYNCRMN
#    define FORT_ADV_FORCING        ADVFORCING  
#    define FORT_ADV_RMN_FORCING    ADVRMNFORCING
#    define FORT_SYNC_ADV_FORCING   SYNCADVFORCING  
#    define FORT_SYNC_RMN_FORCING   SYNCRMNFORCING
#    define FORT_CONSSCALMINMAX     CONSSCALMINMAX
#    define FORT_CONVSCALMINMAX     CONVSCALMINMAX
#    define FORT_UPDATE_AOFS_TF     UPDATEAOFSTF 
#    define FORT_UPDATE_AOFS_TRACER UPDATEAOFSTRACER
#    define FORT_SUM_TF_DIVU_VISC   SUMTFDIVUVISC 
#    define FORT_SUM_TF_DIVU        SUMTFDIVU 
#    define FORT_SET_PARAMS         SETPARAMS
#    define FORT_GET_EIG            GETEIG
#    define FORT_EST_EIG            ESTEIG
#    define FORT_EST_EIG_LIN        ESTEIGLIN
#    define FORT_EST_EIG_TRC        ESTEIGTRC
#    define FORT_EST_EIG_CPL        ESTEIGCPL
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_TEST_UMAC          testumac
#    define FORT_ESTATE_FPU         estatefpu
#    define FORT_ESTATE_LIN         estatelin
#    define FORT_ESTATE_RMN         estatermn
#    define FORT_ESTATE_CPL         estatecpl
#    define FORT_ESTATE_PMR         estatepmr
#    define FORT_ESTATE_TRACER      estatetracer
#    define FORT_SYNC_RMN           syncrmn
#    define FORT_ADV_FORCING        advforcing
#    define FORT_ADV_RMN_FORCING    advrmnforcing
#    define FORT_SYNC_ADV_FORCING   syncadvforcing
#    define FORT_SYNC_RMN_FORCING   syncrmnforcing
#    define FORT_CONSSCALMINMAX     consscalminmax
#    define FORT_CONVSCALMINMAX     convscalminmax
#    define FORT_UPDATE_AOFS_TF     updateaofstf
#    define FORT_UPDATE_AOFS_TRACER updateaofstracer
#    define FORT_SUM_TF_DIVU_VISC   sumtfdivuvisc
#    define FORT_SUM_TF_DIVU        sumtfdivu
#    define FORT_SET_PARAMS         setparams
#    define FORT_GET_EIG            geteig
#    define FORT_EST_EIG            esteig
#    define FORT_EST_EIG_LIN        esteiglin
#    define FORT_EST_EIG_TRC        esteigtrc
#    define FORT_EST_EIG_CPL        esteigcpl
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_TEST_UMAC          testumac_
#    define FORT_ESTATE_FPU         estatefpu_
#    define FORT_ESTATE_LIN         estatelin_
#    define FORT_ESTATE_RMN         estatermn_
#    define FORT_ESTATE_CPL         estatecpl_
#    define FORT_ESTATE_PMR         estatepmr_
#    define FORT_ESTATE_TRACER      estatetracer_
#    define FORT_SYNC_RMN           syncrmn_
#    define FORT_ADV_FORCING        advforcing_
#    define FORT_ADV_RMN_FORCING    advrmnforcing_
#    define FORT_SYNC_ADV_FORCING   syncadvforcing_
#    define FORT_SYNC_RMN_FORCING   syncrmnforcing_
#    define FORT_CONSSCALMINMAX     consscalminmax_
#    define FORT_CONVSCALMINMAX     convscalminmax_
#    define FORT_UPDATE_AOFS_TF     updateaofstf_
#    define FORT_UPDATE_AOFS_TRACER updateaofstracer_
#    define FORT_SUM_TF_DIVU_VISC   sumtfdivuvisc_
#    define FORT_SUM_TF_DIVU        sumtfdivu_
#    define FORT_SET_PARAMS         setparams_
#    define FORT_GET_EIG            geteig_
#    define FORT_EST_EIG            esteig_
#    define FORT_EST_EIG_LIN        esteiglin_
#    define FORT_EST_EIG_TRC        esteigtrc_
#    define FORT_EST_EIG_CPL        esteigcpl_
#  endif

#include <ArrayLim.H>

extern "C"
{
   void FORT_TEST_UMAC(const Real* um, ARLIM_P(ulo), ARLIM_P(uhi),
		       const Real* vm, ARLIM_P(vlo), ARLIM_P(vhi),
#if (BL_SPACEDIM == 3)                            
		       const Real* wm, ARLIM_P(wlo), ARLIM_P(whi),
#endif                                              
		       const Real* rh, ARLIM_P(rlo), ARLIM_P(rhi),
		       const int* lo, const int* hi, const Real* dt, 
       		       const Real* dx, Real* cfl, const Real* u_max);

   void FORT_EST_EIG(const Real* s_dat,     ARLIM_P(s_lo),ARLIM_P(s_hi),
		     const Real* l_dat,     const Real* dl_dat, 
		     const Real* d2l_dat,   ARLIM_P(l_lo), ARLIM_P(l_hi),
		     const Real* rphi_dat,  ARLIM_P(rp_lo),ARLIM_P(rp_hi),
		     const Real* kr_dat,    ARLIM_P(kr_lo),ARLIM_P(kr_hi),
		     const int*  nkr,
		     const Real* uedge_dat, ARLIM_P(ue_lo), ARLIM_P(ue_hi),
		     const Real* kappax,    ARLIM_P(kpx_lo),ARLIM_P(kpx_hi),
		     const Real* vedge_dat, ARLIM_P(ve_lo), ARLIM_P(ve_hi),
		     const Real* kappay,    ARLIM_P(kpy_lo),ARLIM_P(kpy_hi),
#if (BL_SPACEDIM == 3)
		     const Real* wedge_dat, ARLIM_P(we_lo), ARLIM_P(we_hi),
		     const Real* kappaz,    ARLIM_P(kpz_lo),ARLIM_P(kpz_hi),
#endif
		     const int* lo, const int* hi, const Real* dx,
		     const int* bc, Real* eigmax);

    void FORT_EST_EIG_LIN(const Real* rphi_dat,  ARLIM_P(rp_lo),ARLIM_P(rp_hi),
			  const Real* uedge_dat, ARLIM_P(ue_lo), ARLIM_P(ue_hi),
			  const Real* vedge_dat, ARLIM_P(ve_lo), ARLIM_P(ve_hi),
#if (BL_SPACEDIM == 3)
			  const Real* wedge_dat, ARLIM_P(we_lo), ARLIM_P(we_hi),
#endif
			  const int* lo, const int* hi,  Real* eigmax);

    void FORT_EST_EIG_TRC(const Real* s_dat,     ARLIM_P(s_lo),ARLIM_P(s_hi),
			  const int* nc,
			  const Real* rphi_dat,  ARLIM_P(rp_lo),ARLIM_P(rp_hi),
			  const Real* uedge_dat, ARLIM_P(ue_lo), ARLIM_P(ue_hi),
			  const Real* vedge_dat, ARLIM_P(ve_lo), ARLIM_P(ve_hi),
#if (BL_SPACEDIM == 3)
			  const Real* wedge_dat, ARLIM_P(we_lo), ARLIM_P(we_hi),
#endif
			  const int* lo, const int* hi,  Real* eigmax);

   void FORT_EST_EIG_CPL(const Real* s_dat,     ARLIM_P(s_lo),ARLIM_P(s_hi),
			 const Real* l_dat,     const Real* dl_dat, 
			 const Real* d2l_dat,   ARLIM_P(l_lo), ARLIM_P(l_hi),
			 const Real* rphi_dat,  ARLIM_P(rp_lo),ARLIM_P(rp_hi),
			 const Real* kr_dat,    ARLIM_P(kr_lo),ARLIM_P(kr_hi),
			 const int*  nkr,
			 const Real* pc_dat,    ARLIM_P(pc_lo), ARLIM_P(pc_hi),
			 const Real* uedge_dat, ARLIM_P(ue_lo), ARLIM_P(ue_hi),
			 const Real* kappax,    ARLIM_P(kpx_lo),ARLIM_P(kpx_hi),
			 const Real* vedge_dat, ARLIM_P(ve_lo), ARLIM_P(ve_hi),
			 const Real* kappay,    ARLIM_P(kpy_lo),ARLIM_P(kpy_hi),
#if (BL_SPACEDIM == 3)
			 const Real* wedge_dat, ARLIM_P(we_lo), ARLIM_P(we_hi),
			 const Real* kappaz,    ARLIM_P(kpz_lo),ARLIM_P(kpz_hi),
#endif
			 const int* lo, const int* hi, const Real* dx,
			 const int* bc, Real* eigmax);
   
   void FORT_ESTATE_FPU(const Real* s_dat, const Real* tfr_dat, const Real* divu_dat,
			ARLIM_P(s_lo),ARLIM_P(s_hi),
			const Real* xlo_dat, 
			const Real* xhi_dat, const Real* slx_dat, 
			Real* slxscr, Real* stxlo, Real* stxhi,
			const Real* uedge_dat, ARLIM_P(ue_lo), ARLIM_P(ue_hi),
			const Real* stx,       ARLIM_P(stx_lo),ARLIM_P(stx_hi),
			const Real* ylo_dat, 
			const Real* yhi_dat, const Real* sly_dat, 
			Real* slyscr, Real* stylo, Real* styhi,
			const Real* vedge_dat, ARLIM_P(ve_lo), ARLIM_P(ve_hi),
			const Real* sty,       ARLIM_P(sty_lo),ARLIM_P(sty_hi),
#if (BL_SPACEDIM == 3)
			const Real* zlo_dat, 
			const Real* zhi_dat, const Real* slz_dat, 
			Real* slzscr, Real* stzlo, Real* stzhi,
			const Real* wedge_dat, ARLIM_P(we_lo), ARLIM_P(we_hi),
			const Real* stz,       ARLIM_P(stz_lo),ARLIM_P(stz_hi),
#endif
			ARLIM_P(ww_lo),ARLIM_P(ww_hi),
			const int* bc, const int* lo, const int* hi, 
			Real* dt, const Real* dx, int* fort_ind, 
			int* use_forces_in_trans, int* iconserv);
  
   void FORT_ESTATE_LIN(const Real* s_dat,  ARLIM_P(s_lo),  ARLIM_P(s_hi),
			const Real* sn_dat, ARLIM_P(sn_lo), ARLIM_P(sn_hi),
			const Real* tfr_dat, 
			const Real* rphi_dat,  ARLIM_P(rphi_lo),  ARLIM_P(rphi_hi),
			const Real* slx_dat, Real* slxscr,
			const Real* uedge_dat, ARLIM_P(ue_lo), ARLIM_P(ue_hi),
			const Real* stx,       ARLIM_P(stx_lo),ARLIM_P(stx_hi),
			const Real* sly_dat, 
			Real* slyscr,
			const Real* vedge_dat, ARLIM_P(ve_lo), ARLIM_P(ve_hi),
			const Real* sty,       ARLIM_P(sty_lo),ARLIM_P(sty_hi),
#if (BL_SPACEDIM == 3)
			const Real* slz_dat, 
			Real* slzscr,
			const Real* wedge_dat, ARLIM_P(we_lo), ARLIM_P(we_hi),
			const Real* stz,       ARLIM_P(stz_lo),ARLIM_P(stz_hi), 
#endif
			ARLIM_P(ww_lo),ARLIM_P(ww_hi),
			const int* bc, const int* lo, const int* hi, 
			Real* dt, const Real* dx,  int* nscal);

   void FORT_ESTATE_RMN(const Real* s_dat,  ARLIM_P(s_lo),  ARLIM_P(s_hi),
			const Real* sn_dat, ARLIM_P(sn_lo), ARLIM_P(sn_hi),
			const Real* tfr_dat, const Real* divu_dat, 
			const Real* lbd_dat, const Real* dlbd_dat,
			const Real* rphi_dat,  ARLIM_P(rphi_lo),  ARLIM_P(rphi_hi),
			const Real* kappa_dat, ARLIM_P(kappa_lo), ARLIM_P(kappa_hi),
			const Real* slx_dat, Real* slxscr,
			const Real* uedge_dat, ARLIM_P(ue_lo), ARLIM_P(ue_hi),
			const Real* stx,       ARLIM_P(stx_lo),ARLIM_P(stx_hi),
			const Real* kappax,    ARLIM_P(kpx_lo),ARLIM_P(kpx_hi),
			const Real* eigvx, const Real* eiglx, const Real* eigrx,
			const Real* sly_dat, 
			Real* slyscr,
			const Real* vedge_dat, ARLIM_P(ve_lo), ARLIM_P(ve_hi),
			const Real* sty,       ARLIM_P(sty_lo),ARLIM_P(sty_hi),
			const Real* kappay,    ARLIM_P(kpy_lo),ARLIM_P(kpy_hi),
			const Real* eigvy, const Real* eigly, const Real* eigry,
#if (BL_SPACEDIM == 3)
			const Real* slz_dat, 
			Real* slzscr,
			const Real* wedge_dat, ARLIM_P(we_lo), ARLIM_P(we_hi),
			const Real* stz,       ARLIM_P(stz_lo),ARLIM_P(stz_hi),
			const Real* kappaz,    ARLIM_P(kpz_lo),ARLIM_P(kpz_hi),
			const Real* eigvz, const Real* eiglz, const Real* eigrz, 
#endif
			ARLIM_P(ww_lo),ARLIM_P(ww_hi),
			const Real* krcoef,    ARLIM_P(kr_lo), ARLIM_P(kr_hi),
			const int* nkrcoef,
			const int* bc, const int* lo, const int* hi, 
			Real* dt, const Real* dx, 
			int* use_forces_in_trans, int* iconserv, int* nscal);

   void FORT_ESTATE_CPL(const Real* s_dat,  ARLIM_P(s_lo),  ARLIM_P(s_hi),
			const Real* sn_dat, ARLIM_P(sn_lo), ARLIM_P(sn_hi),
			const Real* pc_dat, ARLIM_P(pc_lo), ARLIM_P(pc_hi),
			const Real* tfr_dat, const Real* divu_dat, 
			const Real* lbd_dat, const Real* dlbd_dat,
			const Real* rphi_dat,  ARLIM_P(rphi_lo),  ARLIM_P(rphi_hi),
			const Real* kappa_dat, ARLIM_P(kappa_lo), ARLIM_P(kappa_hi),
			const Real* slx_dat, Real* slxscr,
			const Real* uedge_dat, ARLIM_P(ue_lo),  ARLIM_P(ue_hi),
			const Real* stx,       ARLIM_P(stx_lo), ARLIM_P(stx_hi),
			const Real* kappax,    ARLIM_P(kpx_lo), ARLIM_P(kpx_hi),
			const Real* lambdax,   ARLIM_P(ldx_lo), ARLIM_P(ldx_hi),
			const Real* eigvx, const Real* eiglx, const Real* eigrx,
			const Real* sly_dat, 
			Real* slyscr,
			const Real* vedge_dat, ARLIM_P(ve_lo),  ARLIM_P(ve_hi),
			const Real* sty,       ARLIM_P(sty_lo), ARLIM_P(sty_hi),
			const Real* kappay,    ARLIM_P(kpy_lo), ARLIM_P(kpy_hi),		
			const Real* lambday,   ARLIM_P(ldy_lo), ARLIM_P(ldy_hi),
			const Real* eigvy, const Real* eigly, const Real* eigry,
#if (BL_SPACEDIM == 3)
			const Real* slz_dat, 
			Real* slzscr,
			const Real* wedge_dat, ARLIM_P(we_lo),  ARLIM_P(we_hi),
			const Real* stz,       ARLIM_P(stz_lo), ARLIM_P(stz_hi),
			const Real* kappaz,    ARLIM_P(kpz_lo), ARLIM_P(kpz_hi),		
			const Real* lambdaz,   ARLIM_P(ldz_lo), ARLIM_P(ldz_hi),
			const Real* eigvz, const Real* eiglz, const Real* eigrz, 
#endif
			ARLIM_P(ww_lo),ARLIM_P(ww_hi),
			const Real* krcoef,    ARLIM_P(kr_lo), ARLIM_P(kr_hi),
			const int* nkrcoef,
			const int* bc, const int* lo, const int* hi, 
			Real* dt, const Real* dx, 
			int* use_forces_in_trans, int* iconserv, int* nscal);

   void FORT_ESTATE_PMR(const Real* s_dat, const Real* u_dat, ARLIM_P(s_lo),ARLIM_P(s_hi), 
			const Real* sn_dat, const Real* un_dat, ARLIM_P(sn_lo),ARLIM_P(sn_hi),
			const Real* tfr_dat, const Real* divu_dat,  
			const Real* rphi_dat,  ARLIM_P(rphi_lo),  ARLIM_P(rphi_hi),
			const Real* kappa_dat, ARLIM_P(kappa_lo), ARLIM_P(kappa_hi),
			const Real* slx_dat, 
			Real* slxscr,
			const Real* uedge_dat, ARLIM_P(ue_lo), ARLIM_P(ue_hi),
			const Real* stx,       ARLIM_P(stx_lo),ARLIM_P(stx_hi),
			const Real* kappax,    ARLIM_P(kpx_lo),ARLIM_P(kpx_hi),
			const Real* strcx, const Real* eigvx, 
			const Real* eiglx, const Real* eigrx, 
			const Real* sly_dat, 
			Real* slyscr,
			const Real* vedge_dat, ARLIM_P(ve_lo), ARLIM_P(ve_hi),
			const Real* sty,       ARLIM_P(sty_lo),ARLIM_P(sty_hi),
			const Real* kappay,    ARLIM_P(kpy_lo),ARLIM_P(kpy_hi),
			const Real* strcy, const Real* eigvy, 
			const Real* eigly, const Real* eigry, 
#if (BL_SPACEDIM == 3)
			const Real* slz_dat, 
			Real* slzscr,
			const Real* wedge_dat, ARLIM_P(we_lo), ARLIM_P(we_hi),
			const Real* stz,       ARLIM_P(stz_lo),ARLIM_P(stz_hi),
			const Real* kappaz,    ARLIM_P(kpz_lo),ARLIM_P(kpz_hi),
			const Real* strcz, const Real* eigvz, 
			const Real* eiglz, const Real* eigrz, 
#endif
			ARLIM_P(ww_lo),ARLIM_P(ww_hi),
			const int* bc, const int* lo, const int* hi, 
			Real* dt, const Real* dx, 
			int* use_forces_in_trans, int* iconserv,
			Real* gravity, Real* eigmax);

   void FORT_ESTATE_TRACER(const Real* s_dat,   ARLIM_P(s_lo),ARLIM_P(s_hi),
			   const Real* sn_dat,  ARLIM_P(sn_lo),ARLIM_P(sn_hi),
			   const Real* st_dat,  ARLIM_P(st_lo),ARLIM_P(st_hi),
			   const Real* stn_dat, ARLIM_P(stn_lo),ARLIM_P(stn_hi),
			   const Real* slx_dat, Real* slxscr, 
			   const Real* xlo_dat, const Real* xhi_dat,
			   const Real* uedge_dat, ARLIM_P(ue_lo), ARLIM_P(ue_hi),
			   Real* stx,             ARLIM_P(stx_lo),ARLIM_P(stx_hi),
			   const Real* sly_dat, Real* slyscr,
			   const Real* ylo_dat, const Real* yhi_dat,
			   const Real* vedge_dat, ARLIM_P(ve_lo), ARLIM_P(ve_hi),
			   Real* sty,             ARLIM_P(sty_lo),ARLIM_P(sty_hi),
#if (BL_SPACEDIM == 3)
			   const Real* slz_dat, Real* slzscr,
			   const Real* zlo_dat, const Real* zhi_dat,
			   const Real* wedge_dat, ARLIM_P(we_lo), ARLIM_P(we_hi),
			   Real* stz,             ARLIM_P(stz_lo),ARLIM_P(stz_hi),
#endif
			   const Real* rphi_dat,  ARLIM_P(rphi_lo),  ARLIM_P(rphi_hi),
			   ARLIM_P(ww_lo), ARLIM_P(ww_hi),
			   const int* bc, const int* lo, const int* hi, 
			   Real* dt, const Real* dx, int* nscal);

  void FORT_SYNC_RMN(const Real* s_dat,  ARLIM_P(s_lo),  ARLIM_P(s_hi),
		     const Real* sn_dat, ARLIM_P(sn_lo), ARLIM_P(sn_hi),
		     const Real* tfr_dat, const Real* divu_dat, 
		     const Real* lbd_dat, const Real* dlbd_dat,
		     const Real* rphi_dat,  ARLIM_P(rphi_lo),  ARLIM_P(rphi_hi),
		     const Real* kappa_dat, ARLIM_P(kappa_lo), ARLIM_P(kappa_hi),
		     const Real* slx_dat, Real* slxscr,
		     const Real* uedge_dat, ARLIM_P(ue_lo), ARLIM_P(ue_hi),
		     const Real* stx,       ARLIM_P(stx_lo),ARLIM_P(stx_hi),
		     const Real* kappax,    ARLIM_P(kpx_lo),ARLIM_P(kpx_hi),
		     const Real* eigvx, const Real* eiglx, const Real* eigrx,
		     const Real* sly_dat, Real* slyscr,
		     const Real* vedge_dat, ARLIM_P(ve_lo), ARLIM_P(ve_hi),
		     const Real* sty,       ARLIM_P(sty_lo),ARLIM_P(sty_hi),
		     const Real* kappay,    ARLIM_P(kpy_lo),ARLIM_P(kpy_hi),
		     const Real* eigvy, const Real* eigly, const Real* eigry,
#if (BL_SPACEDIM == 3)
		     const Real* slz_dat, 
		     Real* slzscr,
		     const Real* wedge_dat, ARLIM_P(we_lo), ARLIM_P(we_hi),
		     const Real* stz,       ARLIM_P(stz_lo),ARLIM_P(stz_hi),
		     const Real* kappaz,    ARLIM_P(kpz_lo),ARLIM_P(kpz_hi),
		     const Real* eigvz, const Real* eiglz, const Real* eigrz, 
#endif
		     ARLIM_P(ww_lo),ARLIM_P(ww_hi),
		     const Real* krcoef,    ARLIM_P(kr_lo), ARLIM_P(kr_hi),
		     const int* nkrcoef,
		     const int* bc, const int* lo, const int* hi, 
		     Real* dt, const Real* dx, 
		     int* use_forces_in_trans, int* iconserv, int* nscal);

   void FORT_ADV_FORCING(const Real* aofs_dat, ARLIM_P(a_lo), ARLIM_P(a_hi),
			 const Real* xflux_dat, ARLIM_P(xflux_lo),ARLIM_P(xflux_hi),
                         const Real* uedge_dat, ARLIM_P(uedge_lo),ARLIM_P(uedge_hi),
			 const Real* areax_dat, ARLIM_P(ax_lo), ARLIM_P(ax_hi),
			 const Real* yflux_dat, ARLIM_P(yflux_lo),ARLIM_P(yflux_hi),
                         const Real* vedge_dat, ARLIM_P(vedge_lo),ARLIM_P(vedge_hi),
			 const Real* areay_dat, ARLIM_P(ay_lo), ARLIM_P(ay_hi),
#if (BL_SPACEDIM == 3)                                                    
			 const Real* zflux_dat, ARLIM_P(zflux_lo),ARLIM_P(zflux_hi),
                         const Real* wedge_dat, ARLIM_P(wedge_lo),ARLIM_P(wedge_hi),
			 const Real* areaz_dat, ARLIM_P(az_lo), ARLIM_P(az_hi),
#endif
			 const Real* vol_dat, ARLIM_P(v_lo),ARLIM_P(v_hi),
			 const int* lo, const int* hi, int* iconserv);

   void FORT_ADV_RMN_FORCING(const Real* aofs_dat, ARLIM_P(a_lo), ARLIM_P(a_hi),
			     const Real* xflux_dat, ARLIM_P(xflux_lo),ARLIM_P(xflux_hi),
			     const Real* uedge_dat, ARLIM_P(uedge_lo),ARLIM_P(uedge_hi),
			     const Real* areax_dat, ARLIM_P(ax_lo), ARLIM_P(ax_hi),
			     const Real* yflux_dat, ARLIM_P(yflux_lo),ARLIM_P(yflux_hi),
			     const Real* vedge_dat, ARLIM_P(vedge_lo),ARLIM_P(vedge_hi),
			     const Real* areay_dat, ARLIM_P(ay_lo), ARLIM_P(ay_hi),
#if (BL_SPACEDIM == 3)                                                    
			     const Real* zflux_dat, ARLIM_P(zflux_lo),ARLIM_P(zflux_hi),
			     const Real* wedge_dat, ARLIM_P(wedge_lo),ARLIM_P(wedge_hi),
			     const Real* areaz_dat, ARLIM_P(az_lo), ARLIM_P(az_hi),
#endif
			     const Real* vol_dat, ARLIM_P(v_lo),ARLIM_P(v_hi),
			     const int* lo, const int* hi,  int * nscal);

   
   void FORT_SYNC_ADV_FORCING(const Real* sync_dat,ARLIM_P(s_lo),ARLIM_P(s_hi),
			      const Real* xflux_dat, ARLIM_P(xflux_lo),ARLIM_P(xflux_hi),
                              const Real* ucorr_dat, ARLIM_P(ucorr_lo),ARLIM_P(ucorr_hi),
			      const Real* areax_dat, ARLIM_P(ax_lo), ARLIM_P(ax_hi),
			      const Real* yflux_dat, ARLIM_P(yflux_lo),ARLIM_P(yflux_hi),
                              const Real* vcorr_dat, ARLIM_P(vcorr_lo),ARLIM_P(vcorr_hi),
			      const Real* areay_dat, ARLIM_P(ay_lo), ARLIM_P(ay_hi),
#if (BL_SPACEDIM == 3)                                             
			      const Real* zflux_dat, ARLIM_P(zflux_lo),ARLIM_P(zflux_hi),
                              const Real* wcorr_dat, ARLIM_P(wcorr_lo),ARLIM_P(wcorr_hi),
			      const Real* areaz_dat, ARLIM_P(az_lo), ARLIM_P(az_hi),
#endif 		   
			      const Real* vol_dat, ARLIM_P(v_lo),ARLIM_P(v_hi),
			      const int* lo, const int* hi, int* iconserv);

   void FORT_SYNC_RMN_FORCING(const Real* sync_dat,ARLIM_P(s_lo),ARLIM_P(s_hi),
			      const Real* xflux_dat, ARLIM_P(xflux_lo),ARLIM_P(xflux_hi),
                              const Real* ucorr_dat, ARLIM_P(ucorr_lo),ARLIM_P(ucorr_hi),
			      const Real* areax_dat, ARLIM_P(ax_lo), ARLIM_P(ax_hi),
			      const Real* yflux_dat, ARLIM_P(yflux_lo),ARLIM_P(yflux_hi),
                              const Real* vcorr_dat, ARLIM_P(vcorr_lo),ARLIM_P(vcorr_hi),
			      const Real* areay_dat, ARLIM_P(ay_lo), ARLIM_P(ay_hi),
#if (BL_SPACEDIM == 3)                                             
			      const Real* zflux_dat, ARLIM_P(zflux_lo),ARLIM_P(zflux_hi),
                              const Real* wcorr_dat, ARLIM_P(wcorr_lo),ARLIM_P(wcorr_hi),
			      const Real* areaz_dat, ARLIM_P(az_lo), ARLIM_P(az_hi),
#endif 
			      const Real* krcoef,  ARLIM_P(kr_lo), ARLIM_P(kr_hi),
			      const int* nkrcoef,
			      const Real* vol_dat, ARLIM_P(v_lo),ARLIM_P(v_hi),
			      const int* lo, const int* hi, const int* nscal);
   
   void FORT_CONVSCALMINMAX (Real* Sold_dat, ARLIM_P(slo), ARLIM_P(shi),
                             Real* Snew_dat, ARLIM_P(snlo), ARLIM_P(snhi),
#if (BL_SPACEDIM == 3)
			     const Real* smin_dat, const Real* smax_dat,
			     ARLIM_P(lo),ARLIM_P(hi),
#endif
			     const int* lo, const int* hi, const int* bc);

   void FORT_CONSSCALMINMAX (Real* Sold_dat, Real* Snew_dat, 
                             ARLIM_P(slo), ARLIM_P(shi),
#if (BL_SPACEDIM == 3)
			     const Real* smin_dat, const Real* smax_dat,
			     ARLIM_P(lo),ARLIM_P(hi),
#endif
			     const int* lo, const int* hi, const int* bc);
   
   void  FORT_UPDATE_TF(const Real* SOdat, ARLIM_P(slo), ARLIM_P(shi), 
			const Real* SNdat, ARLIM_P(slo), ARLIM_P(shi),
			const Real* Tfdat, ARLIM_P(tlo), ARLIM_P(thi),
			const int* lo, const int* hi, Real* dt, int* num_comp);
   
   void FORT_UPDATE_AOFS_TF(const Real* SOdat, ARLIM_P(slo), ARLIM_P(shi), 
			    Real*       SNdat, ARLIM_P(slo), ARLIM_P(shi),
			    const Real* AOdat, ARLIM_P(alo), ARLIM_P(ahi),
			    const Real*  pdat, ARLIM_P(plo), ARLIM_P(phi),
			    const Real* tfdat, ARLIM_P(tlo), ARLIM_P(thi),
			    const int* lo, const int* hi, 
			    Real* dt, int* num_comp);

   void FORT_UPDATE_AOFS_TRACER(const Real* COdat, ARLIM_P(colo), ARLIM_P(cohi), 
				Real*       CNdat, ARLIM_P(cnlo), ARLIM_P(cnhi),
				const Real* SOdat, ARLIM_P(solo), ARLIM_P(sohi), 
				const Real* SNdat, ARLIM_P(snlo), ARLIM_P(snhi),
				const Real* AOdat, ARLIM_P(alo), ARLIM_P(ahi),
				const Real* tfdat, ARLIM_P(tlo), ARLIM_P(thi),
				const Real*  pdat, ARLIM_P(plo), ARLIM_P(phi),
				const int* lo, const int* hi, 
				const int* idx_total, const int* idx_total_n,
				const Real* dt, const int* ncomps, const int* ntracers);
   
   void FORT_SUM_TF_DIVU_VISC(const Real* Sdat,  ARLIM_P(slo), ARLIM_P(shi),
			      Real*       Tfdat, ARLIM_P(tlo), ARLIM_P(thi),
			      const Real* DUdat, ARLIM_P(dlo), ARLIM_P(dhi),
			      const Real* VIdat, ARLIM_P(vlo), ARLIM_P(vhi),
			      const int* tlo, const int* thi, 
			      const int* num_comp, const int* iconserv);
   
   void FORT_SUM_TF_DIVU(const Real* Sdat,  ARLIM_P(slo), ARLIM_P(shi),
			 Real*       Tfdat, ARLIM_P(tlo), ARLIM_P(thi),
			 const Real* DUdat, ARLIM_P(dlo), ARLIM_P(dhi),
			 const int* tlo, const int* thi, 
			 const int* num_comp, const int* iconserv);
   
   void FORT_SET_PARAMS(const int& slope_order,
			const int& use_unlimited_slopes);

   void FORT_GET_EIG(const Real* Sdat, ARLIM_P(slo), ARLIM_P(shi),
		     const Real* Fdat, ARLIM_P(wlo), ARLIM_P(whi),
		     int* ncomps, int* idx);
}
#endif

#endif /*_GODUNOV_F_H_*/
