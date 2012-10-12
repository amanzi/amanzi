#ifndef _POROUS_F_H_
#define _POROUS_F_H_

#ifdef BL_LANG_FORT
#    define FORT_INITKAPPA   initkappa
#    define FORT_INITKAPPA2  initkappa2
#    define FORT_INITPHI     initphi
#    define FORT_INITPHI2    initphi2
#    define FORT_CRSENMAT    crsenmat
#    define FORT_INITKR      initkr
#    define FORT_MK_MACCOEF  mkmaccoef
#    define FORT_TURN        turn
#    define FORT_PHIRAND     phirand
#    define FORT_PERIODIC    periodicrand
#else
#  if  defined(BL_FORT_USE_UPPERCASE)
#    define FORT_INITKAPPA   INITKAPPA
#    define FORT_INITKAPPA2  INITKAPPA2
#    define FORT_INITPHI     INITPHI
#    define FORT_INITPHI2    INITPHI2
#    define FORT_CRSENMAT    CRSENKMAT
#    define FORT_INITKR      INITKR
#    define FORT_MK_MACCOEF  MKMACCOEF
#    define FORT_TURN        TURN
#    define FORT_PHIRAND     PHIRAND
#    define FORT_PERIODIC    PERIODICRAND
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_INITKAPPA   initkappa
#    define FORT_INITPHI     initphi
#    define FORT_INITPHI2    initphi2
#    define FORT_CRSENMAT    crsenmat
#    define FORT_INITKR      initkr
#    define FORT_MK_MACCOEF  mkmaccoef
#    define FORT_TURN        turn
#    define FORT_PHIRAND     phirand
#    define FORT_PERIODIC    periodicrand
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_INITKAPPA   initkappa_
#    define FORT_INITKAPPA2  initkappa2_
#    define FORT_INITPHI     initphi_
#    define FORT_INITPHI2    initphi2_
#    define FORT_CRSENMAT    crsenmat_
#    define FORT_INITKR      initkr_
#    define FORT_MK_MACCOEF  mkmaccoef_
#    define FORT_TURN        turn_
#    define FORT_PHIRAND     phirand_
#    define FORT_PERIODIC    periodicrand_
#  endif

#include "ArrayLim.H"

extern "C" {

  void FORT_INITKAPPA (const Real* kappa,ARLIM_P(k_lo), ARLIM_P(k_hi),
		       const Real* kappax,ARLIM_P(kx_lo), ARLIM_P(kx_hi),
		       const Real* kappay,ARLIM_P(ky_lo), ARLIM_P(ky_hi),
#if (BL_SPACEDIM == 3)
		       const Real* kappaz,ARLIM_P(kz_lo), ARLIM_P(kz_hi),
#endif
		       const int* lo, const int* hi,
		       const Real* dx, const Real* probhi,
		       const int* currlev, const int* maxlev, const int* mediumtype, 
		       const Real* kappaval_x, const Real* kappaval_y,
#if (BL_SPACEDIM == 3)  
		       const Real* kappaval_z,
#endif
		       const int* nlayer, int* fratio);

  void FORT_INITKAPPA2 (const Real* kappadat,ARLIM_P(kd_lo), ARLIM_P(kd_hi),
			const Real* kappa,ARLIM_P(k_lo), ARLIM_P(k_hi),
			const Real* kappax,ARLIM_P(kx_lo), ARLIM_P(kx_hi),
			const Real* kappay,ARLIM_P(ky_lo), ARLIM_P(ky_hi),
#if (BL_SPACEDIM == 3)
			const Real* kappaz,ARLIM_P(kz_lo), ARLIM_P(kz_hi),
#endif
			const int* lo, const int* hi,
			const int* currlev, const int* maxlev, int* fratio);


  void FORT_INITPHI (const Real* rock_phi,
		     ARLIM_P(p_lo), ARLIM_P(p_hi),
		     const int* domhi, const Real* dx, 
		     const int* mediumtype, const Real* kappaval, 
		     const int* nlayer);

  void FORT_INITPHI2 (const Real* phidata,
		      ARLIM_P(pd_lo), ARLIM_P(pd_hi),
		      const Real* rock_phi,
		      ARLIM_P(p_lo), ARLIM_P(p_hi),
		      const int* lo, const int* hi,
		      const int* fratio);

  void FORT_CRSENMAT (const Real* kfine, ARLIM_P(kflo), ARLIM_P(kfhi),
                      const Real* kcrse, ARLIM_P(kclo), ARLIM_P(kchi),
                      const int* lo, const int* hi, 
                      const int* fratio, const int* harmdir);

  void FORT_INITKR (const Real* krdata, ARLIM_P(krd_lo), ARLIM_P(krd_hi),
		    const Real* kr, ARLIM_P(kr_lo), ARLIM_P(kr_hi), 
		    const int* nkr, const int* currlev, 
		    const int* maxlev, int* fractio);

    
  void FORT_MK_MACCOEF ( const Real* lx, ARLIM_P(lxlo), ARLIM_P(lxhi),
			 const Real* ly, ARLIM_P(lylo), ARLIM_P(lyhi),
#if (BL_SPACEDIM == 3)
			 const Real* lz, ARLIM_P(lzlo), ARLIM_P(lzhi),
#endif
			 const Real* kappax,ARLIM_P(kx_lo), ARLIM_P(kx_hi),
			 const Real* kappay,ARLIM_P(ky_lo), ARLIM_P(ky_hi),
#if (BL_SPACEDIM == 3)
			 const Real* kappaz,ARLIM_P(kz_lo), ARLIM_P(kz_hi),
#endif
			 const Real* lbd,ARLIM_P(l_lo), ARLIM_P(l_hi),
			 const int* lo, const int* hi,
			 const int* domlo, const int* domhi, const int* bc);

  void FORT_TURN (const Real* u, const int* nx, const int* ny, 
		  const int* nz, const int* iuc);

  void FORT_RESCALE (const Real* kp_dat,   ARLIM_P(kp_lo), ARLIM_P(kp_hi), 
		     const Real* kappaval, const int* nkpval, 
		     const int*  domlo,    const int* domhi);

  void FORT_PHIRAND (const Real* kp_dat,   ARLIM_P(kp_lo), ARLIM_P(kp_hi), 
		     const Real* kappaval, const int* nkpval, const Real* dkappa,
		     const int*  domlo,    const int* domhi, const int* iuc);
  
}
#endif
#endif /*_POROUS_F_H_*/
