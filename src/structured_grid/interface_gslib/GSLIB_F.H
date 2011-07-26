#ifndef _GSLIB_F_H_
#define _GSLIB_F_H_

#ifdef BL_LANG_FORT
#    define FORT_INIT_GSLIB      init_gslib
#    define FORT_INIT_GSLIB2     init_gslib2
#    define FORT_INTERNAL_DATA   internal_data
#    define FORT_SEC_DATA        sec_data
#    define FORT_AVE_DATA        ave_data
#    define FORT_SGSIM           sgsim_gslib
#    define FORT_SGSIM_SETUP     sgsim_gslib_setup
#    define FORT_SGSIM_ITER      sgsim_gslib_iter
#    define FORT_SGSIM_POST      sgsim_gslib_post
#    define FORT_LGNORM          lgnorm
#    define FORT_SGSIM_DEALLOC   sgsim_dealloc
#    define FORT_INTERP_LAYER    interp_layer
#    define FORT_ASSIGN_LAYER    assign_layer
#else
#  if defined(BL_FORT_USE_UPPERCASE)
#    define FORT_INIT_GSLIB      INIT_GSLIB
#    define FORT_INIT_GSLIB2     INIT_GSLIB2
#    define FORT_INTERNAL_DATA   INTERNAL_DATA
#    define FORT_SEC_DATA        SEC_DATA
#    define FORT_AVE_DATA        AVE_DATA
#    define FORT_SGSIM           SGSIM_GSLIB
#    define FORT_SGSIM_SETUP     SGSIM_GSLIBs_SETUP
#    define FORT_SGSIM_ITER      SGSIM_GSLIBs_ITER
#    define FORT_SGSIM_POST      SGSIM_GSLIBs_POST
#    define FORT_LGNORM          LGNORM
#    define FORT_SGSIM_DEALLOC   SGSIM_DEALLOC
#    define FORT_INTERP_LAYER    INTERP_LAYER
#    define FORT_ASSIGN_LAYER    ASSIGN_LAYER
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_INIT_GSLIB      init_gslib
#    define FORT_INIT_GSLIB2     init_gslib2
#    define FORT_INTERNAL_DATA   internal_data
#    define FORT_SEC_DATA        sec_data
#    define FORT_AVE_DATA        ave_data
#    define FORT_SGSIM           sgsim_gslib
#    define FORT_SGSIM_SETUP     sgsim_gslib_setup
#    define FORT_SGSIM_ITER      sgsim_gslib_iter
#    define FORT_SGSIM_POST      sgsim_gslib_post
#    define FORT_LGNORM          lgnorm
#    define FORT_SGSIM_DEALLOC   sgsim_dealloc
#    define FORT_INTERP_LAYER    interp_layer
#    define FORT_ASSIGN_LAYER    assign_layer
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_INIT_GSLIB      init_gslib_
#    define FORT_INIT_GSLIB2     init_gslib2_
#    define FORT_INTERNAL_DATA   internal_data_
#    define FORT_SEC_DATA        sec_data_
#    define FORT_AVE_DATA        ave_data_
#    define FORT_SGSIM           sgsim_gslib_
#    define FORT_SGSIM_SETUP     sgsim_gslib_setup_
#    define FORT_SGSIM_ITER      sgsim_gslib_iter_
#    define FORT_SGSIM_POST      sgsim_gslib_post_
#    define FORT_LGNORM          lgnorm_
#    define FORT_SGSIM_DEALLOC   sgsim_dealloc_
#    define FORT_INTERP_LAYER    interp_layer_
#    define FORT_ASSIGN_LAYER    assign_layer_
#  endif

#include "ArrayLim.H"

  extern "C" {


    void FORT_INIT_GSLIB(const char* paramfl,const int* paramsz,
			 int* nx, Real* xsiz, Real* xmn);

    void FORT_INIT_GSLIB2(const char* paramfl,const int* paramsz,
			  int* c_sz   , int* c_idx,
			  int* real_sz, int* real_idx,
			  int* int_sz , int* int_idx, int* cond_option);

    void FORT_INTERNAL_DATA(const Real* dat, ARLIM_P(dlo), ARLIM_P(dhi),
			    const Real* scratch_c, int* c_sz, int* c_idx,
			    Real* condg, Real* sig, Real* dx,
			    const int* domlo, const int* domhi);

    void FORT_SEC_DATA(const Real* dat, ARLIM_P(dlo), ARLIM_P(dhi),
		       const Real* scratch_r, int* real_sz, int* real_idx,
		       int* ratios, Real* condg, Real* sig, Real* dx);

    void FORT_AVE_DATA(const Real* dat, ARLIM_P(dlo), ARLIM_P(dhi),
		       const Real* sdat, ARLIM_P(slo), ARLIM_P(shi),
		       int* ratios, Real* condg, Real* sig, Real* dx);

    void FORT_SGSIM(const Real* dat, ARLIM_P(dlo), ARLIM_P(dhi),
		    Real* dx);

    void FORT_SGSIM_SETUP(const Real* dat, ARLIM_P(dlo), ARLIM_P(dhi),
			  const int*  order, int* nvalid,
			  const Real* scratch_c, int* c_sz, int* c_idx,
			  const Real* scratch_r, int* real_sz, int* real_idx,
			  const int*  scratch_i, int* int_sz, int* int_idx,
			  const int*  lo, const int* hi,
			  Real* dx, int* iseed);
    void FORT_SGSIM_ITER(const Real* dat, ARLIM_P(dlo), ARLIM_P(dhi),
			 const Real* scratch_c, int* c_sz, int* c_idx,
			 const Real* scratch_r, int* real_sz, int* real_idx,
			 const int*  scratch_i, int* int_sz, int* int_idx,
			 int* iter);
    void FORT_SGSIM_POST(const Real* dat, ARLIM_P(dlo), ARLIM_P(dhi),
			 const Real* scratch_c, int* c_sz, int* c_idx,
			 const Real* scratch_r, int* real_sz, int* real_idx,
			 const int*  scratch_i, int* int_sz, int* int_idx);

    void FORT_LGNORM(const Real* dat, ARLIM_P(dlo), ARLIM_P(dhi),
		     Real* condg, Real* sig);

    void FORT_SGSIM_DEALLOC();

#if (BL_SPACEDIM == 3)

    void FORT_INTERP_LAYER (const Real* zl_dat,   ARLIM_P(zl_lo), ARLIM_P(zl_hi), 
			    const Real* xs, const Real* ys, const Real* zs, 
			    const int* ns, const int* nl, 
			    const Real* dx, const Real* problo, const Real* probhi);

#else

    void FORT_INTERP_LAYER (const Real* zl_dat,   ARLIM_P(zl_lo), ARLIM_P(zl_hi), 
			    const Real* xs, const Real* zs, const int*  ns, const int* nl, 
			    const Real* dx, const Real* problo, const Real* probhi);

#endif

    void FORT_ASSIGN_LAYER (const Real* kp_dat,   ARLIM_P(kp_lo), ARLIM_P(kp_hi), 
			    const Real* zl_dat,   ARLIM_P(zl_lo), ARLIM_P(zl_hi),
			    const Real* kappaval, const int* nkpval, const Real* dkappa,
			    const Real* dx, const Real* problo, const Real* probhi);


  }
#endif
#endif /*_GSLIB_F_H_*/
