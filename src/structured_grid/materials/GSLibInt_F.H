/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#ifndef _GSLibInt_F_H_
#define _GSLibInt_F_H_

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
#    define FORT_PHIRAND         phi_rand
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
#    define FORT_PHIRAND         PHI_RAND
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
#    define FORT_PHIRAND         phi_rand
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
#    define FORT_PHIRAND         phi_rand_
#  endif

#include "ArrayLim.H"

  extern "C" {


    void FORT_INIT_GSLIB(const char* paramfl,const int* paramsz,
			 int* nx, Real* xsiz, Real* xmn,
                         int* rand_seed);

    void FORT_INIT_GSLIB2(const char* paramfl,const int* paramsz,
			  int* c_sz   , int* c_idx, int* c_idx_siz,
			  int* real_sz, int* real_idx, int* real_idx_siz,
			  int* int_sz , int* int_idx, int* int_idx_siz, int* cond_option,
                          int* rand_seed);

    void FORT_INTERNAL_DATA(const Real* scratch_c, int* c_sz, int* c_idx, const int* c_idx_siz);

    void FORT_SEC_DATA(const Real* dat, ARLIM_P(dlo), ARLIM_P(dhi),
		       const Real* scratch_r, int* real_sz, int* real_idx,
		       const int* r_idx_siz, int* ratios, Real* condg, Real* sig, Real* dx);

    void FORT_AVE_DATA(const Real* dat, ARLIM_P(dlo), ARLIM_P(dhi),
		       const Real* sdat, ARLIM_P(slo), ARLIM_P(shi),
		       int* ratios, Real* condg, Real* sig, Real* dx);

    void FORT_SGSIM(const Real* dat, ARLIM_P(dlo), ARLIM_P(dhi),
		    Real* dx);

    void FORT_SGSIM_SETUP(const Real* dat, ARLIM_P(dlo), ARLIM_P(dhi),
			  const int*  order, int* nvalid,
			  const Real* scratch_c, int* c_sz, int* c_idx, const int* c_idx_siz,
			  const Real* scratch_r, int* real_sz, int* real_idx, const int* r_idx_siz,
			  const int*  scratch_i, int* int_sz, int* int_idx, const int* i_idx_siz,
			  const int*  lo, const int* hi,
			  const Real* dx, const Real* plo,
			  const int* dlo, const int* dhi, int* iseed);
    void FORT_SGSIM_ITER(const Real* dat, ARLIM_P(dlo), ARLIM_P(dhi),
			 const Real* scratch_c, int* c_sz, int* c_idx, const int* c_idx_siz,
			 const Real* scratch_r, int* real_sz, int* real_idx, const int* r_idx_siz,
			 const int*  scratch_i, int* int_sz, int* int_idx, const int* i_idx_siz,
			 int* iter);
    void FORT_SGSIM_POST(const Real* dat, ARLIM_P(dlo), ARLIM_P(dhi),
			 const Real* scratch_c, int* c_sz, int* c_idx, const int* c_idx_siz,
			 const Real* scratch_r, int* real_sz, int* real_idx, const int* r_idx_siz,
			 const int*  scratch_i, int* int_sz, int* int_idx,  const int* i_idx_siz);

    void FORT_LGNORM(const Real* dat, ARLIM_P(dlo), ARLIM_P(dhi),
		     const Real* condg, const int* lo, const int* hi);

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

  void FORT_PHIRAND (const Real* kp_dat,   ARLIM_P(kp_lo), ARLIM_P(kp_hi),
		     const Real* kappaval, const int* nkpval, const Real* dkappa,
		     const int*  domlo,    const int* domhi, const int* iuc);

  }
#endif
#endif /*_GSLibInt_F_H_*/
