/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#ifndef _PROB_PM_F_H_
#define _PROB_PM_F_H_

#ifdef BL_LANG_FORT
#    define FORT_INITPARAMS  initparams
#    define FORT_TCRPARAMS   tcrparams
#    define FORT_AUXPARAMS   auxparams
#    define FORT_ONE_N_FILL  onenfill
#    define FORT_ONE_A_FILL  oneafill
#    define FORT_ALL_N_FILL  allnfill
#    define FORT_ALL_T_FILL  alltfill
#    define FORT_ALL_A_FILL  allafill
#    define FORT_REC_N_FILL  recnfill
#    define FORT_ADVFILL     advfill
#    define FORT_XVELFILL    xvelfill
#    define FORT_YVELFILL    yvelfill
#    define FORT_ZVELFILL    zvelfill
#    define FORT_PRESFILL    presfill
#    define FORT_ENTHFILL    enthfill
#    define FORT_DIVUFILL    divufill
#    define FORT_VALGTERROR  valgterr
#    define FORT_DIFFGTERROR diffgterr
#    define FORT_VALLTERROR  vallterr
#    define FORT_KPERROR     kperror
#else
#  if  defined(BL_FORT_USE_UPPERCASE)
#    define FORT_INITPARAMS  INITPARAMS
#    define FORT_TCRPARAMS   TCRPARAMS
#    define FORT_AUXPARAMS   AUXPARAMS
#    define FORT_ONE_N_FILL  ONENFILL
#    define FORT_ONE_A_FILL  ONEAFILL
#    define FORT_ALL_N_FILL  ALLNFILL
#    define FORT_ALL_T_FILL  ALLTFILL
#    define FORT_ALL_A_FILL  ALLAFILL
#    define FORT_REC_N_FILL  RECNFILL
#    define FORT_XVELFILL    XVELFILL
#    define FORT_YVELFILL    YVELFILL
#    define FORT_ZVELFILL    ZVELFILL
#    define FORT_PRESFILL    PRESFILL
#    define FORT_ENTHFILL    ENTHFILL
#    define FORT_DIVUFILL    DIVUFILL
#    define FORT_VALGTERROR  VALGTERR
#    define FORT_VALLTERROR  VALLTERR
#    define FORT_DIFFGTERROR DIFFGTERR
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_INITPARAMS  initparams
#    define FORT_TCRPARAMS   tcrparams
#    define FORT_AUXPARAMS   auxparams
#    define FORT_ONE_N_FILL  onenfill
#    define FORT_ONE_A_FILL  oneafill
#    define FORT_ALL_N_FILL  allnfill
#    define FORT_ALL_T_FILL  alltfill
#    define FORT_ALL_A_FILL  allafill
#    define FORT_REC_N_FILL  recnfill
#    define FORT_XVELFILL    xvelfill
#    define FORT_YVELFILL    yvelfill
#    define FORT_ZVELFILL    zvelfill
#    define FORT_PRESFILL    presfill
#    define FORT_ENTHFILL    enthfill
#    define FORT_DIVUFILL    divufill
#    define FORT_VALGTERROR  valgterr
#    define FORT_VALLTERROR  vallterr
#    define FORT_DIFFGTERROR diffgterr
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_INITPARAMS  initparams_
#    define FORT_TCRPARAMS   tcrparams_
#    define FORT_AUXPARAMS   auxparams_
#    define FORT_ONE_N_FILL  onenfill_
#    define FORT_ONE_A_FILL  oneafill_
#    define FORT_ALL_N_FILL  allnfill_
#    define FORT_ALL_T_FILL  alltfill_
#    define FORT_ALL_A_FILL  allafill_
#    define FORT_REC_N_FILL  recnfill_
#    define FORT_XVELFILL    xvelfill_
#    define FORT_YVELFILL    yvelfill_
#    define FORT_ZVELFILL    zvelfill_
#    define FORT_PRESFILL    presfill_
#    define FORT_ENTHFILL    enthfill_
#    define FORT_DIVUFILL    divufill_
#    define FORT_VALGTERROR  valgterr_
#    define FORT_VALLTERROR  vallterr_
#    define FORT_DIFFGTERROR diffgterr_
#  endif

#include "ArrayLim.H"

  extern "C" {

     void FORT_INITPARAMS (int* ncomps, int* nphases, int* model_in, Real* rhoval,
			   Real* muval, int* phaseidx,Real* gravity, int* gravity_dir);

     void FORT_TCRPARAMS (int* ntracers);

     void FORT_AUXPARAMS (int* naux);

     void FORT_ONE_N_FILL (Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
			   const int* dom_lo, const int* dom_hi,
			   const Real* dx, const Real* grd_lo,
			   const Real* time, const int* bc);

     void FORT_ONE_A_FILL (Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
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

     void FORT_ALL_A_FILL (Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
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

  }
#endif
#endif /*_PROB_PM_F_H_*/
