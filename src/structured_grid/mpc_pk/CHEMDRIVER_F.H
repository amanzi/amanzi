#ifndef _CHEMDRIVER_F_H_
#define _CHEMDRIVER_F_H_

#ifdef BL_LANG_FORT
#    define FORT_CHEM_TOUGH      chem_tough
#    define FORT_INIT_CHEM_TOUGH initchem_tough
#    define FORT_SAT_TOUGH       sat_tough
#    define FORT_INIT_SECONDARY_TOUGH initsecondary_tough
#else
#  if defined(BL_FORT_USE_UPPERCASE)
#    define FORT_CHEM_TOUGH        CHEM_TOUGH
#    define FORT_INIT_CHEM_TOUGH   INITCHEM_TOUGH
#    define FORT_SAT_TOUGH         SAT_TOUGH
#    define FORT_INIT_SECONDARY_TOUGH INITSECONDARY_TOUGH
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_CHEM_TOUGH        chem_tough
#    define FORT_INIT_CHEM_TOUGH   initchem_tough
#    define FORT_SAT_TOUGH         sat_tough
#    define FORT_INIT_SECONDARY_TOUGH initsecondary_tough
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_CHEM_TOUGH        chem_tough_
#    define FORT_INIT_CHEM_TOUGH   initchem_tough_
#    define FORT_SAT_TOUGH         sat_tough_
#    define FORT_INIT_SECONDARY_TOUGH initsecondary_tough_
#  endif

#include "ArrayLim.H"

  extern "C" {

    void FORT_CHEM_TOUGH (const Real* s, ARLIM_P(slo), ARLIM_P(shi),  
			  const Real* press, ARLIM_P(p_lo), ARLIM_P(p_hi),
			  const Real* cr, ARLIM_P(c_lo), ARLIM_P(c_hi),
			  const Real* fr, ARLIM_P(f_lo), ARLIM_P(f_hi),
			  int* nc, const int* lo, const int* hi, Real* dt,
			  int* nd, int* idx_cr, const Real* mweight,
			  const Real* kappa, ARLIM_P(k_lo), ARLIM_P(k_hi),
			  const Real* phi, ARLIM_P(ph_lo), ARLIM_P(ph_hi));

    void FORT_INIT_SECONDARY_TOUGH(
			  const Real* s, ARLIM_P(slo), ARLIM_P(shi),  
			  const Real* press, ARLIM_P(p_lo), ARLIM_P(p_hi),
			  const Real* cr, ARLIM_P(c_lo), ARLIM_P(c_hi),
			  int* nc, const int* lo, const int* hi,
			  int* nd, int* idx_cr, const Real* mweight,
			  const Real* kappa, ARLIM_P(k_lo), ARLIM_P(k_hi),
			  const Real* phi, ARLIM_P(ph_lo), ARLIM_P(ph_hi));

    void FORT_SAT_TOUGH (const Real* s, ARLIM_P(slo), ARLIM_P(shi),
			 const Real* st, ARLIM_P(stlo), ARLIM_P(sthi),
			 const Real* press, ARLIM_P(p_lo), ARLIM_P(p_hi),
			 const Real* cr, ARLIM_P(c_lo), ARLIM_P(c_hi), 
			 int* nc, const int* lo, const int* hi, 
			 int* nd, int* idx_cr, const Real* mweight,
			 const Real* kappa, ARLIM_P(k_lo), ARLIM_P(k_hi),
			 const Real* phi, ARLIM_P(ph_lo), ARLIM_P(ph_hi));

    void FORT_INIT_CHEM_TOUGH();
  }
#endif
#endif /*_CHEMDRIVER_F_H_*/
