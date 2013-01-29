#ifndef _MATFILLERPC_F_H_
#define _MATFILLERPC_F_H_

#ifdef BL_LANG_FORT
#    define FORT_FILCC             filcc
#else
#  if defined(BL_FORT_USE_UPPERCASE)
#    define FORT_FILCC             FILCC
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_FILCC             filcc
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_FILCC             filcc_
#  endif

#include <ArrayLim.H>
#include "PROB_PM_F.H"

extern "C" {
    void FORT_FILCC (Real * q, ARLIM_P(q_lo), ARLIM_P(q_hi),
                     const int * domlo, const int * domhi,
                     const Real * dx_crse, const Real * xlo, 
                     const int * bc);
    
}
    
#endif
#endif /*_MATFILLERPC_F_H_*/
