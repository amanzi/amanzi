#ifndef _MACOUTFLOWBC_F_H
#define _MACOUTFLOWBC_F_H

#ifdef BL_LANG_FORT
#    define FORT_RHOGBC      rhogbc
#else
#  if defined(BL_FORT_USE_UPPERCASE)
#    define FORT_RHOGBC      RHOGBC
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_RHOGBC      rhogbc
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_RHOGBC      rhogbc_
#  endif

#include <ArrayLim.H>

extern "C" {

    void FORT_RHOGBC(const Real* rho, ARLIM_P(rho_lo), ARLIM_P(rho_hi),
                     const Real* phi, ARLIM_P(phi_lo), ARLIM_P(phi_hi),
                     int * face, const Real* gravity, const Real* dx,
                     const int* domlo, const int* domhi,
                     const int* lo_bc, const int* hi_bc);
}

#endif
#endif 
