#ifndef _MATFILLER_F_H_
#define _MATFILLER_F_H_

#ifdef BL_LANG_FORT
#    define FORT_FILCC     filcc
#    define FORT_CRSNHARM  crsnharm
#    define FORT_CRSNARITH crsnarith
#    define FORT_REFINEPC  refinepc
#else
#  if defined(BL_FORT_USE_UPPERCASE)
#    define FORT_FILCC     FILCC
#    define FORT_CRSNHARM  CRSNHARM
#    define FORT_CRSNARITH CRSNARITH
#    define FORT_REFINEPC  REFINEPC
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_FILCC     filcc
#    define FORT_CRSNHARM  crsnharm
#    define FORT_CRSNARITH crsnarith
#    define FORT_REFINEPC  refinepc
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_FILCC     filcc_
#    define FORT_CRSNHARM  crsnharm_
#    define FORT_CRSNARITH crsnarith_
#    define FORT_REFINEPC  refinepc_
#  endif

#include <ArrayLim.H>

extern "C" {
  void FORT_FILCC (Real * q, ARLIM_P(q_lo), ARLIM_P(q_hi),
                   const int * domlo, const int * domhi,
                   const Real * dx_crse, const Real * xlo,
                   const int * bc);

  void FORT_CRSNHARM (const Real* fdat, ARLIM_P(f_lo), ARLIM_P(f_hi),
                      Real*        cdat,ARLIM_P(c_lo), ARLIM_P(c_hi),
                      const int* clo, const int* chi, const int* cratio, const int* nc);

  void FORT_CRSNARITH (const Real* fdat, ARLIM_P(f_lo), ARLIM_P(f_hi),
                       Real*        cdat,ARLIM_P(c_lo), ARLIM_P(c_hi),
                       const int* clo, const int* chi, const int* cratio, const int* nc);

  void FORT_REFINEPC (const Real* cdat,ARLIM_P(c_lo), ARLIM_P(c_hi),
                      Real*       fdat, ARLIM_P(f_lo), ARLIM_P(f_hi),
                      const int* clo, const int* chi, const int* cratio, const int* nc);
    
}
    
#endif
#endif /*_MATFILLER_F_H_*/
