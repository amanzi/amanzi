#ifndef _RICHARD_F_H
#define _RICHARD_F_H

#ifdef BL_LANG_FORT
#    define FORT_RICHARD_GXPA     richardgxpa
#    define FORT_RICHARD_YMULTX   richymultx
#    define FORT_RICHARD_CTE_UPW  richacteupw
#else
#  if  defined(BL_FORT_USE_UPPERCASE)
#    define FORT_RICHARD_GXPA     RICHARDGXPA
#    define FORT_RICHARD_YMULTX   RICHYMULTX
#    define FORT_RICHARD_CTE_UPW  RICHACTEUPW
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_RICHARD_GXPA     richardgxpa
#    define FORT_RICHARD_YMULTX   richymultx
#    define FORT_RICHARD_CTE_UPW  richacteupw
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_RICHARD_GXPA     richardgxpa_
#    define FORT_RICHARD_YMULTX   richymultx_
#    define FORT_RICHARD_CTE_UPW  richacteupw_
#  endif

#include <ArrayLim.H>

extern "C"
{
    void FORT_RICHARD_GXPA(Real* y,       ARLIM_P(ylo), ARLIM_P(yhi),
                           const Real* x, ARLIM_P(xlo), ARLIM_P(xhi),
                           const int*  ccvlo, const int*  ccvhi,
                           const Real* dx, const Real* a, 
                           const int* dir,const int* nc);

    void FORT_RICHARD_YMULTX(Real*       Yfab, ARLIM_P(Ylo), ARLIM_P(Yhi),
                             const Real* Xfab, ARLIM_P(Xlo), ARLIM_P(Xhi),
                             const int* lo, const int* hi, const int* nc);

    void FORT_RICHARD_CTE_UPW (Real* ec,        ARLIM_P(elo), ARLIM_P(ehi),
                               const Real* cc,  ARLIM_P(clo), ARLIM_P(chi),
                               const Real* sgn, ARLIM_P(slo), ARLIM_P(shi),
                               const int*  ccvlo, const int*  ccvhi, 
                               const int* dir, const int* nc);
}

#endif
#endif /*_MFTOWER_F_H*/
