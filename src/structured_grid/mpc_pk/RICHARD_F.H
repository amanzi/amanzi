#ifndef _RICHARD_F_H
#define _RICHARD_F_H

#ifdef BL_LANG_FORT
#    define FORT_CC_TO_EC_UPW     cctoecupw
#else
#  if  defined(BL_FORT_USE_UPPERCASE)
#    define FORT_CC_TO_EC_UPW     CCTOECUPW
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_CC_TO_EC_UPW     cctoecupw
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_CC_TO_EC_UPW     cctoecupw_
#  endif

#include <ArrayLim.H>

extern "C"
{
    void FORT_CC_TO_EC_UPW(Real* coef,    ARLIM_P(clo), ARLIM_P(chi),
                           const Real* l, ARLIM_P(llo), ARLIM_P(lhi),
                           const Real* v, ARLIM_P(vlo), ARLIM_P(vhi),
                           const int*  ccvlo, const int*  ccvhi,
                           const Real* mult, const int* dir,const int* nc);
}

#endif
#endif /*_MFTOWER_F_H*/
