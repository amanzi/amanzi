#ifndef _MFTOWER_F_H
#define _MFTOWER_F_H

#ifdef BL_LANG_FORT
#    define FORT_CC_TO_EC_GRAD    cctoecgrad
#    define FORT_EC_TO_CC_DIV     ectoccdiv
#else
#  if  defined(BL_FORT_USE_UPPERCASE)
#    define FORT_CC_TO_EC_GRAD    CCTOECGRAD
#    define FORT_EC_TO_CC_DIV     ECTOCCDIV
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_CC_TO_EC_GRAD    cctoecgrad
#    define FORT_EC_TO_CC_DIV     ectoccdiv
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_CC_TO_EC_GRAD    cctoecgrad_
#    define FORT_EC_TO_CC_DIV     ectoccdiv_
#  endif

#include <ArrayLim.H>

extern "C"
{
    void FORT_CC_TO_EC_GRAD (Real* ecg,      ARLIM_P(elo), ARLIM_P(ehi),
                             const Real* cc, ARLIM_P(clo), ARLIM_P(chi),
                             const int*  ccvlo, const int*  ccvhi,
                             const Real* dx, const Real* mult, const int* dir,const int* nc);

    void FORT_EC_TO_CC_DIV (Real* ccg,      ARLIM_P(clo), ARLIM_P(chi),
                            const Real* ec, ARLIM_P(elo), ARLIM_P(ehi),
                            const int*  ccvlo, const int*  ccvhi,
                            const Real* dx, const Real* mult, const int* dir,const int* nc);

}

#endif
#endif /*_MFTOWER_F_H*/
