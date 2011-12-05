#ifndef _DERIVE_F_H
#define _DERIVE_F_H

#ifdef BL_LANG_FORT
#    define FORT_DERGRDPX    dergrdpx
#    define FORT_DERGRDPY    dergrdpy
#    define FORT_DERGRDPZ    dergrdpz
#    define FORT_DERGRDN     dergrdn
#else
#  if defined(BL_FORT_USE_UPPERCASE)
#    define FORT_DERGRDPX    DERGRDPX
#    define FORT_DERGRDPY    DERGRDPY
#    define FORT_DERGRDPZ    DERGRDPZ
#    define FORT_DERGRDN     DERGRDN
#  elif defined (BL_FORT_USE_LOWERCASE)
#    define FORT_DERGRDPX    dergrdpx
#    define FORT_DERGRDPY    dergrdpy
#    define FORT_DERGRDPZ    dergrdpz
#    define FORT_DERGRDN     dergrdn

#  elif defined (BL_FORT_USE_UNDERSCORE)
#    define FORT_DERGRDPX    dergrdpx_
#    define FORT_DERGRDPY    dergrdpy_
#    define FORT_DERGRDPZ    dergrdpz_
#    define FORT_DERGRDN     dergrdn_
#  endif

#include <ArrayLim.H>

extern "C" {

     void FORT_DERGRDPX   (Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const Real* delta, const Real* xlo,
			   const Real* time, const Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void FORT_DERGRDPY   (Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const Real* delta, const Real* xlo,
			   const Real* time, const Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void FORT_DERGRDPZ   (Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const Real* delta, const Real* xlo,
			   const Real* time, const Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void FORT_DERGRDN    (Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const Real* delta, const Real* xlo,
			   const Real* time, const Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no);
}
#endif
#endif
