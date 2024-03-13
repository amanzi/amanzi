/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#ifndef _tTensorDiffusion_PK_F_H_
#define _tTensorDiffusion_PK_F_H_

#include <ArrayLim.H>

#if	defined(BL_LANG_FORT)
#  define FORT_FILLCRSE  fillcrse
#  define FORT_FILLFINE  fillfine
#  define FORT_MAKEMU     makemu
#  define FORT_MAKEMU1    makemu1
#  define FORT_FILLRHS    fillrhs
#else
#  if defined(BL_FORT_USE_UPPERCASE)
#    define FORT_FILLCRSE  FILLCRSE
#    define FORT_FILLFINE  FILLFINE
#    define FORT_MAKEMU    MAKEMU
#    define FORT_MAKEMU1   MAKEMU1
#    define FORT_FILLRHS   FILLRHS
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_FILLCRSE  fillcrse_
#    define FORT_FILLFINE  fillfine_
#    define FORT_MAKEMU    makemu_
#    define FORT_MAKEMU1   makemu1_
#    define FORT_FILLRHS   fillrhs_
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_FILLCRSE  fillcrse
#    define FORT_FILLFINE  fillfine
#    define FORT_MAKEMU    makemu
#    define FORT_MAKEMU1   makemu1
#    define FORT_FILLRHS   fillrhs
#  endif

extern "C"
{
    void FORT_FILLCRSE (
	Real* crse, ARLIM_P(crse_lo), ARLIM_P(crse_hi),
	const Real* h, const int* nc,
	const int* dlo, const int* dhi, const Real* plo,
	const int* flag);
    void FORT_FILLFINE (
	Real* fine, ARLIM_P(fine_lo), ARLIM_P(fine_hi),
	const Real* h, const int* nc,
	const int* dlo, const int* dhi, const Real* plo,
	const int* lo, const int* hi,
	const int* flag);
    void FORT_FILLRHS (
	Real* fine, ARLIM_P(fine_lo), ARLIM_P(fine_hi),
	const Real* h, const int* nc,
	const int* dlo, const int* dhi, const Real* plo
	);
    void FORT_MAKEMU(
	Real*mu, ARLIM_P(mulo),ARLIM_P(muhi),
	const Real* h, int&idir,
	const Real* aT,	const Real* aL,
	const Real* u, const Real* v);
    void FORT_MAKEMU1(
	Real*mu1, ARLIM_P(mu1lo),ARLIM_P(mu1hi),
	const Real* h, int&idir,
	const Real* aT,	const Real* aL,
	const Real* u, const Real* v);
};

#endif
#endif
