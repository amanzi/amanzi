/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#ifndef _RichardSolver_F_H
#define _RichardSolver_F_H

#ifdef BL_LANG_FORT
#    define FORT_RS_GXPA       rsardgxpa
#    define FORT_RS_XMULTYZ    rsxmultyz
#    define FORT_RS_CTE_UPW    rsacteupw
#    define FORT_RS_RICHARDRES rsrichardres
#    define FORT_RS_SATURATEDRES rssatres
#    define FORT_RS_ALTUP      rsaltup
#else
#  if  defined(BL_FORT_USE_UPPERCASE)
#    define FORT_RS_GXPA       RSGXPA
#    define FORT_RS_XMULTYZ    RSXMULTYZ
#    define FORT_RS_CTE_UPW    RSACTEUPW
#    define FORT_RS_RICHARDRES RSRICHARDRES
#    define FORT_RS_SATURATEDRES RSSATRES
#    define FORT_RS_ALTUP      RSALTUP
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_RS_GXPA       rsardgxpa
#    define FORT_RS_XMULTYZ    rsxmultyz
#    define FORT_RS_CTE_UPW    rsacteupw
#    define FORT_RS_RICHARDRES rsrichardres
#    define FORT_RS_SATURATEDRES rssatres
#    define FORT_RS_ALTUP      rsaltup
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_RS_GXPA       rsardgxpa_
#    define FORT_RS_XMULTYZ    rsxmultyz_
#    define FORT_RS_CTE_UPW    rsacteupw_
#    define FORT_RS_RICHARDRES rsrichardres_
#    define FORT_RS_SATURATEDRES rssatres_
#    define FORT_RS_ALTUP      rsaltup_
#  endif

#include <ArrayLim.H>

extern "C"
{
  void FORT_RS_GXPA(Real* y,       ARLIM_P(ylo), ARLIM_P(yhi),
		    const Real* x, ARLIM_P(xlo), ARLIM_P(xhi),
		    const int*  ccvlo, const int*  ccvhi,
		    const Real* dx, const Real* a,
		    const int* dir,const int* nc);

  void FORT_RS_XMULTYZ(Real*       Xfab, ARLIM_P(Xlo), ARLIM_P(Xhi),
		       const Real* Yfab, ARLIM_P(Ylo), ARLIM_P(Yhi),
		       const Real* Zfab, ARLIM_P(Zlo), ARLIM_P(Zhi),
		       const int* lo, const int* hi, const int* nc);

  void FORT_RS_CTE_UPW (Real* ec,        ARLIM_P(elo), ARLIM_P(ehi),
			const Real* cc,  ARLIM_P(clo), ARLIM_P(chi),
			const Real* sgn, ARLIM_P(slo), ARLIM_P(shi),
			const int*  ccvlo, const int*  ccvhi,
			const int* dir, const int* nc,
                        const int* bclo, const int* bchi, const int* do_upwind);

  void FORT_RS_SATURATEDRES (Real* res,       ARLIM_P(reslo), ARLIM_P(reshi),
			     const Real* po,  ARLIM_P(polo),  ARLIM_P(pohi),
			     const Real* pn,  ARLIM_P(pnlo),  ARLIM_P(pnhi),
			     const Real* ss,  ARLIM_P(sslo),  ARLIM_P(sshi),
			     const Real* phi, ARLIM_P(philo), ARLIM_P(phihi),
			     const Real* src, ARLIM_P(srclo), ARLIM_P(srchi),
			     const Real* dt, const int*  lo, const int*  hi, const int* nc);

  void FORT_RS_RICHARDRES (Real* res,       ARLIM_P(reslo), ARLIM_P(reshi),
			   const Real* rso, ARLIM_P(rsolo), ARLIM_P(rsohi),
			   const Real* rsn, ARLIM_P(rsnlo), ARLIM_P(rsnhi),
			   const Real* phi, ARLIM_P(philo), ARLIM_P(phihi),
			   const Real* src, ARLIM_P(srclo), ARLIM_P(srchi),
			   const Real* dt, const int*  lo, const int*  hi, const int* nc);

  void FORT_RS_ALTUP(Real* rhosat,     ARLIM_P(rslo),  ARLIM_P(rshi),
                     Real* del_p,      ARLIM_P(dplo),  ARLIM_P(dphi),
                     const Real* al,   ARLIM_P(allo),  ARLIM_P(alhi),
                     const Real* phi,  ARLIM_P(phlo),  ARLIM_P(phhi),
                     const int*  mat,  ARLIM_P(malo),  ARLIM_P(mahi),
                     const Real* ls_factor, const Real* s_thresh,
                     const int* rmID,       const Real* time,
                     const int* lo,         const int* hi);
}

#endif
#endif /*_RichardSolver_F_H*/
