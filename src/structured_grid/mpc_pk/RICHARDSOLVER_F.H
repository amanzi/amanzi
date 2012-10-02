#ifndef _RICHARDSOLVER_F_H
#define _RICHARDSOLVER_F_H

#ifdef BL_LANG_FORT
#    define FORT_RS_GXPA     rsardgxpa
#    define FORT_RS_XMULTYZ  rsxmultyz
#    define FORT_RS_CTE_UPW  rsacteupw
#    define FORT_RS_PDOTRES  rspdotres
#    define FORT_RS_ALTUP    rsaltup
#else
#  if  defined(BL_FORT_USE_UPPERCASE)
#    define FORT_RS_GXPA     RSGXPA
#    define FORT_RS_XMULTYZ  RSXMULTYZ
#    define FORT_RS_CTE_UPW  RSACTEUPW
#    define FORT_RS_PDOTRES  RSPDOTRES
#    define FORT_RS_ALTUP    RSALTUP
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_RS_GXPA     rsardgxpa
#    define FORT_RS_XMULTYZ  rsxmultyz
#    define FORT_RS_CTE_UPW  rsacteupw
#    define FORT_RS_PDOTRES  rspdotres
#    define FORT_RS_ALTUP    rsaltup
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_RS_GXPA     rsardgxpa_
#    define FORT_RS_XMULTYZ  rsxmultyz_
#    define FORT_RS_CTE_UPW  rsacteupw_
#    define FORT_RS_PDOTRES  rspdotres_
#    define FORT_RS_ALTUP    rsaltup_
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

  void FORT_RS_PDOTRES (Real* res,        ARLIM_P(reslo), ARLIM_P(reshi),
			const Real* rso,  ARLIM_P(rsnlo), ARLIM_P(rsnhi),
			const Real* rsn,  ARLIM_P(rsnlo), ARLIM_P(rsnhi),
			const Real* phio, ARLIM_P(phiolo), ARLIM_P(phiohi),
			const Real* phin, ARLIM_P(phinlo), ARLIM_P(phinhi),
			const Real* dt, const int*  lo, const int*  hi, const int* nc);

  void FORT_RS_ALTUP(Real* rhosat,   ARLIM_P(rslo), ARLIM_P(rshi),
                     Real* p,        ARLIM_P(plo),  ARLIM_P(phi),
                     const Real* dp, ARLIM_P(dplo), ARLIM_P(dphi),
                     const Real* k,  ARLIM_P(klo),  ARLIM_P(khi),
                     const Real* c,  ARLIM_P(clo),  ARLIM_P(chi),
                     const int* n_coef, const Real* ls_factor, 
                     const int* lo, const int* hi, const Real* s_thresh);
}

#endif
#endif /*_RSARDSOLVER_F_H*/
