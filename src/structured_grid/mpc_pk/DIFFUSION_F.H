#ifndef _DIFFUSION_F_H
#define _DIFFUSION_F_H

/*
** $Id: DIFFUSION_F.H,v 1.9 2011-07-13 00:16:16 gpau Exp $
*/

#ifdef BL_LANG_FORT
#    define FORT_VISCFLUX_VC      vscflxvc
#    define FORT_VISCFLUX_CC      vscflxcc
#    define FORT_VISCSYNCFLUX     viscsyncflux
#    define FORT_HOOPSRC          hoopsrc
#    define FORT_HOOPRHS          hooprhs
#    define FORT_SETALPHA         setalpha
#    define FORT_CCTOEDGE         cctoedge
#    define FORT_DRHOG_RICHARD    drhogrichard
#    define FORT_FRHOG_RICHARD    frhogrichard
#    define FORT_DFXDOTY          dfxdoty
#    define FORT_DFSXAY           dfsxay 
#    define FORT_ARNOLDI          arnoldi
#else
#  if  defined(BL_FORT_USE_UPPERCASE)
#    define FORT_VISCFLUX_VC      VSCFLXVC
#    define FORT_VISCFLUX_CC      VSCFLXCC
#    define FORT_VISCSYNCFLUX     VISCSYNCFLUX
#    define FORT_HOOPSRC          HOOPSRC
#    define FORT_HOOPRHS          HOOPRHS
#    define FORT_SETALPHA         SETALPHA
#    define FORT_CCTOEDGE         CCTOEDGE
#    define FORT_DRHOG_RICHARD    DRHOGRICHARD
#    define FORT_FRHOG_RICHARD    FRHOGRICHARD
#    define FORT_DFXDOTY          DFXDOTY
#    define FORT_DFSXAY           DFSXAY
#    define FORT_ARNOLDI          ARNOLDI
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_VISCFLUX_VC      vscflxvc
#    define FORT_VISCFLUX_CC      vscflxcc
#    define FORT_VISCSYNCFLUX     viscsyncflux
#    define FORT_HOOPSRC          hoopsrc
#    define FORT_HOOPRHS          hooprhs
#    define FORT_SETALPHA         setalpha
#    define FORT_CCTOEDGE         cctoedge
#    define FORT_DRHOG_RICHARD    drhogrichard
#    define FORT_FRHOG_RICHARD    frhogrichard
#    define FORT_DFXDOTY          dfxdoty
#    define FORT_DFSXAY           dfsxay
#    define FORT_ARNOLDI          arnoldi
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_VISCFLUX_VC      vscflxvc_
#    define FORT_VISCFLUX_CC      vscflxcc_
#    define FORT_VISCSYNCFLUX     viscsyncflux_
#    define FORT_HOOPSRC          hoopsrc_
#    define FORT_HOOPRHS          hooprhs_
#    define FORT_SETALPHA         setalpha_
#    define FORT_CCTOEDGE         cctoedge_
#    define FORT_DRHOG_RICHARD    drhogrichard_
#    define FORT_FRHOG_RICHARD    frhogrichard_
#    define FORT_DFXDOTY          dfxdoty_
#    define FORT_DFSXAY           dfsxay_
#    define FORT_ARNOLDI          arnoldi_
#  endif

#include <ArrayLim.H>

extern "C"
{
    void FORT_VISCFLUX_VC (const Real* s_old, const Real* s_new,
                           ARLIM_P(slo), ARLIM_P(shi),
                           const int*  lo, const int*  hi,
                           const Real* flux, ARLIM_P(f_lo), ARLIM_P(f_hi),
                           const Real* area, ARLIM_P(a_lo), ARLIM_P(a_hi),
                           const Real* betan, const Real* betanp1,
                           ARLIM_P(b_lo), ARLIM_P(b_hi),
                           const Real* dx, const Real* mult, const Real* theta,
                           const int* dir);
    
    void FORT_VISCFLUX_CC (const Real* s_old, const Real* s_new,
                           ARLIM_P(slo), ARLIM_P(shi),
                           const int*  lo, const int*  hi,
                           const Real* flux, ARLIM_P(f_lo), ARLIM_P(f_hi),
                           const Real* area, ARLIM_P(a_lo), ARLIM_P(a_hi),
                           const Real* dx, const Real* mult, const Real* theta,
                           const int* dir);
#if (BL_SPACEDIM == 2)
    void FORT_VISCSYNCFLUX (const Real* ssync, ARLIM_P(slo), ARLIM_P(shi),
                            const int*  lo, const int*  hi,
                            const Real* xflux, ARLIM_P(xflo), ARLIM_P(xfhi),
                            const Real* yflux, ARLIM_P(yflo), ARLIM_P(yfhi),
                            const Real* xarea, ARLIM_P(axlo), ARLIM_P(axhi),
                            const Real* yarea, ARLIM_P(aylo), ARLIM_P(ayhi),
                            const Real* dx, const Real * mult);
    
#endif
#if (BL_SPACEDIM == 3)
    void FORT_VISCSYNCFLUX (const Real* ssync, ARLIM_P(slo), ARLIM_P(shi),
                            const int*  lo, const int*  hi,
                            const Real* xflux, ARLIM_P(xflo), ARLIM_P(xfhi),
                            const Real* yflux, ARLIM_P(yflo), ARLIM_P(yfhi),
                            const Real* zflux, ARLIM_P(zflo), ARLIM_P(zfhi),
                            const Real* xarea, ARLIM_P(axlo), ARLIM_P(axhi),
                            const Real* yarea, ARLIM_P(aylo), ARLIM_P(ayhi),
                            const Real* zarea, ARLIM_P(azlo), ARLIM_P(azhi),
                            const Real* dx, const Real * mult);
    
#endif
    void FORT_DRHOG_RICHARD(const Real* rg,  ARLIM_P(rglo), ARLIM_P(rghi),
			    const Real* cx,  ARLIM_P(cxlo), ARLIM_P(cxhi),
			    const Real* cy,  ARLIM_P(cylo), ARLIM_P(cyhi),
#if (BL_SPACEDIM == 3)
			    const Real* cz,  ARLIM_P(czlo), ARLIM_P(czhi),
#endif
			    const Real* rho, const int*  lo, const int*  hi,
			    const Real* dx,  const Real* grav, const Real* mult);

    void FORT_FRHOG_RICHARD(const Real* fx,  ARLIM_P(fxlo), ARLIM_P(fxhi),
			    const Real* cx,  ARLIM_P(cxlo), ARLIM_P(cxhi),
			    const Real* fy,  ARLIM_P(fylo), ARLIM_P(fyhi),
			    const Real* cy,  ARLIM_P(cylo), ARLIM_P(cyhi),
#if (BL_SPACEDIM == 3)
			    const Real* fz,  ARLIM_P(fzlo), ARLIM_P(fzhi),
			    const Real* cz,  ARLIM_P(czlo), ARLIM_P(czhi),
#endif
			    const Real* rho, const int*  lo, const int*  hi,
			    const Real* dx,  const Real* grav, const Real* mult);

    void FORT_HOOPSRC (ARLIM_P(grid_lo), ARLIM_P(grid_hi),
                       Real *fab, ARLIM_P(fab_lo), ARLIM_P(fab_hi),
                       const Real *u, ARLIM_P(u_lo), ARLIM_P(u_hi),
                       const Real *r, const Real* mu);
    
    void FORT_HOOPRHS (Real *fab, ARLIM_P(fab_lo), ARLIM_P(fab_hi),
                       const Real *u,   ARLIM_P(u_lo), ARLIM_P(u_hi),
                       const Real* r, const Real* b, 
                       const Real *vol, ARLIM_P(vol_lo), ARLIM_P(vol_hi));

    void FORT_SETALPHA (Real *fab, ARLIM_P(fab_lo), ARLIM_P(fab_hi),
			Real *rho, ARLIM_P(rho_lo), ARLIM_P(rho_hi),
                        const int* lo, const int* hi,
                        const Real *vfab, ARLIM_P(v_lo), ARLIM_P(v_hi));
    
    void FORT_CCTOEDGE(const int* lo, const int* hi,
        ARLIM_P(blo), ARLIM_P(bhi), const Real* betadat,
        ARLIM_P(betaxlo), ARLIM_P(betaxhi), Real* betax,
        ARLIM_P(betaylo), ARLIM_P(betayhi), Real* betay
#if (BL_SPACEDIM==3)
        ,
        ARLIM_P(betazlo), ARLIM_P(betazhi), Real* betaz
#endif
      );

    void FORT_DFXDOTY (Real *pw,
		       const Real* phi,  ARLIM_P(phi_lo),  ARLIM_P(phi_hi),
		       const int* np, const int* idxp,
		       const Real* tmpw, ARLIM_P(tmpw_lo), ARLIM_P(tmpw_hi),
		       const int* nw, const int* idxw,
		       const int* lo, const int* hi);

    void FORT_DFSXAY (Real* s,  ARLIM_P(s_lo),  ARLIM_P(s_hi),
		      const Real* x,  ARLIM_P(x_lo),  ARLIM_P(x_hi),
		      const Real* a,
		      const Real* y, ARLIM_P(y_lo), ARLIM_P(y_hi),
		      const int* lo, const int* hi,
		      const int* nc);

  void FORT_ARNOLDI(Real* soln,Real* H, Real* norm_rhs, const int *nc);
}

#endif
#endif /*_DIFFUSION_F_H*/
