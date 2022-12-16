/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#ifndef _TensorOp_F_H_
#define _TensorOp_F_H_

#include <REAL.H>

#if defined(BL_LANG_FORT)
#  define FORT_TOAPPLY            toapply
#  define FORT_TOFLUX             toflux
#  define FORT_TOGSRB             togsrb
#  define FORT_TORESIDL           toresid
#  define FORT_TOAPPLYBC          toappbc
#else
# if defined(BL_FORT_USE_UPPERCASE)
#  define FORT_TOAPPLY            TOAPPLY
#  define FORT_TOFLUX             TOFLUX
#  define FORT_TOGSRB             TOGSRB
#  define FORT_TORESIDL           TORESID
#  define FORT_TOAPPLYBC          TOAPPBC
# elif defined(BL_FORT_USE_LOWERCASE)
#  define FORT_TOAPPLY            toapply
#  define FORT_TOFLUX             toflux
#  define FORT_TOGSRB             togsrb
#  define FORT_TORESIDL           toresid
#  define FORT_TOAPPLYBC          toappbc
# elif defined(BL_FORT_USE_UNDERSCORE)
#  define FORT_TOAPPLY            toapply_
#  define FORT_TOFLUX             toflux_
#  define FORT_TOGSRB             togsrb_
#  define FORT_TORESIDL           toresid_
#  define FORT_TOAPPLYBC          toappbc_
# endif

#include <ArrayLim.H>

extern "C"
{
  void FORT_TORESIDL (
    Real* res      , ARLIM_P(res_lo), ARLIM_P(res_hi),
    const Real* rhs, ARLIM_P(rhs_lo), ARLIM_P(rhs_hi),
    const Real* phi, ARLIM_P(phi_lo), ARLIM_P(phi_hi),
    const int* lo, const int* hi, const int* nc
    );

#if (BL_SPACEDIM == 2)
  void FORT_TOAPPLY(
    const Real*u,         ARLIM_P(ulo),        ARLIM_P(uhi),
    const Real*alpha,     const Real*beta,
    const Real*a,         ARLIM_P(alo),        ARLIM_P(ahi),
    const Real*muX,       ARLIM_P(muXlo),      ARLIM_P(muXhi),
    const Real*mu1X,      ARLIM_P(mu1Xlo),     ARLIM_P(mu1Xhi),
    const Real*muY,       ARLIM_P(muYlo),      ARLIM_P(muYhi),
    const Real*mu1Y,      ARLIM_P(mu1Ylo),     ARLIM_P(mu1Yhi),
    Real*out,             ARLIM_P(outlo),      ARLIM_P(outhi),
    const int* maskn,     ARLIM_P(masknlo),    ARLIM_P(masknhi),
    const int* maske,     ARLIM_P(maskelo),    ARLIM_P(maskehi),
    const int* maskw,     ARLIM_P(maskwlo),    ARLIM_P(maskwhi),
    const int* masks,     ARLIM_P(maskslo),    ARLIM_P(maskshi),
    const Real* trandern, ARLIM_P(trandernlo), ARLIM_P(trandernhi),
    const Real* trandere, ARLIM_P(tranderelo), ARLIM_P(tranderehi),
    const Real* tranderw, ARLIM_P(tranderwlo), ARLIM_P(tranderwhi),
    const Real* tranders, ARLIM_P(tranderslo), ARLIM_P(trandershi),
    const int* lo, const int* hi, const Real* h);
#else
  void FORT_TOAPPLY(
    const Real*u,         ARLIM_P(ulo),        ARLIM_P(uhi),
    const Real*alpha,     const Real*beta,
    const Real*a,         ARLIM_P(alo),        ARLIM_P(ahi),
    const Real*muX,       ARLIM_P(muXlo),      ARLIM_P(muXhi),
    const Real*mu1X,      ARLIM_P(mu1Xlo),     ARLIM_P(mu1Xhi),
    const Real*muY,       ARLIM_P(muYlo),      ARLIM_P(muYhi),
    const Real*mu1Y,      ARLIM_P(mu1Ylo),     ARLIM_P(mu1Yhi),
    const Real*muZ,       ARLIM_P(muZlo),      ARLIM_P(muZhi),
    const Real*mu1Z,      ARLIM_P(mu1Zlo),     ARLIM_P(mu1Zhi),
    Real*out,             ARLIM_P(outlo),      ARLIM_P(outhi),
    const int* maskn,     ARLIM_P(masknlo),    ARLIM_P(masknhi),
    const int* maske,     ARLIM_P(maskelo),    ARLIM_P(maskehi),
    const int* maskw,     ARLIM_P(maskwlo),    ARLIM_P(maskwhi),
    const int* masks,     ARLIM_P(maskslo),    ARLIM_P(maskshi),
    const int* maskt,     ARLIM_P(masktlo),    ARLIM_P(maskthi),
    const int* maskb,     ARLIM_P(maskblo),    ARLIM_P(maskbhi),
    const Real* trandern, ARLIM_P(trandernlo), ARLIM_P(trandernhi),
    const Real* trandere, ARLIM_P(tranderelo), ARLIM_P(tranderehi),
    const Real* tranderw, ARLIM_P(tranderwlo), ARLIM_P(tranderwhi),
    const Real* tranders, ARLIM_P(tranderslo), ARLIM_P(trandershi),
    const Real* trandert, ARLIM_P(trandertlo), ARLIM_P(tranderthi),
    const Real* tranderb, ARLIM_P(tranderblo), ARLIM_P(tranderbhi),
    const int* lo, const int* hi, const Real* h);
#endif


#if BL_SPACEDIM==2
  void FORT_TOFLUX(
    const Real*u,        ARLIM_P(ulo),        ARLIM_P(uhi),
    const Real*muX,      ARLIM_P(muXlo),      ARLIM_P(muXhi),
    const Real*mu1X,     ARLIM_P(mu1Xlo),     ARLIM_P(mu1Xhi),
    const Real*muY,      ARLIM_P(muYlo),      ARLIM_P(muYhi),
    const Real*mu1Y,     ARLIM_P(mu1Ylo),     ARLIM_P(mu1Yhi),
    Real*xflux,          ARLIM_P(xfluxlo),    ARLIM_P(xfluxhi),
    Real*yflux,          ARLIM_P(yfluxlo),    ARLIM_P(yfluxhi),
    const int*maskn,     ARLIM_P(masknlo),    ARLIM_P(masknhi),
    const int*maske,     ARLIM_P(maskelo),    ARLIM_P(maskehi),
    const int*maskw,     ARLIM_P(maskwlo),    ARLIM_P(maskwhi),
    const int*masks,     ARLIM_P(maskslo),    ARLIM_P(maskshi),
    const Real*trandern, ARLIM_P(trandernlo), ARLIM_P(trandernhi),
    const Real*trandere, ARLIM_P(tranderelo), ARLIM_P(tranderehi),
    const Real*tranderw, ARLIM_P(tranderwlo), ARLIM_P(tranderwhi),
    const Real*tranders, ARLIM_P(tranderslo), ARLIM_P(trandershi),
    const int* lo, const int* hi, const Real* h);
#else
  void FORT_TOFLUX(
    const Real*u,        ARLIM_P(ulo),        ARLIM_P(uhi),
    const Real*muX,      ARLIM_P(muXlo),      ARLIM_P(muXhi),
    const Real*mu1X,     ARLIM_P(mu1Xlo),     ARLIM_P(mu1Xhi),
    const Real*muY,      ARLIM_P(muYlo),      ARLIM_P(muYhi),
    const Real*mu1Y,     ARLIM_P(mu1Ylo),     ARLIM_P(mu1Yhi),
    const Real*muZ,      ARLIM_P(muZlo),      ARLIM_P(muZhi),
    const Real*mu1Z,     ARLIM_P(mu1Zlo),     ARLIM_P(mu1Zhi),
    Real*xflux,          ARLIM_P(xfluxlo),    ARLIM_P(xfluxhi),
    Real*yflux,          ARLIM_P(yfluxlo),    ARLIM_P(yfluxhi),
    Real*zflux,          ARLIM_P(zfluxlo),    ARLIM_P(zfluxhi),
    const int*maskn,     ARLIM_P(masknlo),    ARLIM_P(masknhi),
    const int*maske,     ARLIM_P(maskelo),    ARLIM_P(maskehi),
    const int*maskw,     ARLIM_P(maskwlo),    ARLIM_P(maskwhi),
    const int*masks,     ARLIM_P(maskslo),    ARLIM_P(maskshi),
    const int*maskt,     ARLIM_P(masktlo),    ARLIM_P(maskthi),
    const int*maskb,     ARLIM_P(maskblo),    ARLIM_P(maskbhi),
    const Real*trandern, ARLIM_P(trandernlo), ARLIM_P(trandernhi),
    const Real*trandere, ARLIM_P(tranderelo), ARLIM_P(tranderehi),
    const Real*tranderw, ARLIM_P(tranderwlo), ARLIM_P(tranderwhi),
    const Real*tranders, ARLIM_P(tranderslo), ARLIM_P(trandershi),
    const Real*trandert, ARLIM_P(trandertlo), ARLIM_P(tranderthi),
    const Real*tranderb, ARLIM_P(tranderblo), ARLIM_P(tranderbhi),
    const int* lo, const int* hi, const Real* h);
#endif


#if BL_SPACEDIM==2
  void FORT_TOGSRB(
    Real* u,              ARLIM_P(ulo),        ARLIM_P(uhi),
    const Real* rhs,      ARLIM_P(rhslo),      ARLIM_P(rhshi),  Real* alpha, Real* beta,
    const Real* a,        ARLIM_P(alo),        ARLIM_P(ahi),
    const Real* muX,      ARLIM_P(muXlo),      ARLIM_P(muXhi),
    const Real* mu1X,     ARLIM_P(mu1Xlo),     ARLIM_P(mu1Xhi),
    const Real* muY,      ARLIM_P(muYlo),      ARLIM_P(muYhi),
    const Real* mu1Y,     ARLIM_P(mu1Ylo),     ARLIM_P(mu1Yhi),
    const int* maskn,     ARLIM_P(masknlo),    ARLIM_P(masknhi),
    const Real* fn,       ARLIM_P(fnlo),       ARLIM_P(fnhi),
    const Real* ftn,      ARLIM_P(ftnlo),      ARLIM_P(ftnhi),
    const int* maske,     ARLIM_P(maskelo),    ARLIM_P(maskehi),
    const Real* fe,       ARLIM_P(felo),       ARLIM_P(fehi),
    const Real* fte,      ARLIM_P(ftelo),      ARLIM_P(ftehi),
    const int* maskw,     ARLIM_P(maskwlo),    ARLIM_P(maskwhi),
    const Real* fw,       ARLIM_P(fwlo),       ARLIM_P(fwhi),
    const Real* ftw,      ARLIM_P(ftwlo),      ARLIM_P(ftwhi),
    const int* masks,     ARLIM_P(maskslo),    ARLIM_P(maskshi),
    const Real* fs,       ARLIM_P(fslo),       ARLIM_P(fshi),
    const Real* fts,      ARLIM_P(ftslo),      ARLIM_P(ftshi),
    const Real* trandern, ARLIM_P(trandernlo), ARLIM_P(trandernhi),
    const Real* trandere, ARLIM_P(tranderelo), ARLIM_P(tranderehi),
    const Real* tranderw, ARLIM_P(tranderwlo), ARLIM_P(tranderwhi),
    const Real* tranders, ARLIM_P(tranderslo), ARLIM_P(trandershi),
    const int* lo, const int* hi, const Real* h, int &phaseflag);
#else
  void FORT_TOGSRB(
    Real* u,              ARLIM_P(ulo),        ARLIM_P(uhi),
    const Real* rhs,      ARLIM_P(rhslo),      ARLIM_P(rhshi),  Real* alpha, Real* beta,
    const Real* a,        ARLIM_P(alo),        ARLIM_P(ahi),
    const Real* muX,      ARLIM_P(muXlo),      ARLIM_P(muXhi),
    const Real* mu1X,     ARLIM_P(mu1Xlo),     ARLIM_P(mu1Xhi),
    const Real* muY,      ARLIM_P(muYlo),      ARLIM_P(muYhi),
    const Real* mu1Y,     ARLIM_P(mu1Ylo),     ARLIM_P(mu1Yhi),
    const Real* muZ,      ARLIM_P(muZlo),      ARLIM_P(muZhi),
    const Real* mu1Z,     ARLIM_P(mu1Zlo),     ARLIM_P(mu1Zhi),
    const int* maskn,     ARLIM_P(masknlo),    ARLIM_P(masknhi),
    const Real* fn,       ARLIM_P(fnlo),       ARLIM_P(fnhi),
    const Real* ftn,      ARLIM_P(ftnlo),      ARLIM_P(ftnhi),
    const int* maske,     ARLIM_P(maskelo),    ARLIM_P(maskehi),
    const Real* fe,       ARLIM_P(felo),       ARLIM_P(fehi),
    const Real* fte,      ARLIM_P(ftelo),      ARLIM_P(ftehi),
    const int* maskw,     ARLIM_P(maskwlo),    ARLIM_P(maskwhi),
    const Real* fw,       ARLIM_P(fwlo),       ARLIM_P(fwhi),
    const Real* ftw,      ARLIM_P(ftwlo),      ARLIM_P(ftwhi),
    const int* masks,     ARLIM_P(maskslo),    ARLIM_P(maskshi),
    const Real* fs,       ARLIM_P(fslo),       ARLIM_P(fshi),
    const Real* fts,      ARLIM_P(ftslo),      ARLIM_P(ftshi),
    const int* maskt,     ARLIM_P(masktlo),    ARLIM_P(maskthi),
    const Real* ft,       ARLIM_P(ftlo),       ARLIM_P(fthi),
    const Real* ftt,      ARLIM_P(fttlo),      ARLIM_P(ftthi),
    const int* maskb,     ARLIM_P(maskblo),    ARLIM_P(maskbhi),
    const Real* fb,       ARLIM_P(fblo),       ARLIM_P(fbhi),
    const Real* ftb,      ARLIM_P(ftblo),      ARLIM_P(ftbhi),
    const Real* trandern, ARLIM_P(trandernlo), ARLIM_P(trandernhi),
    const Real* trandere, ARLIM_P(tranderelo), ARLIM_P(tranderehi),
    const Real* tranderw, ARLIM_P(tranderwlo), ARLIM_P(tranderwhi),
    const Real* tranders, ARLIM_P(tranderslo), ARLIM_P(trandershi),
    const Real* trandert, ARLIM_P(trandertlo), ARLIM_P(tranderthi),
    const Real* tranderb, ARLIM_P(tranderblo), ARLIM_P(tranderbhi),
    const int* lo, const int* hi, const Real* h, int &phaseflag);
#endif


    void FORT_TOAPPLYBC(
		    const int *flagden,  // 1 if want values in den
		    const int *flagbc,   // 1 for inhomogeneous
		    const int *maxorder, // maxorder of interpolant
		    Real *phi, ARLIM_P(phi_lo), ARLIM_P(phi_hi),
		    const int *cdir,
		    const int *bct,
		    const Real *bcl,
		    const Real *bcval, ARLIM_P(bcval_lo), ARLIM_P(bcval_hi),
#if BL_SPACEDIM==2
		    const int *mask,   ARLIM_P(mask_lo),  ARLIM_P(mask_hi),
		    const int *maskp,  ARLIM_P(maskp_lo),  ARLIM_P(maskp_hi),
		    const int *maskm,  ARLIM_P(maskm_lo),  ARLIM_P(maskm_hi),
#elif BL_SPACEDIM==3
		    const int *maskn,  ARLIM_P(maskn_lo),  ARLIM_P(maskn_hi),
		    const int *maske,  ARLIM_P(maske_lo),  ARLIM_P(maske_hi),
		    const int *maskw,  ARLIM_P(maskw_lo),  ARLIM_P(maskw_hi),
		    const int *masks,  ARLIM_P(masks_lo),  ARLIM_P(masks_hi),
		    const int *maskt,  ARLIM_P(maskt_lo),  ARLIM_P(maskt_hi),
		    const int *maskb,  ARLIM_P(maskb_lo),  ARLIM_P(maskb_hi),
#endif
		    Real *den,         ARLIM_P(den_lo),   ARLIM_P(den_hi),
		    Real *dent,        ARLIM_P(dent_lo),  ARLIM_P(dent_hi),
		    const Real *exttd, ARLIM_P(exttd_lo), ARLIM_P(exttd_hi),
		    Real *tander,      ARLIM_P(tander_lo),ARLIM_P(tander_hi),
		    const int *lo, const int *hi,
		    const int *nc,
		    const Real *h
		    );
}

#endif

#endif /* _TensorOp_F_H_ */
