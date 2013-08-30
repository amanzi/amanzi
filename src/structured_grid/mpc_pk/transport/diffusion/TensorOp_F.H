#ifndef _TensorOp_F_H_
#define _TensorOp_F_H_

#include <REAL.H>

#if defined(BL_LANG_FORT)
#  define FORT_APPLY            toapply
#  define FORT_FLUX             toflux
#  define FORT_GSRB             togsrb
#  define FORT_RESIDL           toresid
#  define FORT_AVERAGECC        toaveragecc
#  define FORT_AVERAGEEC        toaverageec
#  define FORT_HARMONIC_AVERAGEEC toharaverageec
#else
# if defined(BL_FORT_USE_UPPERCASE)
#  define FORT_APPLY            TOAPPLY
#  define FORT_FLUX             TOFLUX
#  define FORT_GSRB             GSRBVISC
#  define FORT_RESIDL           TORESID
#  define FORT_AVERAGECC        TOAVERAGECC
#  define FORT_AVERAGEEC        TOAVERAGEEC
#  define FORT_HARMONIC_AVERAGEEC TOHARAVERAGEEC
# elif defined(BL_FORT_USE_LOWERCASE)
#  define FORT_APPLY            toapply
#  define FORT_FLUX             toflux
#  define FORT_GSRB             togsrb
#  define FORT_RESIDL           toresid
#  define FORT_AVERAGECC        toaveragecc
#  define FORT_AVERAGEEC        toaverageec
#  define FORT_HARMONIC_AVERAGEEC toharaverageec
# elif defined(BL_FORT_USE_UNDERSCORE)
#  define FORT_APPLY            toapply_
#  define FORT_FLUX             toflux_
#  define FORT_GSRB             togsrb_
#  define FORT_RESIDL           toresid_
#  define FORT_AVERAGECC        toaveragecc_
#  define FORT_AVERAGEEC        toaverageec_
#  define FORT_HARMONIC_AVERAGEEC toharaverageec_
# endif

#include <ArrayLim.H>

extern "C"
{
  void FORT_AVERAGECC (
    Real* crseX,       ARLIM_P(crseX_lo), ARLIM_P(crseX_hi),
    const Real* fineX, ARLIM_P(fineX_lo), ARLIM_P(fineX_hi),
    const int *tlo, const int *thi, const int *nc
    );

  void FORT_RESIDL (
    Real* res      , ARLIM_P(res_lo), ARLIM_P(res_hi),
    const Real* rhs, ARLIM_P(rhs_lo), ARLIM_P(rhs_hi),
    const Real* phi, ARLIM_P(phi_lo), ARLIM_P(phi_hi),
    const int* lo, const int* hi, const int* nc
    );
  void FORT_AVERAGEEC (
    Real* crseX,       ARLIM_P(crseX_lo), ARLIM_P(crseX_hi),
    const Real* fineX, ARLIM_P(fineX_lo), ARLIM_P(fineX_hi),
    const int *tlo, const int *thi, const int *nc,
    const int *axis
    );

  void FORT_HARMONIC_AVERAGEEC (
    Real* crseX,       ARLIM_P(crseX_lo), ARLIM_P(crseX_hi),
    const Real* fineX, ARLIM_P(fineX_lo), ARLIM_P(fineX_hi),
    const int *tlo, const int *thi, const int *nc,
    const int *axis
    );

  void FORT_APPLY(
    const Real*u, ARLIM_P(ulo), ARLIM_P(uhi),
    const Real*alpha, const Real*beta,
    const Real*a, ARLIM_P(alo),   ARLIM_P(ahi),
    const Real*muX, ARLIM_P(muXlo), ARLIM_P(muXhi),
    const Real*mu1X, ARLIM_P(mu1Xlo), ARLIM_P(mu1Xhi),
    const Real*muY, ARLIM_P(muYlo), ARLIM_P(muYhi),
    const Real*mu1Y, ARLIM_P(mu1Ylo), ARLIM_P(mu1Yhi),
#if (BL_SPACEDIM > 2)
    const Real*muZ, ARLIM_P(muZlo), ARLIM_P(muZhi),
    const Real*mu1Z, ARLIM_P(mu1Zlo), ARLIM_P(mu1Zhi),
#endif
    Real*out, ARLIM_P(outlo), ARLIM_P(outhi),
    const int*maskn, ARLIM_P(masknlo), ARLIM_P(masknhi),
    const int*maske, ARLIM_P(maskelo), ARLIM_P(maskehi),
    const int*maskw, ARLIM_P(maskwlo), ARLIM_P(maskwhi),
    const int*masks, ARLIM_P(maskslo), ARLIM_P(maskshi),
#if (BL_SPACEDIM > 2)
    const int*maskt, ARLIM_P(masktlo), ARLIM_P(maskthi),
    const int*maskb, ARLIM_P(maskblo), ARLIM_P(maskbhi),
#endif
    const Real*trandern,
    ARLIM_P(trandernlo),ARLIM_P(trandernhi),
    const Real*trandere,
    ARLIM_P(tranderelo),ARLIM_P(tranderehi),
    const Real*tranderw,
    ARLIM_P(tranderwlo),ARLIM_P(tranderwhi),
    const Real*tranders,
    ARLIM_P(tranderslo),ARLIM_P(trandershi),
#if (BL_SPACEDIM > 2)
    const Real*trandert,
    ARLIM_P(trandertlo),ARLIM_P(tranderthi),
    const Real*tranderb,
    ARLIM_P(tranderblo),ARLIM_P(tranderbhi),
#endif
    const int* lo, const int* hi, const Real* h);

  void FORT_FLUX(
    const Real*u, ARLIM_P(ulo), ARLIM_P(uhi),
    const Real*muX, ARLIM_P(muXlo), ARLIM_P(muXhi),
    const Real*mu1X, ARLIM_P(mu1Xlo), ARLIM_P(mu1Xhi),
    const Real*muY, ARLIM_P(muYlo), ARLIM_P(muYhi),
    const Real*mu1Y, ARLIM_P(mu1Ylo), ARLIM_P(mu1Yhi),
#if BL_SPACEDIM>2
    const Real*muZ, ARLIM_P(muZlo), ARLIM_P(muZhi),
    const Real*mu1Z, ARLIM_P(mu1Zlo), ARLIM_P(mu1Zhi),
#endif
    Real*xflux, ARLIM_P(xfluxlo), ARLIM_P(xfluxhi),
    Real*yflux, ARLIM_P(yfluxlo), ARLIM_P(yfluxhi),
#if BL_SPACEDIM>2
    Real*zflux, ARLIM_P(zfluxlo), ARLIM_P(zfluxhi),
#endif
    const int*maskn, ARLIM_P(masknlo), ARLIM_P(masknhi),
    const int*maske, ARLIM_P(maskelo), ARLIM_P(maskehi),
    const int*maskw, ARLIM_P(maskwlo), ARLIM_P(maskwhi),
    const int*masks, ARLIM_P(maskslo), ARLIM_P(maskshi),
#if BL_SPACEDIM>2
    const int*maskt, ARLIM_P(masktlo), ARLIM_P(maskthi),
    const int*maskb, ARLIM_P(maskblo), ARLIM_P(maskbhi),
#endif
    const Real*trandern,
    ARLIM_P(trandernlo),ARLIM_P(trandernhi),
    const Real*trandere,
    ARLIM_P(tranderelo),ARLIM_P(tranderehi),
    const Real*tranderw,
    ARLIM_P(tranderwlo),ARLIM_P(tranderwhi),
    const Real*tranders,
    ARLIM_P(tranderslo),ARLIM_P(trandershi),
#if BL_SPACEDIM>2
    const Real*trandert,
    ARLIM_P(trandertlo),ARLIM_P(tranderthi),
    const Real*tranderb,
    ARLIM_P(tranderblo),ARLIM_P(tranderbhi),
#endif
    const int* lo, const int* hi, const Real* h);

  void FORT_GSRB(
    Real* u, ARLIM_P(ulo), ARLIM_P(uhi),
    const Real* rhs, ARLIM_P(rhslo), ARLIM_P(rhshi),
    Real* alpha, Real* beta,
    const Real* a, ARLIM_P(alo), ARLIM_P(ahi),
    const Real* muX, ARLIM_P(muXlo), ARLIM_P(muXhi),
    const Real* mu1X, ARLIM_P(mu1Xlo), ARLIM_P(mu1Xhi),
    const Real* muY, ARLIM_P(muYlo), ARLIM_P(muYhi),
    const Real* mu1Y, ARLIM_P(mu1Ylo), ARLIM_P(mu1Yhi),
#if BL_SPACEDIM>2
    const Real* muZ, ARLIM_P(muZlo), ARLIM_P(muZhi),
    const Real* mu1Z, ARLIM_P(mu1Zlo), ARLIM_P(mu1Zhi),
#endif
    const int* maskn, ARLIM_P(masknlo), ARLIM_P(masknhi),
    const Real* fn, ARLIM_P(fnlo), ARLIM_P(fnhi),
    const int* maske, ARLIM_P(maskelo), ARLIM_P(maskehi),
    const Real* fe, ARLIM_P(felo), ARLIM_P(fehi),
    const int* maskw, ARLIM_P(maskwlo), ARLIM_P(maskwhi),
    const Real* fw, ARLIM_P(fwlo), ARLIM_P(fwhi),
    const int* masks, ARLIM_P(maskslo), ARLIM_P(maskshi),
    const Real* fs, ARLIM_P(fslo), ARLIM_P(fshi),
#if BL_SPACEDIM>2
    const int* maskt, ARLIM_P(masktlo), ARLIM_P(maskthi),
    const Real* ft, ARLIM_P(ftlo), ARLIM_P(fthi),
    const int* maskb, ARLIM_P(maskblo), ARLIM_P(maskbhi),
    const Real* fb, ARLIM_P(fblo), ARLIM_P(fbhi),
#endif
    const Real* trandern, 
    ARLIM_P(trandernlo), ARLIM_P(trandernhi),
    const Real* trandere, 
    ARLIM_P(tranderelo), ARLIM_P(tranderehi),
    const Real* tranderw, 
    ARLIM_P(tranderwlo), ARLIM_P(tranderwhi),
    const Real* tranders, 
    ARLIM_P(tranderslo), ARLIM_P(trandershi),
#if BL_SPACEDIM>2
    const Real* trandert, 
    ARLIM_P(trandertlo), ARLIM_P(tranderthi),
    const Real* tranderb, 
    ARLIM_P(tranderblo), ARLIM_P(tranderbhi),
#endif
    const int* lo, const int* hi,
    const Real* h, int &phaseflag);
		 
}

#endif

#endif /* _TensorOp_F_H_ */
