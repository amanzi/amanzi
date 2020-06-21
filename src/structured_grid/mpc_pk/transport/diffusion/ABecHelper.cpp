#include <ABecHelper.H>
#include <LO_F.H>
#include <ABec_F.H>

void 
ABecHelper::compFlux (D_DECL(MultiFab &xflux, MultiFab &yflux, MultiFab &zflux),
                      MultiFab& in, const BC_Mode& bc_mode,
                      int src_comp, int dst_comp, int num_comp, int bnd_comp) 
{
  compFlux(D_DECL(xflux, yflux, zflux), in, true, bc_mode, src_comp, dst_comp, num_comp, bnd_comp);
}
  
void 
ABecHelper::compFlux (D_DECL(MultiFab &xflux, MultiFab &yflux, MultiFab &zflux),
                    MultiFab& in, bool do_ApplyBC, const BC_Mode& bc_mode,
                    int src_comp, int dst_comp, int num_comp, int bnd_comp)
{
  int bndryComp = (bnd_comp < 0 ? default_bndryComp : bnd_comp);
  int alphaComp = default_betaComp;
  int betaComp = default_betaComp;

  const int level = 0;
  BL_ASSERT(num_comp==1);
  
  if (do_ApplyBC)
      applyBC(in,src_comp,num_comp,level,bc_mode,bnd_comp);
  
  const MultiFab& a = aCoefficients(level);
  
  D_TERM(const MultiFab& bX = bCoefficients(0,level);,
         const MultiFab& bY = bCoefficients(1,level);,
         const MultiFab& bZ = bCoefficients(2,level););
  
  Real _alpha = get_alpha();
  Real _beta = get_beta();

  const bool tiling = true;

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter inmfi(in,tiling); inmfi.isValid(); ++inmfi)
  {
    D_TERM(const Box& xbx   = inmfi.nodaltilebox(0);,
	   const Box& ybx   = inmfi.nodaltilebox(1);,
	   const Box& zbx   = inmfi.nodaltilebox(2););

    FArrayBox& infab = in[inmfi];
    
    D_TERM(const FArrayBox& bxfab = bX[inmfi];,
           const FArrayBox& byfab = bY[inmfi];,
           const FArrayBox& bzfab = bZ[inmfi];);
    
    D_TERM(FArrayBox& xfluxfab = xflux[inmfi];,
           FArrayBox& yfluxfab = yflux[inmfi];,
           FArrayBox& zfluxfab = zflux[inmfi];);
    
    FORT_FLUX(infab.dataPtr(src_comp),
              ARLIM(infab.loVect()), ARLIM(infab.hiVect()),
              &_alpha, &_beta, a[inmfi].dataPtr(alphaComp),
              ARLIM(a[inmfi].loVect()), ARLIM(a[inmfi].hiVect()),
              bxfab.dataPtr(betaComp), 
              ARLIM(bxfab.loVect()), ARLIM(bxfab.hiVect()),
#if (BL_SPACEDIM >= 2)
              byfab.dataPtr(betaComp), 
              ARLIM(byfab.loVect()), ARLIM(byfab.hiVect()),
#if (BL_SPACEDIM == 3)
              bzfab.dataPtr(betaComp), 
              ARLIM(bzfab.loVect()), ARLIM(bzfab.hiVect()),
#endif
#endif
	      xbx.loVect(), xbx.hiVect(), 
#if (BL_SPACEDIM >= 2)
	      ybx.loVect(), ybx.hiVect(), 
#if (BL_SPACEDIM == 3)
	      zbx.loVect(), zbx.hiVect(), 
#endif
#endif
              &num_comp,
              h[level],
              xfluxfab.dataPtr(dst_comp),
              ARLIM(xfluxfab.loVect()), ARLIM(xfluxfab.hiVect())
#if (BL_SPACEDIM >= 2)
              ,yfluxfab.dataPtr(dst_comp),
              ARLIM(yfluxfab.loVect()), ARLIM(yfluxfab.hiVect())
#endif
#if (BL_SPACEDIM == 3)
              ,zfluxfab.dataPtr(dst_comp),
              ARLIM(zfluxfab.loVect()), ARLIM(zfluxfab.hiVect())
#endif
      );
  }
}

void
ABecHelper::apply (MultiFab&      out,
                   MultiFab&      in,
                   int            level,
                   LinOp::BC_Mode bc_mode,
                   bool           local,
                   int            src_comp,
                   int            dst_comp,
                   int            num_comp,
                   int            bndry_comp)
{
  int bndryComp = (bndry_comp < 0 ? default_bndryComp : bndry_comp);
  applyBC(in,src_comp,num_comp,level,bc_mode,local,bndryComp);
  Fapply(out,dst_comp,in,src_comp,num_comp,level);
}

void
ABecHelper::applyBC (MultiFab&      inout,
                     int            src_comp,
                     int            num_comp,
                     int            level,
                     LinOp::BC_Mode bc_mode,
                     bool           local,
                     int            bndry_comp)
{
  int bndryComp = (bndry_comp < 0 ? default_bndryComp : bndry_comp);
  LinOp::applyBC(inout,src_comp,num_comp,level,bc_mode,local,bndryComp);
}

void
ABecHelper::residual (MultiFab&       residL,
                      const MultiFab& rhsL,
                      MultiFab&       solnL,
                      int             level,
                      LinOp::BC_Mode  bc_mode,
                      bool            local)
{
  int srcComp = 0;
  int destComp = 0;
  int numComp = 1;
  int rhsComp = 0;
  int bndryComp = default_bndryComp;
  apply(residL, solnL, level, bc_mode, local, srcComp, destComp, numComp, bndryComp);
  
  const bool tiling = true;

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter solnLmfi(solnL,tiling); solnLmfi.isValid(); ++solnLmfi)
  {
    const int nc = residL.nComp();
    //
    // Only single-component solves supported (verified) by this class.
    //
    BL_ASSERT(nc == 1);
    
    const Box& tbx = solnLmfi.tilebox();
    
    BL_ASSERT(gbox[level][solnLmfi.index()] == solnLmfi.validbox());
    
    FArrayBox& residfab = residL[solnLmfi];
    const FArrayBox& rhsfab = rhsL[solnLmfi];
    
    FORT_RESIDL(
      residfab.dataPtr(destComp), 
      ARLIM(residfab.loVect()), ARLIM(residfab.hiVect()),
      rhsfab.dataPtr(rhsComp), 
      ARLIM(rhsfab.loVect()), ARLIM(rhsfab.hiVect()),
      residfab.dataPtr(), 
      ARLIM(residfab.loVect()), ARLIM(residfab.hiVect()),
      tbx.loVect(), tbx.hiVect(), &nc);
  }
}

void 
ABecHelper::smooth (MultiFab&       solnL,
                    const MultiFab& rhsL,
                    int             level,
                    LinOp::BC_Mode  bc_mode)
{
  for (int redBlackFlag = 0; redBlackFlag < 2; redBlackFlag++) {
    applyBC(solnL, 0, 1, level, bc_mode, false, default_bndryComp);
    Fsmooth(solnL, rhsL, level, redBlackFlag);
  }
}

void
ABecHelper::jacobi_smooth (MultiFab&       solnL,
                           const MultiFab& rhsL,
                           int             level,
                           LinOp::BC_Mode  bc_mode)
{
  Fsmooth_jacobi(solnL, rhsL, level);
}

Real
ABecHelper::norm (int nm, int level, const bool local)
{
  int alphaComp = default_alphaComp;
  int betaComp = default_betaComp;
  int bndryComp = default_bndryComp;

  BL_ASSERT(nm == 0);
  const MultiFab& a   = aCoefficients(level);

  D_TERM(const MultiFab& bX  = bCoefficients(0,level);,
         const MultiFab& bY  = bCoefficients(1,level);,
         const MultiFab& bZ  = bCoefficients(2,level););

  //const int nc = a.nComp(); // FIXME: This LinOp only really support single-component
  const int nc = 1;
  Real res = 0.0;
  Real _alpha = get_alpha();
  Real _beta = get_beta();

  const bool tiling = true;

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    Real pres = 0.0;

    for (MFIter amfi(a,tiling); amfi.isValid(); ++amfi)
    {

      Real tres;

      const Box&       tbx  = amfi.tilebox();
      const FArrayBox& afab = a[amfi];

      D_TERM(const FArrayBox& bxfab = bX[amfi];,
	     const FArrayBox& byfab = bY[amfi];,
	     const FArrayBox& bzfab = bZ[amfi];);
 
#if (BL_SPACEDIM==2)
      FORT_NORMA(&tres,
		 &_alpha, &_beta,
		 afab.dataPtr(alphaComp), ARLIM(afab.loVect()),  ARLIM(afab.hiVect()),
		 bxfab.dataPtr(betaComp), ARLIM(bxfab.loVect()), ARLIM(bxfab.hiVect()),
		 byfab.dataPtr(betaComp), ARLIM(byfab.loVect()), ARLIM(byfab.hiVect()),
		 tbx.loVect(), tbx.hiVect(), &nc,
		 h[level]);

#elif (BL_SPACEDIM==3)
      FORT_NORMA(&tres,
		 &_alpha, &_beta,
		 afab.dataPtr(alphaComp),  ARLIM(afab.loVect()), ARLIM(afab.hiVect()),
		 bxfab.dataPtr(betaComp), ARLIM(bxfab.loVect()), ARLIM(bxfab.hiVect()),
		 byfab.dataPtr(betaComp), ARLIM(byfab.loVect()), ARLIM(byfab.hiVect()),
		 bzfab.dataPtr(betaComp), ARLIM(bzfab.loVect()), ARLIM(bzfab.hiVect()),
		 tbx.loVect(), tbx.hiVect(), &nc,
		 h[level]);
#endif
      pres = std::max(pres, tres);
    }
    
#ifdef _OPENMP
#pragma omp critical(ABec_norm)
#endif
    { // reduction(max:) only became available in OpenMP 3.1, so we still use omp critical here
      res = std::max(res, pres);
    }
  }

  if (!local)
    ParallelDescriptor::ReduceRealMax(res);
  return res;
}

void 
ABecHelper::Fapply (MultiFab&       out,
                    const MultiFab& in,
                    int             level)
{
  int num_comp = 1;
  int src_comp = 0;
  int dst_comp = 0;
  Fapply(out,dst_comp,in,src_comp,num_comp,level);
}

void
ABecHelper::Fapply (MultiFab&       out,
                    int             dst_comp,
                    const MultiFab& in,
                    int             src_comp,
                    int             num_comp,
                    int             level)
{
  int alphaComp = default_alphaComp;
  int betaComp = default_betaComp;
  int bndryComp = default_bndryComp;

  BL_ASSERT(out.nComp()>=dst_comp+num_comp);
  BL_ASSERT(in.nComp()>=src_comp+num_comp);

  const MultiFab& a   = aCoefficients(level);

  D_TERM(const MultiFab& bX  = bCoefficients(0,level);,
         const MultiFab& bY  = bCoefficients(1,level);,
         const MultiFab& bZ  = bCoefficients(2,level););

  const bool tiling = true;

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(out,tiling); mfi.isValid(); ++mfi)
  {
    const Box&       tbx  = mfi.tilebox();
    FArrayBox&       outfab = out[mfi];
    const FArrayBox& infab = in[mfi];
    const FArrayBox& afab = a[mfi];

    D_TERM(const FArrayBox& bxfab = bX[mfi];,
           const FArrayBox& byfab = bY[mfi];,
           const FArrayBox& bzfab = bZ[mfi];);

    Real _alpha = get_alpha();
    Real _beta = get_beta();

#if (BL_SPACEDIM == 2)
    FORT_ADOTX(outfab.dataPtr(dst_comp),
               ARLIM(outfab.loVect()),ARLIM(outfab.hiVect()),
               infab.dataPtr(src_comp),
               ARLIM(infab.loVect()), ARLIM(infab.hiVect()),
               &_alpha, &_beta, afab.dataPtr(alphaComp), 
               ARLIM(afab.loVect()), ARLIM(afab.hiVect()),
               bxfab.dataPtr(betaComp), 
               ARLIM(bxfab.loVect()), ARLIM(bxfab.hiVect()),
               byfab.dataPtr(betaComp), 
               ARLIM(byfab.loVect()), ARLIM(byfab.hiVect()),
               tbx.loVect(), tbx.hiVect(), &num_comp,
               h[level]);
#endif
#if (BL_SPACEDIM ==3)
    FORT_ADOTX(outfab.dataPtr(dst_comp),
               ARLIM(outfab.loVect()),ARLIM(outfab.hiVect()),
               infab.dataPtr(src_comp),
               ARLIM(infab.loVect()), ARLIM(infab.hiVect()),
               &_alpha, &_beta, afab.dataPtr(alphaComp), 
               ARLIM(afab.loVect()), ARLIM(afab.hiVect()),
               bxfab.dataPtr(betaComp), 
               ARLIM(bxfab.loVect()), ARLIM(bxfab.hiVect()),
               byfab.dataPtr(betaComp), 
               ARLIM(byfab.loVect()), ARLIM(byfab.hiVect()),
               bzfab.dataPtr(betaComp), 
               ARLIM(bzfab.loVect()), ARLIM(bzfab.hiVect()),
               tbx.loVect(), tbx.hiVect(), &num_comp,
               h[level]);
#endif
  }
}

void
ABecHelper::Fsmooth (MultiFab&       solnL,
                     const MultiFab& rhsL,
                     int             level,
                     int             redBlackFlag)
{
  int alphaComp = default_alphaComp;
  int betaComp = default_betaComp;
  int bndryComp = default_bndryComp;
  int rhsComp = 0;
  int solnComp = 0;

  OrientationIter oitr;

  const FabSet& f0 = (*undrrelxr[level])[oitr()]; oitr++;
  const FabSet& f1 = (*undrrelxr[level])[oitr()]; oitr++;
  const FabSet& f2 = (*undrrelxr[level])[oitr()]; oitr++;
  const FabSet& f3 = (*undrrelxr[level])[oitr()]; oitr++;
#if (BL_SPACEDIM > 2)
  const FabSet& f4 = (*undrrelxr[level])[oitr()]; oitr++;
  const FabSet& f5 = (*undrrelxr[level])[oitr()]; oitr++;
#endif    
  const MultiFab& a = aCoefficients(level);

  D_TERM(const MultiFab& bX = bCoefficients(0,level);,
         const MultiFab& bY = bCoefficients(1,level);,
         const MultiFab& bZ = bCoefficients(2,level););

  //const int nc = solnL.nComp(); // FIXME: This LinOp only really supports single-component
  const int nc = 1;

  const bool tiling = true;

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter solnLmfi(solnL,tiling); solnLmfi.isValid(); ++solnLmfi)
  {
    OrientationIter oitr;

    const int gn = solnLmfi.index();

    const Mask& m0 = maskvals[level][oitr()][solnLmfi]; oitr++;
    const Mask& m1 = maskvals[level][oitr()][solnLmfi]; oitr++;
    const Mask& m2 = maskvals[level][oitr()][solnLmfi]; oitr++;
    const Mask& m3 = maskvals[level][oitr()][solnLmfi]; oitr++;
#if (BL_SPACEDIM > 2)
    const Mask& m4 = maskvals[level][oitr()][solnLmfi]; oitr++;
    const Mask& m5 = maskvals[level][oitr()][solnLmfi]; oitr++;
#endif
    const Box&       tbx     = solnLmfi.tilebox();
    const Box&       vbx     = solnLmfi.validbox();
    FArrayBox&       solnfab = solnL[gn];
    const FArrayBox& rhsfab  = rhsL[gn];
    const FArrayBox& afab    = a[gn];

    D_TERM(const FArrayBox& bxfab = bX[gn];,
           const FArrayBox& byfab = bY[gn];,
           const FArrayBox& bzfab = bZ[gn];);

    const FArrayBox& f0fab = f0[gn];
    const FArrayBox& f1fab = f1[gn];
    const FArrayBox& f2fab = f2[gn];
    const FArrayBox& f3fab = f3[gn];
#if (BL_SPACEDIM > 2)
    const FArrayBox& f4fab = f4[gn];
    const FArrayBox& f5fab = f5[gn];
#endif

    Real _alpha = get_alpha();
    Real _beta = get_beta();

#if (BL_SPACEDIM == 2)
    FORT_GSRB(solnfab.dataPtr(solnComp), ARLIM(solnfab.loVect()),ARLIM(solnfab.hiVect()),
              rhsfab.dataPtr(rhsComp), ARLIM(rhsfab.loVect()), ARLIM(rhsfab.hiVect()),
              &_alpha, &_beta,
              afab.dataPtr(alphaComp), ARLIM(afab.loVect()),    ARLIM(afab.hiVect()),
              bxfab.dataPtr(betaComp), ARLIM(bxfab.loVect()),   ARLIM(bxfab.hiVect()),
              byfab.dataPtr(betaComp), ARLIM(byfab.loVect()),   ARLIM(byfab.hiVect()),
              f0fab.dataPtr(), ARLIM(f0fab.loVect()),   ARLIM(f0fab.hiVect()),
              m0.dataPtr(), ARLIM(m0.loVect()),   ARLIM(m0.hiVect()),
              f1fab.dataPtr(), ARLIM(f1fab.loVect()),   ARLIM(f1fab.hiVect()),
              m1.dataPtr(), ARLIM(m1.loVect()),   ARLIM(m1.hiVect()),
              f2fab.dataPtr(), ARLIM(f2fab.loVect()),   ARLIM(f2fab.hiVect()),
              m2.dataPtr(), ARLIM(m2.loVect()),   ARLIM(m2.hiVect()),
              f3fab.dataPtr(), ARLIM(f3fab.loVect()),   ARLIM(f3fab.hiVect()),
              m3.dataPtr(), ARLIM(m3.loVect()),   ARLIM(m3.hiVect()),
              tbx.loVect(), tbx.hiVect(), vbx.loVect(), vbx.hiVect(),
              &nc, h[level], &redBlackFlag);
#endif

#if (BL_SPACEDIM == 3)
    FORT_GSRB(solnfab.dataPtr(solnComp), ARLIM(solnfab.loVect()),ARLIM(solnfab.hiVect()),
              rhsfab.dataPtr(rhsComp), ARLIM(rhsfab.loVect()), ARLIM(rhsfab.hiVect()),
              &_alpha, &_beta,
              afab.dataPtr(alphaComp), ARLIM(afab.loVect()), ARLIM(afab.hiVect()),
              bxfab.dataPtr(betaComp), ARLIM(bxfab.loVect()), ARLIM(bxfab.hiVect()),
              byfab.dataPtr(betaComp), ARLIM(byfab.loVect()), ARLIM(byfab.hiVect()),
              bzfab.dataPtr(betaComp), ARLIM(bzfab.loVect()), ARLIM(bzfab.hiVect()),
              f0fab.dataPtr(), ARLIM(f0fab.loVect()), ARLIM(f0fab.hiVect()),
              m0.dataPtr(), ARLIM(m0.loVect()), ARLIM(m0.hiVect()),
              f1fab.dataPtr(), ARLIM(f1fab.loVect()), ARLIM(f1fab.hiVect()),
              m1.dataPtr(), ARLIM(m1.loVect()), ARLIM(m1.hiVect()),
              f2fab.dataPtr(), ARLIM(f2fab.loVect()), ARLIM(f2fab.hiVect()),
              m2.dataPtr(), ARLIM(m2.loVect()), ARLIM(m2.hiVect()),
              f3fab.dataPtr(), ARLIM(f3fab.loVect()), ARLIM(f3fab.hiVect()),
              m3.dataPtr(), ARLIM(m3.loVect()), ARLIM(m3.hiVect()),
              f4fab.dataPtr(), ARLIM(f4fab.loVect()), ARLIM(f4fab.hiVect()),
              m4.dataPtr(), ARLIM(m4.loVect()), ARLIM(m4.hiVect()),
              f5fab.dataPtr(), ARLIM(f5fab.loVect()), ARLIM(f5fab.hiVect()),
              m5.dataPtr(), ARLIM(m5.loVect()), ARLIM(m5.hiVect()),
              tbx.loVect(), tbx.hiVect(), vbx.loVect(), vbx.hiVect(),
              &nc, h[level], &redBlackFlag);
#endif
  }
}

void
ABecHelper::Fsmooth_jacobi (MultiFab&       solnL,
                            const MultiFab& rhsL,
                            int             level)
{
  int alphaComp = default_alphaComp;
  int betaComp = default_betaComp;
  int bndryComp = default_bndryComp;
  int rhsComp = 0;
  int solnComp = 0;

  OrientationIter oitr;

  const FabSet& f0 = (*undrrelxr[level])[oitr()]; oitr++;
  const FabSet& f1 = (*undrrelxr[level])[oitr()]; oitr++;
  const FabSet& f2 = (*undrrelxr[level])[oitr()]; oitr++;
  const FabSet& f3 = (*undrrelxr[level])[oitr()]; oitr++;
#if (BL_SPACEDIM > 2)
  const FabSet& f4 = (*undrrelxr[level])[oitr()]; oitr++;
  const FabSet& f5 = (*undrrelxr[level])[oitr()]; oitr++;
#endif    
  const MultiFab& a = aCoefficients(level);

  D_TERM(const MultiFab& bX = bCoefficients(0,level);,
         const MultiFab& bY = bCoefficients(1,level);,
         const MultiFab& bZ = bCoefficients(2,level););

  //const int nc = solnL.nComp(); // FIXME: This LinOp only really supports single-component
  const int nc = 1;

  for (MFIter solnLmfi(solnL); solnLmfi.isValid(); ++solnLmfi)
  {
    oitr.rewind();

    const int gn = solnLmfi.index();

    const Mask& m0 = maskvals[level][oitr()][solnLmfi]; oitr++;
    const Mask& m1 = maskvals[level][oitr()][solnLmfi]; oitr++;
    const Mask& m2 = maskvals[level][oitr()][solnLmfi]; oitr++;
    const Mask& m3 = maskvals[level][oitr()][solnLmfi]; oitr++;
#if (BL_SPACEDIM > 2)
    const Mask& m4 = maskvals[level][oitr()][solnLmfi]; oitr++;
    const Mask& m5 = maskvals[level][oitr()][solnLmfi]; oitr++;
#endif
    const Box&       vbx     = solnLmfi.validbox();
    FArrayBox&       solnfab = solnL[gn];
    const FArrayBox& rhsfab  = rhsL[gn];
    const FArrayBox& afab    = a[gn];

    D_TERM(const FArrayBox& bxfab = bX[gn];,
           const FArrayBox& byfab = bY[gn];,
           const FArrayBox& bzfab = bZ[gn];);

    const FArrayBox& f0fab = f0[gn];
    const FArrayBox& f1fab = f1[gn];
    const FArrayBox& f2fab = f2[gn];
    const FArrayBox& f3fab = f3[gn];
#if (BL_SPACEDIM > 2)
    const FArrayBox& f4fab = f4[gn];
    const FArrayBox& f5fab = f5[gn];
#endif

    Real _alpha = get_alpha();
    Real _beta = get_beta();

#if (BL_SPACEDIM == 2)
    FORT_JACOBI(solnfab.dataPtr(solnComp), ARLIM(solnfab.loVect()),ARLIM(solnfab.hiVect()),
                rhsfab.dataPtr(rhsComp), ARLIM(rhsfab.loVect()), ARLIM(rhsfab.hiVect()),
                &_alpha, &_beta,
                afab.dataPtr(alphaComp), ARLIM(afab.loVect()),    ARLIM(afab.hiVect()),
                bxfab.dataPtr(betaComp), ARLIM(bxfab.loVect()),   ARLIM(bxfab.hiVect()),
                byfab.dataPtr(betaComp), ARLIM(byfab.loVect()),   ARLIM(byfab.hiVect()),
                f0fab.dataPtr(), ARLIM(f0fab.loVect()),   ARLIM(f0fab.hiVect()),
                m0.dataPtr(), ARLIM(m0.loVect()),   ARLIM(m0.hiVect()),
                f1fab.dataPtr(), ARLIM(f1fab.loVect()),   ARLIM(f1fab.hiVect()),
                m1.dataPtr(), ARLIM(m1.loVect()),   ARLIM(m1.hiVect()),
                f2fab.dataPtr(), ARLIM(f2fab.loVect()),   ARLIM(f2fab.hiVect()),
                m2.dataPtr(), ARLIM(m2.loVect()),   ARLIM(m2.hiVect()),
                f3fab.dataPtr(), ARLIM(f3fab.loVect()),   ARLIM(f3fab.hiVect()),
                m3.dataPtr(), ARLIM(m3.loVect()),   ARLIM(m3.hiVect()),
                vbx.loVect(), vbx.hiVect(),
                &nc, h[level]);
#endif

#if (BL_SPACEDIM == 3)
    FORT_JACOBI(solnfab.dataPtr(), ARLIM(solnfab.loVect()),ARLIM(solnfab.hiVect()),
                rhsfab.dataPtr(), ARLIM(rhsfab.loVect()), ARLIM(rhsfab.hiVect()),
                &_alpha, &_beta,
                afab.dataPtr(alphaComp), ARLIM(afab.loVect()), ARLIM(afab.hiVect()),
                bxfab.dataPtr(betaComp), ARLIM(bxfab.loVect()), ARLIM(bxfab.hiVect()),
                byfab.dataPtr(betaComp), ARLIM(byfab.loVect()), ARLIM(byfab.hiVect()),
                bzfab.dataPtr(betaComp), ARLIM(bzfab.loVect()), ARLIM(bzfab.hiVect()),
                f0fab.dataPtr(), ARLIM(f0fab.loVect()), ARLIM(f0fab.hiVect()),
                m0.dataPtr(), ARLIM(m0.loVect()), ARLIM(m0.hiVect()),
                f1fab.dataPtr(), ARLIM(f1fab.loVect()), ARLIM(f1fab.hiVect()),
                m1.dataPtr(), ARLIM(m1.loVect()), ARLIM(m1.hiVect()),
                f2fab.dataPtr(), ARLIM(f2fab.loVect()), ARLIM(f2fab.hiVect()),
                m2.dataPtr(), ARLIM(m2.loVect()), ARLIM(m2.hiVect()),
                f3fab.dataPtr(), ARLIM(f3fab.loVect()), ARLIM(f3fab.hiVect()),
                m3.dataPtr(), ARLIM(m3.loVect()), ARLIM(m3.hiVect()),
                f4fab.dataPtr(), ARLIM(f4fab.loVect()), ARLIM(f4fab.hiVect()),
                m4.dataPtr(), ARLIM(m4.loVect()), ARLIM(m4.hiVect()),
                f5fab.dataPtr(), ARLIM(f5fab.loVect()), ARLIM(f5fab.hiVect()),
                m5.dataPtr(), ARLIM(m5.loVect()), ARLIM(m5.hiVect()),
                vbx.loVect(), vbx.hiVect(),
                &nc, h[level]);
#endif
  }
}

