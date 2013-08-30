#include <winstd.H>
#include <TensorOp.H>
#include <TensorOp_F.H>

Real TensorOp::a_def     = 0.0;
Real TensorOp::b_def     = 1.0;
Real TensorOp::b1_def    = 0.0;
Real TensorOp::alpha_def = 1.0;
Real TensorOp::beta_def  = 1.0;

int
TensorOp::numberComponents ()
{
  return 1;
}

int
TensorOp::numberPhases ()
{
  return BL_SPACEDIM==2 ? 4 : 8;
}

TensorOp::TensorOp (const BndryData& _bd,
                    Real             _h,
                    int              _nc)
  :
  MCLinOp(_bd, _h, _nc),
  alpha(alpha_def),
  beta(beta_def)
{
  Real __h[BL_SPACEDIM];

  D_TERM(__h[0]=_h;, __h[1]=_h;, __h[2]=_h;);

  initConstruct(__h);
}

TensorOp::TensorOp (const BndryData& _bd,
                    const Real*      _h,
                    int              _nc)
  :
  MCLinOp(_bd, _h, _nc),
  alpha(alpha_def),
  beta(beta_def)
{
  initConstruct(_h);
}

void
TensorOp::initConstruct (const Real* _h)
{
  const int level       = 0;

  initCoefficients(gbox[level]);

  numphase = numberPhases();     // wyc

  undrrelxr.resize(1);
  undrrelxr[level] = new BndryRegister(gbox[level], 1, 0, 0, nComp());
  tangderiv.resize(1);
#if BL_SPACEDIM==2
  tangderiv[level] = new BndryRegister(gbox[level], 0, 1, 0, nComp());
#elif BL_SPACEDIM==3
  tangderiv[level] = new BndryRegister(gbox[level], 0, 1, 0, nComp()*(1+3));
#else
# error "BL_SPACEDIME must be 2 or 3"
#endif
}

TensorOp::~TensorOp ()
{
  clearToLevel(-1);
}

void
TensorOp::setScalars (Real _alpha,
                      Real _beta)
{
  alpha = _alpha;
  beta  = _beta;
}

void
TensorOp::clearToLevel (int level)
{
  BL_ASSERT(level >= -1);

  for (int i = level+1; i < numLevels(); ++i)
  {
    delete acoefs[i];
    for (int j = 0; j < BL_SPACEDIM; ++j)
    {
      delete bcoefs[i][j];
      delete b1coefs[i][j];
    }
  }
}

void
TensorOp::prepareForLevel (int level)
{
  MCLinOp::prepareForLevel(level);
  if (level == 0)
    return;
  prepareForLevel(level-1);
  //
  // If coefficients were marked invalid, or if not yet made, make new ones
  // (Note: makeCoefficients is a MCLinOp routine, and it allocates AND
  // fills coefficients.  A more efficient implementation would allocate
  // and fill in separate steps--we could then use the a_valid bool
  // along with the length of a_valid to separately determine whether to
  // fill or allocate the coefficient MultiFabs.
  //
  if (level >= a_valid.size() || a_valid[level] == false)
  {
    if (acoefs.size() < level+1)
    {
      acoefs.resize(level+1);
      acoefs[level] = new MultiFab;
    }
    else
    {
      delete acoefs[level];
      acoefs[level] = new MultiFab;
    }
    makeCoefficients(*acoefs[level], *acoefs[level-1], level);
    a_valid.resize(level+1);
    a_valid[level] = true;
  }
    
  if (level >= b_valid.size() || b_valid[level] == false)
  {
    if (bcoefs.size() < level+1)
    {
      bcoefs.resize(level+1);
      for (int i = 0; i < BL_SPACEDIM; ++i)
        bcoefs[level][i] = new MultiFab;
    }
    else
    {
      for (int i = 0; i < BL_SPACEDIM; ++i)
      {
        delete bcoefs[level][i];
        bcoefs[level][i] = new MultiFab;
      }
    }
    for (int i = 0; i < BL_SPACEDIM; ++i)
      makeCoefficients(*bcoefs[level][i], *bcoefs[level-1][i], level);

    b_valid.resize(level+1);
    b_valid[level] = true;
  }

  if (level >= b1_valid.size() || b1_valid[level] == false)
  {
    if (b1coefs.size() < level+1)
    {
      b1coefs.resize(level+1);
      for (int i = 0; i < BL_SPACEDIM; ++i)
        b1coefs[level][i] = new MultiFab;
    }
    else
    {
      for (int i = 0; i < BL_SPACEDIM; ++i)
      {
        delete b1coefs[level][i];
        b1coefs[level][i] = new MultiFab;
      }
    }
    for (int i = 0; i < BL_SPACEDIM; ++i)
      makeCoefficients(*b1coefs[level][i], *b1coefs[level-1][i], level);

    b1_valid.resize(level+1);
    b1_valid[level] = true;
  }
}

void
TensorOp::initCoefficients (const BoxArray &_ba)
{
  const int nGrow = 0;
  const int level = 0;

  acoefs.resize(1);
  bcoefs.resize(1);
  b1coefs.resize(1);

#ifndef NDEBUG
  if (BL_SPACEDIM == 3)
    BL_ASSERT(geomarray[level].IsCartesian());
#endif

  acoefs[level] = new MultiFab(_ba, nComp(), nGrow);
  acoefs[level]->setVal(a_def);
  a_valid.resize(1);
  a_valid[level] = true;

  for (int i = 0; i < BL_SPACEDIM; ++i)
  {
    BoxArray edge_boxes(_ba);
    edge_boxes.surroundingNodes(i);
    bcoefs[level][i] = new MultiFab(edge_boxes, nComp(), nGrow);
    b1coefs[level][i] = new MultiFab(edge_boxes, nComp(), nGrow);
    bcoefs[level][i]->setVal(b_def);
    b1coefs[level][i]->setVal(b1_def);
  }
  b_valid.resize(1);
  b1_valid.resize(1);
  b_valid[level] = true;
  b1_valid[level] = true;
}

void
TensorOp::invalidate_a_to_level (int lev)
{
  lev = (lev >= 0 ? lev : 0);
  for (int i = lev; i < numLevels(); i++)
    a_valid[i]=false;
}

void
TensorOp::invalidate_b_to_level (int lev)
{
  lev = (lev >= 0 ? lev : 0);
  for (int i = lev; i < numLevels(); i++)
    b_valid[i]=false;
}

void
TensorOp::invalidate_b1_to_level (int lev)
{
  lev = (lev >= 0 ? lev : 0);
  for (int i = lev; i < numLevels(); i++)
    b1_valid[i]=false;
}


void
TensorOp::aCoefficients (const MultiFab& _a)
{
  BL_ASSERT(_a.ok());
  BL_ASSERT(_a.boxArray() == (acoefs[0])->boxArray());
  BL_ASSERT(_a.nComp() >= nComp());
  invalidate_a_to_level(0);
  (*acoefs[0]).copy(_a,0,0,nComp());
}

void
TensorOp::bCoefficients (const MultiFab& _b,
                         int             dir)
{
  BL_ASSERT(_b.ok());
  BL_ASSERT(_b.boxArray() == (bcoefs[0][dir])->boxArray());
  BL_ASSERT(_b.nComp() >= nComp());
  invalidate_b_to_level(0);
  (*bcoefs[0][dir]).copy(_b,0,0,nComp());
}

void
TensorOp::b1Coefficients (const MultiFab& _b1,
                          int             dir)
{
  BL_ASSERT(_b1.ok());
  BL_ASSERT(_b1.boxArray() == (b1coefs[0][dir])->boxArray());
  BL_ASSERT(_b1.nComp() >= nComp());
  invalidate_b1_to_level(0);
  (*b1coefs[0][dir]).copy(_b1,0,0,nComp());
}

const MultiFab&
TensorOp::aCoefficients (int level)
{
  prepareForLevel(level);
  return *acoefs[level];
}

const MultiFab&
TensorOp::bCoefficients (int dir,
                         int level)
{
  prepareForLevel(level);
  return *bcoefs[level][dir];
}

const MultiFab&
TensorOp::b1Coefficients (int dir,
                          int level)
{
  prepareForLevel(level);
  return *b1coefs[level][dir];
}

//
// Must be defined for MultiGrid/CGSolver to work.
//
void
TensorOp::Fsmooth (MultiFab&       solnL,
                   const MultiFab& rhsL,
                   int             level,
                   int             phaseflag)
{
  OrientationIter oitr;

  const FabSet& fw  = (*undrrelxr[level])[oitr()]; 
  const FabSet& tdw = (*tangderiv[level])[oitr()];
  oitr++;
  const FabSet& fs  = (*undrrelxr[level])[oitr()]; 
  const FabSet& tds = (*tangderiv[level])[oitr()];
  oitr++;
#if BL_SPACEDIM>2
  const FabSet& fb  = (*undrrelxr[level])[oitr()]; 
  const FabSet& tdb = (*tangderiv[level])[oitr()];
  oitr++;
#endif
  const FabSet& fe  = (*undrrelxr[level])[oitr()]; 
  const FabSet& tde = (*tangderiv[level])[oitr()];
  oitr++;
  const FabSet& fn  = (*undrrelxr[level])[oitr()]; 
  const FabSet& tdn = (*tangderiv[level])[oitr()];
  oitr++;
#if BL_SPACEDIM>2
  const FabSet& ft  = (*undrrelxr[level])[oitr()]; 
  const FabSet& tdt = (*tangderiv[level])[oitr()];
  oitr++;
#endif
  const MultiFab& a  = aCoefficients(level);

  D_TERM(const MultiFab& bX = bCoefficients(0,level);,
         const MultiFab& bY = bCoefficients(1,level);,
         const MultiFab& bZ = bCoefficients(2,level););

  D_TERM(const MultiFab& b1X = b1Coefficients(0,level);,
         const MultiFab& b1Y = b1Coefficients(1,level);,
         const MultiFab& b1Z = b1Coefficients(2,level););

  int nc = solnL.nComp();
  BL_ASSERT(nc>=nComp());
  BL_ASSERT(nc==1);

  for (MFIter solnLmfi(solnL); solnLmfi.isValid(); ++solnLmfi)
  {
    oitr.rewind();

    const int gn = solnLmfi.index();

    const MCLinOp::MaskTuple& mtuple = maskvals[level][gn];

    D_TERM(const Mask& mw = *mtuple[oitr()]; oitr++;,
           const Mask& ms = *mtuple[oitr()]; oitr++;,
           const Mask& mb = *mtuple[oitr()]; oitr++;);

    D_TERM(const Mask& me = *mtuple[oitr()]; oitr++;,
           const Mask& mn = *mtuple[oitr()]; oitr++;,
           const Mask& mt = *mtuple[oitr()]; oitr++;);

    FArrayBox&       solfab = solnL[gn];
    const FArrayBox& rhsfab = rhsL[gn];
    const FArrayBox& afab   = a[gn];
    const FArrayBox& fnfab  = fn[gn];
    const FArrayBox& fefab  = fe[gn];
    const FArrayBox& fwfab  = fw[gn];
    const FArrayBox& fsfab  = fs[gn];
    const FArrayBox& tdnfab = tdn[gn];
    const FArrayBox& tdefab = tde[gn];
    const FArrayBox& tdwfab = tdw[gn];
    const FArrayBox& tdsfab = tds[gn];

#if BL_SPACEDIM>2
    const FArrayBox& ftfab  = ft[gn];
    const FArrayBox& fbfab  = fb[gn];
    const FArrayBox& tdtfab = tdt[gn];
    const FArrayBox& tdbfab = tdb[gn];
#endif

    D_TERM(const FArrayBox& bxfab = bX[gn];,
           const FArrayBox& byfab = bY[gn];,
           const FArrayBox& bzfab = bZ[gn];);

    D_TERM(const FArrayBox& b1xfab = b1X[gn];,
           const FArrayBox& b1yfab = b1Y[gn];,
           const FArrayBox& b1zfab = b1Z[gn];);

    FORT_GSRB(
      solfab.dataPtr(), 
      ARLIM(solfab.loVect()),ARLIM(solfab.hiVect()),
      rhsfab.dataPtr(),
      ARLIM(rhsfab.loVect()), ARLIM(rhsfab.hiVect()),
      &alpha, &beta,
      afab.dataPtr(),
      ARLIM(afab.loVect()),    ARLIM(afab.hiVect()),
      bxfab.dataPtr(),
      ARLIM(bxfab.loVect()),   ARLIM(bxfab.hiVect()),
      b1xfab.dataPtr(),
      ARLIM(b1xfab.loVect()),   ARLIM(b1xfab.hiVect()),
      byfab.dataPtr(),
      ARLIM(byfab.loVect()),   ARLIM(byfab.hiVect()),
      b1yfab.dataPtr(),
      ARLIM(b1yfab.loVect()),   ARLIM(b1yfab.hiVect()),
#if BL_SPACEDIM>2
      bzfab.dataPtr(),
      ARLIM(bzfab.loVect()),   ARLIM(bzfab.hiVect()),
      b1zfab.dataPtr(),
      ARLIM(b1zfab.loVect()),   ARLIM(b1zfab.hiVect()),
#endif
      mn.dataPtr(),
      ARLIM(mn.loVect()),ARLIM(mn.hiVect()),
      fnfab.dataPtr(),
      ARLIM(fnfab.loVect()),   ARLIM(fnfab.hiVect()),
      me.dataPtr(),
      ARLIM(me.loVect()),ARLIM(me.hiVect()),
      fefab.dataPtr(),
      ARLIM(fefab.loVect()),   ARLIM(fefab.hiVect()),
      mw.dataPtr(),
      ARLIM(mw.loVect()),ARLIM(mw.hiVect()),
      fwfab.dataPtr(),
      ARLIM(fwfab.loVect()),   ARLIM(fwfab.hiVect()),
      ms.dataPtr(),
      ARLIM(ms.loVect()),ARLIM(ms.hiVect()),
      fsfab.dataPtr(),
      ARLIM(fsfab.loVect()),   ARLIM(fsfab.hiVect()),
#if BL_SPACEDIM>2
      mt.dataPtr(),
      ARLIM(mt.loVect()),ARLIM(mt.hiVect()),
      ftfab.dataPtr(),
      ARLIM(ftfab.loVect()),   ARLIM(ftfab.hiVect()),
      mb.dataPtr(),
      ARLIM(mb.loVect()),ARLIM(mb.hiVect()),
      fbfab.dataPtr(),
      ARLIM(fbfab.loVect()),   ARLIM(fbfab.hiVect()),
#endif
      tdnfab.dataPtr(),
      ARLIM(tdnfab.loVect()),ARLIM(tdnfab.hiVect()),
      tdefab.dataPtr(),
      ARLIM(tdefab.loVect()),ARLIM(tdefab.hiVect()),
      tdwfab.dataPtr(),
      ARLIM(tdwfab.loVect()),ARLIM(tdwfab.hiVect()),
      tdsfab.dataPtr(),
      ARLIM(tdsfab.loVect()),ARLIM(tdsfab.hiVect()),
#if BL_SPACEDIM>2
      tdtfab.dataPtr(),
      ARLIM(tdtfab.loVect()),ARLIM(tdtfab.hiVect()),
      tdbfab.dataPtr(),
      ARLIM(tdbfab.loVect()),ARLIM(tdbfab.hiVect()),
#endif
      solnLmfi.validbox().loVect(), solnLmfi.validbox().hiVect(),
      h[level], phaseflag);
  }
}


void
TensorOp::compFlux (D_DECL(MultiFab &xflux, MultiFab &yflux, MultiFab &zflux),
                    MultiFab& x, const MCBC_Mode& bc_mode,
                    int sComp, int dComp, int nComp, int bndComp)
{
  const int level   = 0;
  applyBC(x,level,bc_mode);
  BL_ASSERT(x.nComp()>=sComp+nComp);
  BL_ASSERT(xflux.nComp()>=dComp+nComp);
  BL_ASSERT(yflux.nComp()>=dComp+nComp);
#if (BL_SPACEDIM>2)
  BL_ASSERT(zflux.nComp()>=dComp+nComp);
#endif
  BL_ASSERT(nComp==1);
    
  D_TERM(const MultiFab& bX = bCoefficients(0,level);,
         const MultiFab& bY = bCoefficients(1,level);,
         const MultiFab& bZ = bCoefficients(2,level););

  D_TERM(const MultiFab& b1X = b1Coefficients(0,level);,
         const MultiFab& b1Y = b1Coefficients(1,level);,
         const MultiFab& b1Z = b1Coefficients(2,level););

  OrientationIter oitr;

  D_TERM(const FabSet& tdw = (*tangderiv[level])[oitr()]; oitr++;,
         const FabSet& tds = (*tangderiv[level])[oitr()]; oitr++;,
         const FabSet& tdb = (*tangderiv[level])[oitr()]; oitr++;);

  D_TERM(const FabSet& tde = (*tangderiv[level])[oitr()]; oitr++;,
         const FabSet& tdn = (*tangderiv[level])[oitr()]; oitr++;,
         const FabSet& tdt = (*tangderiv[level])[oitr()]; oitr++;);

  for (MFIter xmfi(x); xmfi.isValid(); ++xmfi)
  {
    oitr.rewind();

    const int gn = xmfi.index();

    const MCLinOp::MaskTuple& mtuple = maskvals[level][gn];

    D_TERM(const Mask& mw = *mtuple[oitr()]; oitr++;,
           const Mask& ms = *mtuple[oitr()]; oitr++;,
           const Mask& mb = *mtuple[oitr()]; oitr++;);

    D_TERM(const Mask& me = *mtuple[oitr()]; oitr++;,
           const Mask& mn = *mtuple[oitr()]; oitr++;,
           const Mask& mt = *mtuple[oitr()]; oitr++;);

    FArrayBox&       xfab = x[gn];
    const FArrayBox& tdnfab = tdn[gn];
    const FArrayBox& tdefab = tde[gn];
    const FArrayBox& tdwfab = tdw[gn];
    const FArrayBox& tdsfab = tds[gn];

#if BL_SPACEDIM>2
    const FArrayBox& tdtfab = tdt[gn];
    const FArrayBox& tdbfab = tdb[gn];
#endif

    D_TERM(const FArrayBox& bxfab = bX[gn];,
           const FArrayBox& byfab = bY[gn];,
           const FArrayBox& bzfab = bZ[gn];);

    D_TERM(const FArrayBox& b1xfab = b1X[gn];,
           const FArrayBox& b1yfab = b1Y[gn];,
           const FArrayBox& b1zfab = b1Z[gn];);

    D_TERM(FArrayBox& xfluxfab = xflux[gn];,
           FArrayBox& yfluxfab = yflux[gn];,
           FArrayBox& zfluxfab = zflux[gn];);

    FORT_FLUX(
      xfab.dataPtr(sComp), 
      ARLIM(xfab.loVect()), ARLIM(xfab.hiVect()),
      bxfab.dataPtr(bndComp), 
      ARLIM(bxfab.loVect()), ARLIM(bxfab.hiVect()),
      b1xfab.dataPtr(bndComp), 
      ARLIM(b1xfab.loVect()), ARLIM(b1xfab.hiVect()),
      byfab.dataPtr(bndComp), 
      ARLIM(byfab.loVect()), ARLIM(byfab.hiVect()),
      b1yfab.dataPtr(bndComp), 
      ARLIM(b1yfab.loVect()), ARLIM(b1yfab.hiVect()),
#if BL_SPACEDIM>2
      bzfab.dataPtr(bndComp), 
      ARLIM(bzfab.loVect()), ARLIM(bzfab.hiVect()),
      b1zfab.dataPtr(bndComp), 
      ARLIM(b1zfab.loVect()), ARLIM(b1zfab.hiVect()),
#endif
      xfluxfab.dataPtr(dComp), 
      ARLIM(xfluxfab.loVect()), ARLIM(xfluxfab.hiVect()),
      yfluxfab.dataPtr(dComp), 
      ARLIM(yfluxfab.loVect()), ARLIM(yfluxfab.hiVect()),
#if BL_SPACEDIM>2
      zfluxfab.dataPtr(dComp), 
      ARLIM(zfluxfab.loVect()), ARLIM(zfluxfab.hiVect()),
#endif
      mn.dataPtr(),
      ARLIM(mn.loVect()),ARLIM(mn.hiVect()),
      me.dataPtr(),
      ARLIM(me.loVect()),ARLIM(me.hiVect()),
      mw.dataPtr(),
      ARLIM(mw.loVect()),ARLIM(mw.hiVect()),
      ms.dataPtr(),
      ARLIM(ms.loVect()),ARLIM(ms.hiVect()),
#if BL_SPACEDIM>2
      mt.dataPtr(),
      ARLIM(mt.loVect()),ARLIM(mt.hiVect()),
      mb.dataPtr(),
      ARLIM(mb.loVect()),ARLIM(mb.hiVect()),
#endif
      tdnfab.dataPtr(bndComp),
      ARLIM(tdnfab.loVect()),ARLIM(tdnfab.hiVect()),
      tdefab.dataPtr(bndComp),
      ARLIM(tdefab.loVect()),ARLIM(tdefab.hiVect()),
      tdwfab.dataPtr(bndComp),
      ARLIM(tdwfab.loVect()),ARLIM(tdwfab.hiVect()),
      tdsfab.dataPtr(bndComp),
      ARLIM(tdsfab.loVect()),ARLIM(tdsfab.hiVect()),
#if BL_SPACEDIM>2
      tdtfab.dataPtr(bndComp),
      ARLIM(tdtfab.loVect()),ARLIM(tdtfab.hiVect()),
      tdbfab.dataPtr(bndComp),
      ARLIM(tdbfab.loVect()),ARLIM(tdbfab.hiVect()),
#endif
      xmfi.validbox().loVect(), xmfi.validbox().hiVect(),
      h[level]);
  }
}

void
TensorOp::Fapply (MultiFab&       y,
                  const MultiFab& x,
                  int             level)
{
  BL_ASSERT(y.nComp()>=nComp());
  BL_ASSERT(x.nComp()>=nComp());
  BL_ASSERT(nComp()==1);

  const MultiFab& a = aCoefficients(level);

  D_TERM(const MultiFab& bX = bCoefficients(0,level);,
         const MultiFab& bY = bCoefficients(1,level);,
         const MultiFab& bZ = bCoefficients(2,level););

  D_TERM(const MultiFab& b1X = b1Coefficients(0,level);,
         const MultiFab& b1Y = b1Coefficients(1,level);,
         const MultiFab& b1Z = b1Coefficients(2,level););

  OrientationIter oitr;

  D_TERM(const FabSet& tdw = (*tangderiv[level])[oitr()]; oitr++;,
         const FabSet& tds = (*tangderiv[level])[oitr()]; oitr++;,
         const FabSet& tdb = (*tangderiv[level])[oitr()]; oitr++;);

  D_TERM(const FabSet& tde = (*tangderiv[level])[oitr()]; oitr++;,
         const FabSet& tdn = (*tangderiv[level])[oitr()]; oitr++;,
         const FabSet& tdt = (*tangderiv[level])[oitr()]; oitr++;);

  for (MFIter xmfi(x); xmfi.isValid(); ++xmfi)
  {
    oitr.rewind();

    const int gn = xmfi.index();

    const MCLinOp::MaskTuple& mtuple = maskvals[level][gn];

    D_TERM(const Mask& mw = *mtuple[oitr()]; oitr++;,
           const Mask& ms = *mtuple[oitr()]; oitr++;,
           const Mask& mb = *mtuple[oitr()]; oitr++;);

    D_TERM(const Mask& me = *mtuple[oitr()]; oitr++;,
           const Mask& mn = *mtuple[oitr()]; oitr++;,
           const Mask& mt = *mtuple[oitr()]; oitr++;);

    FArrayBox&       yfab = y[gn];
    const FArrayBox& xfab = x[gn];

    const FArrayBox& afab = a[gn];
    const FArrayBox& tdnfab = tdn[gn];
    const FArrayBox& tdefab = tde[gn];
    const FArrayBox& tdwfab = tdw[gn];
    const FArrayBox& tdsfab = tds[gn];

#if BL_SPACEDIM>2
    const FArrayBox& tdtfab = tdt[gn];
    const FArrayBox& tdbfab = tdb[gn];
#endif
    D_TERM(const FArrayBox& bxfab = bX[gn];,
           const FArrayBox& byfab = bY[gn];,
           const FArrayBox& bzfab = bZ[gn];);

    D_TERM(const FArrayBox& b1xfab = b1X[gn];,
           const FArrayBox& b1yfab = b1Y[gn];,
           const FArrayBox& b1zfab = b1Z[gn];);

    FORT_APPLY(
      xfab.dataPtr(), 
      ARLIM(xfab.loVect()), ARLIM(xfab.hiVect()),
      &alpha, &beta,
      afab.dataPtr(), 
      ARLIM(afab.loVect()), ARLIM(afab.hiVect()),
      bxfab.dataPtr(), 
      ARLIM(bxfab.loVect()), ARLIM(bxfab.hiVect()),
      b1xfab.dataPtr(), 
      ARLIM(b1xfab.loVect()), ARLIM(b1xfab.hiVect()),
      byfab.dataPtr(), 
      ARLIM(byfab.loVect()), ARLIM(byfab.hiVect()),
      b1yfab.dataPtr(), 
      ARLIM(b1yfab.loVect()), ARLIM(b1yfab.hiVect()),
#if BL_SPACEDIM>2
      bzfab.dataPtr(), 
      ARLIM(bzfab.loVect()), ARLIM(bzfab.hiVect()),
      b1zfab.dataPtr(), 
      ARLIM(b1zfab.loVect()), ARLIM(b1zfab.hiVect()),
#endif
      yfab.dataPtr(), 
      ARLIM(yfab.loVect()), ARLIM(yfab.hiVect()),
      mn.dataPtr(),
      ARLIM(mn.loVect()),ARLIM(mn.hiVect()),
      me.dataPtr(),
      ARLIM(me.loVect()),ARLIM(me.hiVect()),
      mw.dataPtr(),
      ARLIM(mw.loVect()),ARLIM(mw.hiVect()),
      ms.dataPtr(),
      ARLIM(ms.loVect()),ARLIM(ms.hiVect()),
#if BL_SPACEDIM>2
      mt.dataPtr(),
      ARLIM(mt.loVect()),ARLIM(mt.hiVect()),
      mb.dataPtr(),
      ARLIM(mb.loVect()),ARLIM(mb.hiVect()),
#endif
      tdnfab.dataPtr(),
      ARLIM(tdnfab.loVect()),ARLIM(tdnfab.hiVect()),
      tdefab.dataPtr(),
      ARLIM(tdefab.loVect()),ARLIM(tdefab.hiVect()),
      tdwfab.dataPtr(),
      ARLIM(tdwfab.loVect()),ARLIM(tdwfab.hiVect()),
      tdsfab.dataPtr(),
      ARLIM(tdsfab.loVect()),ARLIM(tdsfab.hiVect()),
#if BL_SPACEDIM>2
      tdtfab.dataPtr(),
      ARLIM(tdtfab.loVect()),ARLIM(tdtfab.hiVect()),
      tdbfab.dataPtr(),
      ARLIM(tdbfab.loVect()),ARLIM(tdbfab.hiVect()),
#endif
      xmfi.validbox().loVect(), xmfi.validbox().hiVect(),
      h[level]);
  }
}
