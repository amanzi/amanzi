#include <winstd.H>
#include <TensorOp.H>
#include <TensorOp_F.H>
#include <MCLO_F.H>
#include <VisMF.H>
#include <Utility.H>

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
  beta(beta_def),
  default_bndryComp(0), default_alphaComp(0), default_betaComp(0), default_beta1Comp(0)
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
  beta(beta_def),
  default_bndryComp(0), default_alphaComp(0), default_betaComp(0), default_beta1Comp(0)
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
  undrrelxrt.resize(1);
  undrrelxrt[level] = new BndryRegister(gbox[level], 1, 0, 0, nComp()); // FIXME: For 3D probably need n*(1+3)...
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

  MCLinOp::clearToLevel(level);
  undrrelxrt.resize(level+1);
  for (int i = level+1; i < numLevels(); ++i)
  {
    delete acoefs[i]; acoefs[i] = 0;
    for (int j = 0; j < BL_SPACEDIM; ++j)
    {
      delete bcoefs[i][j]; bcoefs[i][j] = 0;
      delete b1coefs[i][j]; b1coefs[i][j] = 0;
    }
  }
}

void
TensorOp::prepareForLevel (int level)
{
  MCLinOp::prepareForLevel(level);
  if (level == 0)
    return;

  TensorOp::prepareForLevel(level-1);
  //
  // Add the BndryRegister of relax values to the new coarser level.
  //
  if (undrrelxrt.size() <= level) {
    undrrelxrt.resize(level+1);
    undrrelxrt[level] = new BndryRegister(gbox[level], 1, 0, 0, numcomp);
  }
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
      acoefs.resize(level+1, 0);
      if (acoefs[level] != 0) {
        delete acoefs[level];
      }
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
    if (bcoefs.size() > level) {
      for (int lev=level; lev<bcoefs.size(); ++lev) {
        for (int i = 0; i < BL_SPACEDIM; ++i) {
          delete bcoefs[lev][i];
        }
      }
    }
    bcoefs.resize(level+1);
    
    for (int i = 0; i < BL_SPACEDIM; ++i) {
      bcoefs[level][i] = new MultiFab;
    }

    for (int i = 0; i < BL_SPACEDIM; ++i)
      makeCoefficients(*bcoefs[level][i], *bcoefs[level-1][i], level);

    b_valid.resize(level+1);
    b_valid[level] = true;
  }

  if (level >= b1_valid.size() || b1_valid[level] == false)
  {
    if (b1coefs.size() > level) {
      for (int lev=level; lev<b1coefs.size(); ++lev) {
        for (int i = 0; i < BL_SPACEDIM; ++i) {
          delete b1coefs[lev][i];
        }
      }
    }
    b1coefs.resize(level+1);
    
    for (int i = 0; i < BL_SPACEDIM; ++i) {
      b1coefs[level][i] = new MultiFab;
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

#ifndef NDEBUG
  if (BL_SPACEDIM == 3)
    BL_ASSERT(geomarray[level].IsCartesian());
#endif

  if (acoefs.size() > 0) {
    for (int lev=0; lev<acoefs.size(); ++lev) {
      delete acoefs[lev];
    }
  }
  acoefs.resize(1,0);
  acoefs[level] = new MultiFab(_ba, nComp(), nGrow);
  acoefs[level]->setVal(a_def);
  a_valid.resize(1);
  a_valid[level] = true;


  if (bcoefs.size() > 0) {
    for (int lev=0; lev<bcoefs.size(); ++lev) {
      for (int i = 0; i < BL_SPACEDIM; ++i) {
        delete bcoefs[lev][i];
      }
    }
  }
  if (b1coefs.size() > 0) {
    for (int lev=0; lev<b1coefs.size(); ++lev) {
      for (int i = 0; i < BL_SPACEDIM; ++i) {
        delete b1coefs[lev][i];
      }
    }
  }
  bcoefs.resize(1);
  b1coefs.resize(1);

  for (int i = 0; i < BL_SPACEDIM; ++i)
  {
    BoxArray edge_boxes(_ba);
    edge_boxes.surroundingNodes(i);
    bcoefs[level][i] = new MultiFab(edge_boxes, nComp(), nGrow);
    bcoefs[level][i]->setVal(b_def);

    b1coefs[level][i] = new MultiFab(edge_boxes, nComp(), nGrow);
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

void
TensorOp::apply (MultiFab& out,
                 MultiFab& in,
                 int       level,
                 MCBC_Mode bc_mode)
{
  bool local = false;
  int src_comp = 0;
  int dst_comp = 0;
  int num_comp = 1;
  int bndry_comp = default_bndryComp;
  apply(out,in,level,bc_mode,local,src_comp,dst_comp,num_comp,bndry_comp);
}

void
TensorOp::apply (MultiFab&      out,
                 MultiFab&      in,
                 int            level,
                 MCBC_Mode      bc_mode,
                 bool           local,
                 int            src_comp,
                 int            dst_comp,
                 int            num_comp,
                 int            bndry_comp)
{
  int bndryComp = (bndry_comp < 0 ? default_bndryComp : bndry_comp);
  applyBC (in,src_comp,num_comp,level,bc_mode,local,bndryComp);
  Fapply(out,dst_comp,in,src_comp,num_comp,level);
}

void
TensorOp::smooth (MultiFab&       solnL,
                  const MultiFab& rhsL,
                  int             level,
                  MCBC_Mode       bc_mode)
{
  for (int phaseflag = 0; phaseflag < numphase; phaseflag++) {
    applyBC(solnL, level, bc_mode);
    Fsmooth(solnL, rhsL, level, phaseflag);
  }
}

//
// Fills level boundary cells using BC_mode flag, int. BC data if reqd.
//
void
TensorOp::applyBC (MultiFab& inout,
                   int       level,
                   MCBC_Mode bc_mode)
{
  int sComp = 0;
  int nComp = 1;
  bool local = false;
  int bndryComp = default_bndryComp;
  applyBC (inout,sComp,nComp,level,bc_mode,local,bndryComp);
}

void
TensorOp::applyBC (MultiFab&      inout,
                   int            src_comp,
                   int            num_comp,
                   int            level,
                   MCBC_Mode      bc_mode,
                   bool           local,
                   int            bndry_comp)
{
  int bndryComp = (bndry_comp < 0 ? default_bndryComp : bndry_comp);
  //
  // No coarsened boundary values, cannot apply inhomog at lev>0.
  //
  BL_ASSERT(!(level>0 && bc_mode == MCInhomogeneous_BC));
    
  int flagden = 1;	// fill in the bndry data and undrrelxr
  int flagbc  = 1;	// with values
  if (bc_mode == MCHomogeneous_BC)
    flagbc = 0; // nodata if homog
  BL_ASSERT(num_comp==1);

  const bool cross = false;
  const bool do_corners = true;

  inout.setBndry(1.e30,src_comp,num_comp);
  inout.FillBoundary(src_comp,num_comp,local,cross);
  prepareForLevel(level);
  geomarray[level].FillPeriodicBoundary(inout,src_comp,num_comp,do_corners,local);

  //
  // Fill boundary cells.
  //
  const int N = inout.IndexArray().size();

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i < N; i++)
  {
    const int gn = inout.IndexArray()[i];

    BL_ASSERT(gbox[level][gn] == inout.box(gn));

    const BndryData::RealTuple&      bdl = bgb.bndryLocs(gn);
    const Array< Array<BoundCond> >& bdc = bgb.bndryConds(gn);
    const MaskTuple&                 msk = maskvals[level][gn];

    for (OrientationIter oitr; oitr; ++oitr)
    {
      const Orientation face = oitr();
      FabSet& f  = (*undrrelxr[level])[face];
      FabSet& ft = (*undrrelxrt[level])[face];
      FabSet& td = (*tangderiv[level])[face];
      int cdr(face);
      const FabSet& fs = bgb.bndryValues(face);
      Real bcl = bdl[face];
      const Array<BoundCond>& bc = bdc[face];
      const int *bct = (const int*) bc.dataPtr();
      const FArrayBox& fsfab = fs[gn];
      const Real* bcvalptr = fsfab.dataPtr(bndryComp);
      //
      // Way external derivs stored.
      //
      const Real* exttdptr = fsfab.dataPtr(numcomp); 
      const int* fslo      = fsfab.loVect();
      const int* fshi      = fsfab.hiVect();
      FArrayBox& inoutfab  = inout[gn];
      FArrayBox& denfab    = f[gn];
      FArrayBox& dentfab   = ft[gn];
      FArrayBox& tdfab     = td[gn];
#if BL_SPACEDIM==2
      int cdir = face.coordDir(), perpdir = -1;
      if (cdir == 0)
        perpdir = 1;
      else if (cdir == 1)
        perpdir = 0;
      else
        BoxLib::Abort("TensorOp::applyBC(): bad logic");

      const Mask& m    = *msk[face];
      const Mask& mphi = *msk[Orientation(perpdir,Orientation::high)];
      const Mask& mplo = *msk[Orientation(perpdir,Orientation::low)];
      FORT_TOAPPLYBC( &flagden, &flagbc, &maxorder,
        inoutfab.dataPtr(src_comp),  ARLIM(inoutfab.loVect()), ARLIM(inoutfab.hiVect()), &cdr, bct, &bcl,
        bcvalptr,                    ARLIM(fslo),              ARLIM(fshi),
        m.dataPtr(),                 ARLIM(m.loVect()),        ARLIM(m.hiVect()),
        mphi.dataPtr(),              ARLIM(mphi.loVect()),     ARLIM(mphi.hiVect()),
        mplo.dataPtr(),              ARLIM(mplo.loVect()),     ARLIM(mplo.hiVect()),
        denfab.dataPtr(),            ARLIM(denfab.loVect()),   ARLIM(denfab.hiVect()),
        dentfab.dataPtr(),           ARLIM(dentfab.loVect()),  ARLIM(dentfab.hiVect()),
        exttdptr,                    ARLIM(fslo),              ARLIM(fshi),
        tdfab.dataPtr(bndryComp),    ARLIM(tdfab.loVect()),    ARLIM(tdfab.hiVect()),
        inout.box(gn).loVect(), inout.box(gn).hiVect(),
        &num_comp, h[level]);
#elif BL_SPACEDIM==3
      const Mask& mn = *msk[Orientation(1,Orientation::high)];
      const Mask& me = *msk[Orientation(0,Orientation::high)];
      const Mask& mw = *msk[Orientation(0,Orientation::low)];
      const Mask& ms = *msk[Orientation(1,Orientation::low)];
      const Mask& mt = *msk[Orientation(2,Orientation::high)];
      const Mask& mb = *msk[Orientation(2,Orientation::low)];
      FORT_TOAPPLYBC( &flagden, &flagbc, &maxorder,
        inoutfab.dataPtr(src_comp),  ARLIM(inoutfab.loVect()), ARLIM(inoutfab.hiVect()), &cdr, bct, &bcl,
        bcvalptr,                    ARLIM(fslo),              ARLIM(fshi),
        mn.dataPtr(),                ARLIM(mn.loVect()),       ARLIM(mn.hiVect()),
        me.dataPtr(),                ARLIM(me.loVect()),       ARLIM(me.hiVect()),
        mw.dataPtr(),                ARLIM(mw.loVect()),       ARLIM(mw.hiVect()),
        ms.dataPtr(),                ARLIM(ms.loVect()),       ARLIM(ms.hiVect()),
        mt.dataPtr(),                ARLIM(mt.loVect()),       ARLIM(mt.hiVect()),
        mb.dataPtr(),                ARLIM(mb.loVect()),       ARLIM(mb.hiVect()),
        denfab.dataPtr(),            ARLIM(denfab.loVect()),   ARLIM(denfab.hiVect()),
        dentfab.dataPtr(),           ARLIM(dentfab.loVect()),  ARLIM(dentfab.hiVect()),
        exttdptr,                    ARLIM(fslo),              ARLIM(fshi),
        tdfab.dataPtr(bndryComp),    ARLIM(tdfab.loVect()),    ARLIM(tdfab.hiVect()),
        inout.box(gn).loVect(), inout.box(gn).hiVect(),
        &num_comp, h[level]);
#endif
    }
  }
#if 1
  // Clean up corners:
  // The problem here is that APPLYBC fills only grow cells normal to the boundary.
  // As a result, any corner cell on the boundary (either coarse-fine or fine-fine)
  // is not filled.  For coarse-fine, the operator adjusts itself, sliding away from
  // the box edge to avoid referencing that corner point.  On the physical boundary
  // though, the corner point is needed.  Particularly if a fine-fine boundary intersects
  // the physical boundary, since we want the stencil to be independent of the box
  // blocking.  FillBoundary operations wont fix the problem because the "good"
  // data we need is living in the grow region of adjacent fabs.  So, here we play
  // the usual games to treat the newly filled grow cells as "valid" data.

  // Note that we only need to do something where the grids touch the physical boundary.

  const Geometry& geomlev = geomarray[level];
  const BoxArray& grids = inout.boxArray();
  const Box& domain = geomlev.Domain();
  int nGrow = 1;

  // Lets do a quick check to see if we need to do anything at all here
  BoxArray BIGba = BoxArray(grids).grow(nGrow);

  if (! (domain.contains(BIGba.minimalBox())) ) {

    BoxArray boundary_pieces;
    Array<int> proc_idxs;
    Array<Array<int> > old_to_new(grids.size());
    const DistributionMapping& dmap=inout.DistributionMap();

    for (int d=0; d<BL_SPACEDIM; ++d) {
      if (! (geomlev.isPeriodic(d)) ) {

        BoxArray gba = BoxArray(grids).grow(d,nGrow);
        for (int i=0; i<gba.size(); ++i) {
          BoxArray new_pieces = BoxLib::boxComplement(gba[i],domain);
          int size_new = new_pieces.size();
          if (size_new>0) {
            int size_old = boundary_pieces.size();
            boundary_pieces.resize(size_old+size_new);
            proc_idxs.resize(boundary_pieces.size());
            for (int j=0; j<size_new; ++j) {
              boundary_pieces.set(size_old+j,new_pieces[j]);
              proc_idxs[size_old+j] = dmap[i];
              old_to_new[i].push_back(size_old+j);
            }
          }
        }
      }
    }

    proc_idxs.push_back(ParallelDescriptor::MyProc());

    MultiFab boundary_data(boundary_pieces,num_comp,nGrow,
                           DistributionMapping(proc_idxs));

    for (MFIter mfi(inout); mfi.isValid(); ++mfi) {
      const FArrayBox& src_fab = inout[mfi];
      for (int j=0; j<old_to_new[mfi.index()].size(); ++j) {
        int new_box_idx = old_to_new[mfi.index()][j];
        boundary_data[new_box_idx].copy(src_fab,src_comp,0,num_comp);
      }
    }

    boundary_data.FillBoundary();

    // Use a hacked Geometry object to handle the periodic intersections for us.
    // Here, the "domain" is the plane of cells on non-periodic boundary faces.
    // and there may be cells over the periodic boundary in the remaining directions.
    // We do a Geometry::PFB on each non-periodic face to sync these up.
    if (geomlev.isAnyPeriodic()) {
      Array<int> is_per(BL_SPACEDIM,0);
      for (int d=0; d<BL_SPACEDIM; ++d) {
        is_per[d] = geomlev.isPeriodic(d);
      }
      for (int d=0; d<BL_SPACEDIM; ++d) {
        if (! is_per[d]) {
          Box tmpLo = BoxLib::adjCellLo(geomlev.Domain(),d,1);
          Geometry tmpGeomLo(tmpLo,&(geomlev.ProbDomain()),(int)geomlev.Coord(),is_per.dataPtr());
          tmpGeomLo.FillPeriodicBoundary(boundary_data);

          Box tmpHi = BoxLib::adjCellHi(geomlev.Domain(),d,1);
          Geometry tmpGeomHi(tmpHi,&(geomlev.ProbDomain()),(int)geomlev.Coord(),is_per.dataPtr());
          tmpGeomHi.FillPeriodicBoundary(boundary_data);
        }
      }
    }

    for (MFIter mfi(inout); mfi.isValid(); ++mfi) {
      int idx = mfi.index();
      FArrayBox& dst_fab = inout[mfi];
      for (int j=0; j<old_to_new[idx].size(); ++j) {
        int new_box_idx = old_to_new[mfi.index()][j];
        const FArrayBox& src_fab = boundary_data[new_box_idx];
        const Box& src_box = src_fab.box();

        BoxArray pieces_outside_domain = BoxLib::boxComplement(src_box,domain);
        for (int k=0; k<pieces_outside_domain.size(); ++k) {
          const Box& outside = pieces_outside_domain[k] & dst_fab.box();
          if (outside.ok()) {
            dst_fab.copy(src_fab,outside,0,outside,src_comp,num_comp);
          }
        }
      }
    }
  }
#endif
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
  int alphaComp = default_alphaComp;
  int betaComp = default_betaComp;
  int beta1Comp = default_beta1Comp;
  int bndryComp = default_bndryComp;
  int rhsComp = 0;
  int solnComp = 0;

  OrientationIter oitr;

  const FabSet& fw  = (*undrrelxr[level])[oitr()]; 
  const FabSet& ftw = (*undrrelxrt[level])[oitr()]; 
  const FabSet& tdw = (*tangderiv[level])[oitr()];
  oitr++;
  const FabSet& fs  = (*undrrelxr[level])[oitr()]; 
  const FabSet& fts = (*undrrelxrt[level])[oitr()]; 
  const FabSet& tds = (*tangderiv[level])[oitr()];
  oitr++;
#if BL_SPACEDIM>2
  const FabSet& fb  = (*undrrelxr[level])[oitr()]; 
  const FabSet& ftb = (*undrrelxrt[level])[oitr()]; 
  const FabSet& tdb = (*tangderiv[level])[oitr()];
  oitr++;
#endif
  const FabSet& fe  = (*undrrelxr[level])[oitr()]; 
  const FabSet& fte = (*undrrelxrt[level])[oitr()]; 
  const FabSet& tde = (*tangderiv[level])[oitr()];
  oitr++;
  const FabSet& fn  = (*undrrelxr[level])[oitr()]; 
  const FabSet& ftn = (*undrrelxrt[level])[oitr()]; 
  const FabSet& tdn = (*tangderiv[level])[oitr()];
  oitr++;
#if BL_SPACEDIM>2
  const FabSet& ft  = (*undrrelxr[level])[oitr()]; 
  const FabSet& ftt = (*undrrelxrt[level])[oitr()]; 
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

  for (MFIter mfi(solnL); mfi.isValid(); ++mfi)
  {
    oitr.rewind();

    const int gn = mfi.index();

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

    D_TERM(const FArrayBox& fwfab  = fw[gn];,
	   const FArrayBox& fsfab  = fs[gn];,
	   const FArrayBox& fbfab  = fb[gn];);

    D_TERM(const FArrayBox& ftwfab = ftw[gn];,
	   const FArrayBox& ftsfab = fts[gn];,
	   const FArrayBox& ftbfab = ftb[gn];);

    D_TERM(const FArrayBox& tdwfab  = tdw[gn];,
	   const FArrayBox& tdsfab  = tds[gn];,
	   const FArrayBox& tdbfab  = tdb[gn];);

    D_TERM(const FArrayBox& fefab  = fe[gn];,
	   const FArrayBox& fnfab  = fn[gn];,
	   const FArrayBox& ftfab  = ft[gn];);

    D_TERM(const FArrayBox& ftefab = fte[gn];,
	   const FArrayBox& ftnfab = ftn[gn];,
	   const FArrayBox& fttfab = ftt[gn];);

    D_TERM(const FArrayBox& tdefab  = tde[gn];,
	   const FArrayBox& tdnfab  = tdn[gn];,
	   const FArrayBox& tdtfab  = tdt[gn];);

    D_TERM(const FArrayBox& bxfab = bX[gn];,
           const FArrayBox& byfab = bY[gn];,
           const FArrayBox& bzfab = bZ[gn];);

    D_TERM(const FArrayBox& b1xfab = b1X[gn];,
           const FArrayBox& b1yfab = b1Y[gn];,
           const FArrayBox& b1zfab = b1Z[gn];);

#if BL_SPACEDIM==2
    FORT_TOGSRB(
      solfab.dataPtr(solnComp),  ARLIM(solfab.loVect()), ARLIM(solfab.hiVect()),
      rhsfab.dataPtr(rhsComp),   ARLIM(rhsfab.loVect()), ARLIM(rhsfab.hiVect()),
      &alpha, &beta,
      afab.dataPtr(alphaComp),   ARLIM(afab.loVect()),   ARLIM(afab.hiVect()),
      bxfab.dataPtr(betaComp),   ARLIM(bxfab.loVect()),  ARLIM(bxfab.hiVect()),
      b1xfab.dataPtr(beta1Comp), ARLIM(b1xfab.loVect()), ARLIM(b1xfab.hiVect()),
      byfab.dataPtr(betaComp),   ARLIM(byfab.loVect()),  ARLIM(byfab.hiVect()),
      b1yfab.dataPtr(beta1Comp), ARLIM(b1yfab.loVect()), ARLIM(b1yfab.hiVect()),
      mn.dataPtr(),              ARLIM(mn.loVect()),     ARLIM(mn.hiVect()),
      fnfab.dataPtr(),           ARLIM(fnfab.loVect()),  ARLIM(fnfab.hiVect()),
      ftnfab.dataPtr(),          ARLIM(ftnfab.loVect()), ARLIM(ftnfab.hiVect()),
      me.dataPtr(),              ARLIM(me.loVect()),     ARLIM(me.hiVect()),
      fefab.dataPtr(),           ARLIM(fefab.loVect()),  ARLIM(fefab.hiVect()),
      ftefab.dataPtr(),          ARLIM(ftefab.loVect()), ARLIM(ftefab.hiVect()),
      mw.dataPtr(),              ARLIM(mw.loVect()),     ARLIM(mw.hiVect()),
      fwfab.dataPtr(),           ARLIM(fwfab.loVect()),  ARLIM(fwfab.hiVect()),
      ftwfab.dataPtr(),          ARLIM(ftwfab.loVect()), ARLIM(ftwfab.hiVect()),
      ms.dataPtr(),              ARLIM(ms.loVect()),     ARLIM(ms.hiVect()),
      fsfab.dataPtr(),           ARLIM(fsfab.loVect()),  ARLIM(fsfab.hiVect()),
      ftsfab.dataPtr(),          ARLIM(ftsfab.loVect()), ARLIM(ftsfab.hiVect()),
      tdnfab.dataPtr(bndryComp), ARLIM(tdnfab.loVect()), ARLIM(tdnfab.hiVect()),
      tdefab.dataPtr(bndryComp), ARLIM(tdefab.loVect()), ARLIM(tdefab.hiVect()),
      tdwfab.dataPtr(bndryComp), ARLIM(tdwfab.loVect()), ARLIM(tdwfab.hiVect()),
      tdsfab.dataPtr(bndryComp), ARLIM(tdsfab.loVect()), ARLIM(tdsfab.hiVect()),
      mfi.validbox().loVect(), mfi.validbox().hiVect(),
      h[level], phaseflag);
#else
    FORT_TOGSRB(
      solfab.dataPtr(solnComp),  ARLIM(solfab.loVect()), ARLIM(solfab.hiVect()),
      rhsfab.dataPtr(rhsComp),   ARLIM(rhsfab.loVect()), ARLIM(rhsfab.hiVect()),
      &alpha, &beta,
      afab.dataPtr(alphaComp),   ARLIM(afab.loVect()),   ARLIM(afab.hiVect()),
      bxfab.dataPtr(betaComp),   ARLIM(bxfab.loVect()),  ARLIM(bxfab.hiVect()),
      b1xfab.dataPtr(beta1Comp), ARLIM(b1xfab.loVect()), ARLIM(b1xfab.hiVect()),
      byfab.dataPtr(betaComp),   ARLIM(byfab.loVect()),  ARLIM(byfab.hiVect()),
      b1yfab.dataPtr(beta1Comp), ARLIM(b1yfab.loVect()), ARLIM(b1yfab.hiVect()),
      bzfab.dataPtr(betaComp),   ARLIM(bzfab.loVect()),  ARLIM(bzfab.hiVect()),
      b1zfab.dataPtr(beta1Comp), ARLIM(b1zfab.loVect()), ARLIM(b1zfab.hiVect()),
      mn.dataPtr(),              ARLIM(mn.loVect()),     ARLIM(mn.hiVect()),
      fnfab.dataPtr(),           ARLIM(fnfab.loVect()),  ARLIM(fnfab.hiVect()),
      ftnfab.dataPtr(),          ARLIM(ftnfab.loVect()), ARLIM(ftnfab.hiVect()),
      me.dataPtr(),              ARLIM(me.loVect()),     ARLIM(me.hiVect()),
      fefab.dataPtr(),           ARLIM(fefab.loVect()),  ARLIM(fefab.hiVect()),
      ftefab.dataPtr(),          ARLIM(ftefab.loVect()), ARLIM(ftefab.hiVect()),
      mw.dataPtr(),              ARLIM(mw.loVect()),     ARLIM(mw.hiVect()),
      fwfab.dataPtr(),           ARLIM(fwfab.loVect()),  ARLIM(fwfab.hiVect()),
      ftwfab.dataPtr(),          ARLIM(ftwfab.loVect()), ARLIM(ftwfab.hiVect()),
      ms.dataPtr(),              ARLIM(ms.loVect()),     ARLIM(ms.hiVect()),
      fsfab.dataPtr(),           ARLIM(fsfab.loVect()),  ARLIM(fsfab.hiVect()),
      ftsfab.dataPtr(),          ARLIM(ftsfab.loVect()), ARLIM(ftsfab.hiVect()),
      mt.dataPtr(),              ARLIM(mt.loVect()),     ARLIM(mt.hiVect()),
      ftfab.dataPtr(),           ARLIM(ftfab.loVect()),  ARLIM(ftfab.hiVect()),
      fttfab.dataPtr(),          ARLIM(fttfab.loVect()), ARLIM(fttfab.hiVect()),
      mb.dataPtr(),              ARLIM(mb.loVect()),     ARLIM(mb.hiVect()),
      fbfab.dataPtr(),           ARLIM(fbfab.loVect()),  ARLIM(fbfab.hiVect()),
      ftbfab.dataPtr(),          ARLIM(ftbfab.loVect()), ARLIM(ftbfab.hiVect()),
      tdnfab.dataPtr(bndryComp), ARLIM(tdnfab.loVect()), ARLIM(tdnfab.hiVect()),
      tdefab.dataPtr(bndryComp), ARLIM(tdefab.loVect()), ARLIM(tdefab.hiVect()),
      tdwfab.dataPtr(bndryComp), ARLIM(tdwfab.loVect()), ARLIM(tdwfab.hiVect()),
      tdsfab.dataPtr(bndryComp), ARLIM(tdsfab.loVect()), ARLIM(tdsfab.hiVect()),
      tdtfab.dataPtr(bndryComp), ARLIM(tdtfab.loVect()), ARLIM(tdtfab.hiVect()),
      tdbfab.dataPtr(bndryComp), ARLIM(tdbfab.loVect()), ARLIM(tdbfab.hiVect()),
      mfi.validbox().loVect(), mfi.validbox().hiVect(),
      h[level], phaseflag);
#endif
  }
}


void
TensorOp::compFlux (D_DECL(MultiFab &xflux, MultiFab &yflux, MultiFab &zflux),
                    MultiFab& x, const MCBC_Mode& bc_mode,
                    int sComp, int dComp, int nComp, int bndComp)
{
  const int level   = 0;
  bool local = false;
  applyBC(x,sComp,nComp,level,bc_mode,local,bndComp);

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

#if BL_SPACEDIM==2
    FORT_TOFLUX(
      xfab.dataPtr(sComp),     ARLIM(xfab.loVect()),     ARLIM(xfab.hiVect()),
      bxfab.dataPtr(bndComp),  ARLIM(bxfab.loVect()),    ARLIM(bxfab.hiVect()),
      b1xfab.dataPtr(bndComp), ARLIM(b1xfab.loVect()),   ARLIM(b1xfab.hiVect()),
      byfab.dataPtr(bndComp),  ARLIM(byfab.loVect()),    ARLIM(byfab.hiVect()),
      b1yfab.dataPtr(bndComp), ARLIM(b1yfab.loVect()),   ARLIM(b1yfab.hiVect()),
      xfluxfab.dataPtr(dComp), ARLIM(xfluxfab.loVect()), ARLIM(xfluxfab.hiVect()),
      yfluxfab.dataPtr(dComp), ARLIM(yfluxfab.loVect()), ARLIM(yfluxfab.hiVect()),
      mn.dataPtr(),            ARLIM(mn.loVect()),       ARLIM(mn.hiVect()),
      me.dataPtr(),            ARLIM(me.loVect()),       ARLIM(me.hiVect()),
      mw.dataPtr(),            ARLIM(mw.loVect()),       ARLIM(mw.hiVect()),
      ms.dataPtr(),            ARLIM(ms.loVect()),       ARLIM(ms.hiVect()),
      tdnfab.dataPtr(bndComp), ARLIM(tdnfab.loVect()),   ARLIM(tdnfab.hiVect()),
      tdefab.dataPtr(bndComp), ARLIM(tdefab.loVect()),   ARLIM(tdefab.hiVect()),
      tdwfab.dataPtr(bndComp), ARLIM(tdwfab.loVect()),   ARLIM(tdwfab.hiVect()),
      tdsfab.dataPtr(bndComp), ARLIM(tdsfab.loVect()),   ARLIM(tdsfab.hiVect()),
      xmfi.validbox().loVect(), xmfi.validbox().hiVect(),
      h[level]);

#else
    FORT_TOFLUX(
      xfab.dataPtr(sComp),     ARLIM(xfab.loVect()),     ARLIM(xfab.hiVect()),
      bxfab.dataPtr(bndComp),  ARLIM(bxfab.loVect()),    ARLIM(bxfab.hiVect()),
      b1xfab.dataPtr(bndComp), ARLIM(b1xfab.loVect()),   ARLIM(b1xfab.hiVect()),
      byfab.dataPtr(bndComp),  ARLIM(byfab.loVect()),    ARLIM(byfab.hiVect()),
      b1yfab.dataPtr(bndComp), ARLIM(b1yfab.loVect()),   ARLIM(b1yfab.hiVect()),
      bzfab.dataPtr(bndComp),  ARLIM(bzfab.loVect()),    ARLIM(bzfab.hiVect()),
      b1zfab.dataPtr(bndComp), ARLIM(b1zfab.loVect()),   ARLIM(b1zfab.hiVect()),
      xfluxfab.dataPtr(dComp), ARLIM(xfluxfab.loVect()), ARLIM(xfluxfab.hiVect()),
      yfluxfab.dataPtr(dComp), ARLIM(yfluxfab.loVect()), ARLIM(yfluxfab.hiVect()),
      zfluxfab.dataPtr(dComp), ARLIM(zfluxfab.loVect()), ARLIM(zfluxfab.hiVect()),
      mn.dataPtr(),            ARLIM(mn.loVect()),       ARLIM(mn.hiVect()),
      me.dataPtr(),            ARLIM(me.loVect()),       ARLIM(me.hiVect()),
      mw.dataPtr(),            ARLIM(mw.loVect()),       ARLIM(mw.hiVect()),
      ms.dataPtr(),            ARLIM(ms.loVect()),       ARLIM(ms.hiVect()),
      mt.dataPtr(),            ARLIM(mt.loVect()),       ARLIM(mt.hiVect()),
      mb.dataPtr(),            ARLIM(mb.loVect()),       ARLIM(mb.hiVect()),
      tdnfab.dataPtr(bndComp), ARLIM(tdnfab.loVect()),   ARLIM(tdnfab.hiVect()),
      tdefab.dataPtr(bndComp), ARLIM(tdefab.loVect()),   ARLIM(tdefab.hiVect()),
      tdwfab.dataPtr(bndComp), ARLIM(tdwfab.loVect()),   ARLIM(tdwfab.hiVect()),
      tdsfab.dataPtr(bndComp), ARLIM(tdsfab.loVect()),   ARLIM(tdsfab.hiVect()),
      tdtfab.dataPtr(bndComp), ARLIM(tdtfab.loVect()),   ARLIM(tdtfab.hiVect()),
      tdbfab.dataPtr(bndComp), ARLIM(tdbfab.loVect()),   ARLIM(tdbfab.hiVect()),
      xmfi.validbox().loVect(), xmfi.validbox().hiVect(),
      h[level]);
#endif
  }
}

void
TensorOp::Fapply (MultiFab&       y,
                  const MultiFab& x,
                  int             level)
{
  int num_comp = 1;
  int src_comp = 0;
  int dst_comp = 0;
  Fapply(y,dst_comp,x,src_comp,num_comp,level);
}

void
TensorOp::Fapply (MultiFab&       y,
                  int             dst_comp,
                  const MultiFab& x,
                  int             src_comp,
                  int             num_comp,
                  int             level)
{
  BL_ASSERT(y.nComp()>=dst_comp+y.nComp());
  BL_ASSERT(x.nComp()>=src_comp+x.nComp());
  BL_ASSERT(nComp()==1);
  int alphaComp = default_alphaComp;
  int betaComp = default_betaComp;
  int beta1Comp = default_beta1Comp;
  int bndryComp = default_bndryComp;

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

#if BL_SPACEDIM==2
    FORT_TOAPPLY(
      xfab.dataPtr(src_comp),    ARLIM(xfab.loVect()),   ARLIM(xfab.hiVect()),
      &alpha, &beta,
      afab.dataPtr(alphaComp),   ARLIM(afab.loVect()),   ARLIM(afab.hiVect()),
      bxfab.dataPtr(betaComp),   ARLIM(bxfab.loVect()),  ARLIM(bxfab.hiVect()),
      b1xfab.dataPtr(beta1Comp), ARLIM(b1xfab.loVect()), ARLIM(b1xfab.hiVect()),
      byfab.dataPtr(betaComp),   ARLIM(byfab.loVect()),  ARLIM(byfab.hiVect()),
      b1yfab.dataPtr(beta1Comp), ARLIM(b1yfab.loVect()), ARLIM(b1yfab.hiVect()),
      yfab.dataPtr(dst_comp),    ARLIM(yfab.loVect()),   ARLIM(yfab.hiVect()),
      mn.dataPtr(),              ARLIM(mn.loVect()),     ARLIM(mn.hiVect()),
      me.dataPtr(),              ARLIM(me.loVect()),     ARLIM(me.hiVect()),
      mw.dataPtr(),              ARLIM(mw.loVect()),     ARLIM(mw.hiVect()),
      ms.dataPtr(),              ARLIM(ms.loVect()),     ARLIM(ms.hiVect()),
      tdnfab.dataPtr(bndryComp), ARLIM(tdnfab.loVect()), ARLIM(tdnfab.hiVect()),
      tdefab.dataPtr(bndryComp), ARLIM(tdefab.loVect()), ARLIM(tdefab.hiVect()),
      tdwfab.dataPtr(bndryComp), ARLIM(tdwfab.loVect()), ARLIM(tdwfab.hiVect()),
      tdsfab.dataPtr(bndryComp), ARLIM(tdsfab.loVect()), ARLIM(tdsfab.hiVect()),
      xmfi.validbox().loVect(), xmfi.validbox().hiVect(),
      h[level]);
#else
    FORT_TOAPPLY(
      xfab.dataPtr(src_comp),    ARLIM(xfab.loVect()),   ARLIM(xfab.hiVect()),
      &alpha, &beta,
      afab.dataPtr(alphaComp),   ARLIM(afab.loVect()),   ARLIM(afab.hiVect()),
      bxfab.dataPtr(betaComp),   ARLIM(bxfab.loVect()),  ARLIM(bxfab.hiVect()),
      b1xfab.dataPtr(beta1Comp), ARLIM(b1xfab.loVect()), ARLIM(b1xfab.hiVect()),
      byfab.dataPtr(betaComp),   ARLIM(byfab.loVect()),  ARLIM(byfab.hiVect()),
      b1yfab.dataPtr(beta1Comp), ARLIM(b1yfab.loVect()), ARLIM(b1yfab.hiVect()),
      bzfab.dataPtr(betaComp),   ARLIM(bzfab.loVect()),  ARLIM(bzfab.hiVect()),
      b1zfab.dataPtr(beta1Comp), ARLIM(b1zfab.loVect()), ARLIM(b1zfab.hiVect()),
      yfab.dataPtr(dst_comp),    ARLIM(yfab.loVect()),   ARLIM(yfab.hiVect()),
      mn.dataPtr(),              ARLIM(mn.loVect()),     ARLIM(mn.hiVect()),
      me.dataPtr(),              ARLIM(me.loVect()),     ARLIM(me.hiVect()),
      mw.dataPtr(),              ARLIM(mw.loVect()),     ARLIM(mw.hiVect()),
      ms.dataPtr(),              ARLIM(ms.loVect()),     ARLIM(ms.hiVect()),
      mt.dataPtr(),              ARLIM(mt.loVect()),     ARLIM(mt.hiVect()),
      mb.dataPtr(),              ARLIM(mb.loVect()),     ARLIM(mb.hiVect()),
      tdnfab.dataPtr(bndryComp), ARLIM(tdnfab.loVect()), ARLIM(tdnfab.hiVect()),
      tdefab.dataPtr(bndryComp), ARLIM(tdefab.loVect()), ARLIM(tdefab.hiVect()),
      tdwfab.dataPtr(bndryComp), ARLIM(tdwfab.loVect()), ARLIM(tdwfab.hiVect()),
      tdsfab.dataPtr(bndryComp), ARLIM(tdsfab.loVect()), ARLIM(tdsfab.hiVect()),
      tdtfab.dataPtr(bndryComp), ARLIM(tdtfab.loVect()), ARLIM(tdtfab.hiVect()),
      tdbfab.dataPtr(bndryComp), ARLIM(tdbfab.loVect()), ARLIM(tdbfab.hiVect()),
      xmfi.validbox().loVect(), xmfi.validbox().hiVect(),
      h[level]);
#endif
  }
}
