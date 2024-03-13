/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <MFTower.H>
#include <Layout.H>
#include <Utility.H>
#include <VisMF.H>
#include <MFTower_F.H>
#include <WritePlotfile.H>

MFTower::MFTower(const Layout&    _layout,
                 const IndexType& t,
                 int              _nComp,
                 int              _nGrow,
		 int              _nLevs)
    : layout(_layout), iType(t), nComp(_nComp), nGrow(_nGrow), nLevs(_nLevs)
{
  BL_PROFILE("MFTower::MFTower()");
  if (nLevs<0) nLevs = layout.NumLevels();
  BL_ASSERT(nLevs <= layout.NumLevels());
  define_alloc();
}

MFTower::MFTower(Layout&           _layout,
                 PArray<MultiFab>& pamf,
		 int               _nLevs)
    : layout(_layout), nComp(pamf[0].nComp()), nGrow(pamf[0].nGrow()), nLevs(_nLevs)
{
  BL_PROFILE("MFTower::MFTower1()");
  if (nLevs<0) nLevs = layout.NumLevels();
  BL_ASSERT(nLevs <= layout.NumLevels());
  define_noalloc(pamf);
}

MultiFab& MFTower::operator[](int i) {return mft[i];}
const MultiFab& MFTower::operator[](int i) const {return mft[i];}

const Layout&
MFTower::GetLayout() const
{
    return layout;
}

int
MFTower::NComp() const
{
    return nComp;
}

int
MFTower::NGrow() const
{
    return nGrow;
}

int
MFTower::NumLevels() const
{
   return nLevs;
}

IndexType MFTower::EC[BL_SPACEDIM] = {D_DECL(IndexType(BoxLib::BASISV(0)),
                                             IndexType(BoxLib::BASISV(1)),
                                             IndexType(BoxLib::BASISV(2)))};
IndexType MFTower::CC(IntVect::TheZeroVector());
IndexType MFTower::NC(IntVect::TheUnitVector());

const IndexType&
MFTower::ixType() const
{
    return iType;
}


void
MFTower::define_alloc()
{
  BL_PROFILE("MFTower::define_alloc()");
  mft.resize(nLevs,PArrayManage);
  const Array<BoxArray>& gridArray = layout.GridArray();
  for (int lev=0; lev<nLevs; ++lev)
  {
    BoxArray ba = gridArray[lev];
    if (iType!=CC) {
      if (iType == NC) {
        ba.surroundingNodes();
      }
      else {
        for (int d=0; d<BL_SPACEDIM; ++d) {
          if (iType == EC[d]) {
            ba.surroundingNodes(d);
          }
        }
      }
    }
    mft.set(lev,new MultiFab(ba,nComp,nGrow));
  }
}

void
MFTower::define_noalloc(PArray<MultiFab>& pamf)
{
  BL_PROFILE("MFTower::define_noalloc()");
    iType = pamf[0].boxArray()[0].ixType();
    mft.resize(nLevs,PArrayNoManage);
    const Array<BoxArray>& gridArray = layout.GridArray();
    for (int lev=0; lev<nLevs; ++lev)
    {
        BoxArray ba = pamf[lev].boxArray();
        if (iType!=CC) {
            ba.enclosedCells();
        }
        BL_ASSERT(ba == gridArray[lev]);
        BL_ASSERT(pamf[lev].nComp() == nComp);
        BL_ASSERT(pamf[lev].nGrow() == nGrow);
        mft.set(lev,&(pamf[lev]));
    }
}

Real
MFTower::norm(int numLevs) const // currently only max norm supported
{
  BL_PROFILE("MFTower::norm()");
    if (numLevs<0) numLevs=nLevs;
    BL_ASSERT(numLevs<=nLevs);
    FArrayBox fab;
    Real norm = 0;
    const Array<IntVect>& refRatio = layout.RefRatio();
    const Array<BoxArray>& gridArray = layout.GridArray();
    for (int lev=0; lev<numLevs; ++lev)
    {
        for (MFIter mfi(mft[lev]); mfi.isValid(); ++mfi)
        {
            Box vbox = mfi.validbox();
            fab.resize(vbox,1);
            fab.copy(mft[lev][mfi]);
            if (lev<numLevs-1) {
                BoxArray cfba = BoxArray(gridArray[lev-1]).coarsen(refRatio[lev-1]);
                std::vector< std::pair<int,Box> > isects = cfba.intersections(vbox);
                for (int i = 0; i < isects.size(); i++)
                {
                    Box ovlp = isects[i].second & vbox;
                    if (ovlp.ok()) {
                        fab.setVal(0);
                    }
                }
            }
            norm = std::max(norm, fab.norm(2,0,1));
        }
    }
    ParallelDescriptor::ReduceRealMax(norm);
    return norm;
}

void
MFTower::SetValCovered(Real value)
{
  BL_PROFILE("MFTower::SetValCovered()");
  const Array<IntVect>& refRatio = layout.RefRatio();
  const Array<BoxArray>& gridArray = layout.GridArray();
  for (int lev=0; lev<nLevs-1; ++lev) {
    for (MFIter mfi(mft[lev]); mfi.isValid(); ++mfi) {
      Box vbox = mfi.validbox();
      BoxArray cfba = BoxArray(gridArray[lev+1]).coarsen(refRatio[lev]);
      std::vector< std::pair<int,Box> > isects = cfba.intersections(vbox);
      for (int i = 0; i < isects.size(); i++) {
	Box ovlp = isects[i].second & vbox;
	if (ovlp.ok()) {
	  mft[lev].setVal(value,ovlp,0,nComp);
	}
      }
    }
  }
}

bool
MFTower::IsCompatible(const MFTower& rhs) const
{
  BL_PROFILE("MFTower::IsCompatible()");
    int numLevs = rhs.NumLevels();
    bool isok = numLevs<=rhs.GetLayout().NumLevels();
    const Array<IntVect>& refRatio = layout.RefRatio();
    const Array<BoxArray>& gridArray = layout.GridArray();
    for (int lev=0; lev<numLevs && isok; ++lev)
    {
        isok &= gridArray[lev] == rhs.GetLayout().GridArray()[lev];
        if (lev < numLevs-1) {
            isok &= rhs.GetLayout().RefRatio()[lev]==refRatio[lev];
        }
    }
    return isok;
}

void
MFTower::CCtoECgrad(PArray<MFTower>& mfte,
                    const MFTower&   mftc,
                    Real             mult,
                    int              sComp,
                    int              dComp,
                    int              nComp,
		    int              numLevs)
{
  BL_PROFILE("MFTower::CCtoECgrad()");
    // Note: Assumes that grow cells of mftc have been filled "properly":
    //    f-f: COI
    //    c-f: Parallel interp of coarse data to plane parallel to c-f, then perp interp to fine cc in grow
    //   phys: Dirichlet data at wall extrapolated to fine cc in grow
    //
    int numLevs_tmp = mftc.NumLevels();
    const Layout& theLayout = mftc.GetLayout();
    for (int d=0; d<BL_SPACEDIM; ++d) {
        if (numLevs<0) numLevs_tmp = std::min(numLevs_tmp,mfte[d].NumLevels());
        BL_ASSERT(theLayout.IsCompatible(mfte[d]));
    }
    if (numLevs<0) numLevs = numLevs_tmp;
    BL_ASSERT(numLevs<=mftc.NumLevels());
    const Array<Geometry>& geomArray = theLayout.GeomArray();
    for (int lev=0; lev<numLevs; ++lev) {
        const MultiFab& mfc = mftc[lev];
        BL_ASSERT(mfc.nGrow()>=1);
        BL_ASSERT(sComp+nComp<=mfc.nComp());

        const Real* dx = geomArray[lev].CellSize();
        for (MFIter mfi(mfc); mfi.isValid(); ++mfi) {
            const FArrayBox& cfab = mfc[mfi];
            const Box& vcbox = mfi.validbox();

            for (int d=0; d<BL_SPACEDIM; ++d) {
                FArrayBox& efab = mfte[d][lev][mfi];
                BL_ASSERT(dComp+nComp<=efab.nComp());
                efab.setVal(0);
                BL_ASSERT(Box(vcbox).surroundingNodes(d).contains(efab.box()));

                FORT_CC_TO_EC_GRAD(efab.dataPtr(dComp),ARLIM(efab.loVect()), ARLIM(efab.hiVect()),
                                   cfab.dataPtr(sComp),ARLIM(cfab.loVect()), ARLIM(cfab.hiVect()),
                                   vcbox.loVect(), vcbox.hiVect(),dx,&mult,&d,&nComp);
            }
        }
    }
}

void
MFTower::CCtoECavg(PArray<MFTower>& mfte,
                   const MFTower&   mftc,
                   Real             mult,
                   int              sComp,
                   int              dComp,
                   int              nComp,
                   int              do_harmonic,
                   int              numLevs)
{
  BL_PROFILE("MFTower::CCtoECavg()");
    // Note: Assumes that grow cells of mftc have been filled "properly":
    //    f-f: COI
    //    c-f: Parallel interp of coarse data to plane parallel to c-f, then perp interp to fine cc in grow
    //   phys: Dirichlet data at wall extrapolated to fine cc in grow
    //
    int numLevs_tmp = mftc.NumLevels();
    const Layout& theLayout = mftc.GetLayout();
    for (int d=0; d<BL_SPACEDIM; ++d) {
        if (numLevs<0) numLevs_tmp = std::min(numLevs_tmp,mfte[d].NumLevels());
        BL_ASSERT(theLayout.IsCompatible(mfte[d]));
    }
    if (numLevs<0) numLevs = numLevs_tmp;
    BL_ASSERT(numLevs<=mftc.NumLevels());
    const Array<Geometry>& geomArray = theLayout.GeomArray();
    for (int lev=0; lev<numLevs; ++lev) {
        const MultiFab& mfc = mftc[lev];
        BL_ASSERT(mfc.nGrow()>=1);

        for (MFIter mfi(mfc); mfi.isValid(); ++mfi) {
            const FArrayBox& cfab = mfc[mfi];
            const Box& vcbox = mfi.validbox();

            for (int d=0; d<BL_SPACEDIM; ++d) {
                FArrayBox& efab = mfte[d][lev][mfi];
                BL_ASSERT(dComp<efab.nComp());
                efab.setVal(0);
                BL_ASSERT(Box(vcbox).surroundingNodes(d).contains(efab.box()));

                int src_comp = (nComp<0 ? sComp + d : sComp);
                int num_comp = (nComp<0 ? 1 : nComp);
                BL_ASSERT(src_comp+num_comp<=mfc.nComp());
                FORT_CC_TO_EC_AVG(efab.dataPtr(dComp),   ARLIM(efab.loVect()), ARLIM(efab.hiVect()),
                                  cfab.dataPtr(src_comp),ARLIM(cfab.loVect()), ARLIM(cfab.hiVect()),
                                  vcbox.loVect(), vcbox.hiVect(),&mult,&d,&num_comp,&do_harmonic);
            }
        }
    }
}

void
MFTower::ECtoCCdiv(MFTower&               mftc,
                   const PArray<MFTower>& mfte,
                   const Array<Real>&     mult,
                   int                    sComp,
                   int                    dComp,
                   int                    nComp,
		   int                    numLevs)
{
  BL_PROFILE("MFTower::ECtoCCdiv()");
    if (numLevs<0) numLevs = mftc.NumLevels();
    const Layout& theLayout = mftc.GetLayout();
    BL_ASSERT(theLayout.NumLevels()>=numLevs);
    BL_ASSERT(mftc.NumLevels()>=numLevs);
    for (int d=0; d<BL_SPACEDIM; ++d) {
        BL_ASSERT(theLayout.IsCompatible(mfte[d]));
        BL_ASSERT(mfte[d].NumLevels()>=numLevs);
    }

    const Array<Geometry>& geomArray = theLayout.GeomArray();
    for (int lev=0; lev<numLevs; ++lev) {
        MultiFab& mfc = mftc[lev];
        BL_ASSERT(dComp+nComp<=mfc.nComp());

        const MultiFab& Vol = theLayout.Volume(lev);

        for (MFIter mfi(mfc); mfi.isValid(); ++mfi) {
            FArrayBox& cfab = mfc[mfi];
            const FArrayBox& vol = Vol[mfi];
            const Box& vcbox = mfi.validbox();

            cfab.setVal(0);

            for (int d=0; d<BL_SPACEDIM; ++d) {
                const FArrayBox& efab = mfte[d][lev][mfi];
                const FArrayBox& area = theLayout.Area(lev,d)[mfi];
                BL_ASSERT(sComp+nComp<=efab.nComp());
                BL_ASSERT(Box(vcbox).surroundingNodes(d).contains(efab.box()));

                FORT_EC_TO_CC_DIV(cfab.dataPtr(dComp),ARLIM(cfab.loVect()), ARLIM(cfab.hiVect()),
                                  efab.dataPtr(sComp),ARLIM(efab.loVect()), ARLIM(efab.hiVect()),
                                  area.dataPtr(),     ARLIM(area.loVect()), ARLIM(area.hiVect()),
                                  vol.dataPtr(),      ARLIM(vol.loVect()),  ARLIM(vol.hiVect()),
                                  vcbox.loVect(), vcbox.hiVect(), mult.dataPtr(), &d, &nComp);
            }
        }
    }
}

void
MFTower::AverageDown(MFTower& mft,
                     int      sComp,
                     int      nComp,
		     int      numLevs)
{
  BL_PROFILE("MFTower::AverageDown()");
    if (numLevs<0) numLevs = mft.NumLevels();
    BL_ASSERT(mft.NumLevels()>=numLevs);

    int dir = -1;
    for (int d=0; d<BL_SPACEDIM; ++d) {
        if (mft.ixType() == EC[d]) {
            dir = d;
        }
    }
    if (dir<0) {
        BoxLib::Error("MFTower::AverageDown not implemented for this centering yet");
    }

    const Layout& theLayout = mft.GetLayout();
    const Array<IntVect>& refRatio = theLayout.RefRatio();
    const Array<BoxArray>& gridArray = theLayout.GridArray();
    FArrayBox cfab;
    for (int lev=numLevs-1; lev>0; --lev) {
        MultiFab& mfe = mft[lev];
        BL_ASSERT(sComp+nComp<=mfe.nComp());
        const IntVect& ratio = refRatio[lev-1];

        // Build a container to hold the coarsened data on the fine distribution
        BoxArray eba = BoxArray(gridArray[lev]).coarsen(ratio).surroundingNodes(dir);
        MultiFab cemf(eba,nComp,0);
        cemf.setVal(0);

        for (MFIter mfi(mfe); mfi.isValid(); ++mfi) {
            const FArrayBox& fefab = mfe[mfi];
            FArrayBox& cefab = cemf[mfi];
            const Box& cebox = cefab.box();

            FORT_COARSEN_EC(cefab.dataPtr(sComp),ARLIM(cefab.loVect()), ARLIM(cefab.hiVect()),
                            fefab.dataPtr(sComp),ARLIM(fefab.loVect()), ARLIM(fefab.hiVect()),
                            cebox.loVect(), cebox.hiVect(), ratio.getVect(), &dir, &nComp);

        }

        mft[lev-1].copy(cemf); // parallel copy
    }
}

MFTFillPatch::MFTFillPatch(Layout& _layout)
    : layout(_layout), nLevs(_layout.NumLevels())
{
}



/*
     BuildInterpCoefs:

     This routine returns the Lagrange interpolating coefficients for a
     polynomial through N points, evaluated at xInt=-1 (see Numerical Recipes,
     v2, p102, e.g.):

             (x-x2)(x-x3)...(x-xN)              (x-x1)(x-x2)...(x-x(N-1))
     P(x) = ----------------------- y1  + ... + ------------------------  yN
            (x1-x2)(x1-x3)...(x1-xN)            (x1-x2)(x1-x3)...(x1-xN)

      P(xInt) = sum_(i=1)^(N) y[i]*c[i]
*/

void
BuildInterpCoefs(Real xVal, Array<Real>& coefs)
{
    int N = coefs.size();
    Array<Real> x(N);
    x[0] = -0.5-xVal; // xVal is location of off-grid Dirichlet value
    for (int i=1; i<N; ++i) {
        x[i] = i-1;
    }
    for (int j=0; j<N; ++j) {
        Real num = 1;
        Real den = 1;
        for (int i=0; i<N; ++i) {
            if (i!=j) {
                num *= -1 - x[i];
                den *= x[j] - x[i];
            }
        }
        coefs[j] = num/den;
    }
}

void
MFTFillPatch::BuildCFParallelInterpStencil()
{
    if (myhash.maxorder>=0) {
        return; // We have already made this stencil and it is independent of maxorder
    }

    // Some handy intvects
    Array<IntVect> ivp(BL_SPACEDIM), ivpp(BL_SPACEDIM), ivm(BL_SPACEDIM), ivmm(BL_SPACEDIM);
    for (int d=0; d<BL_SPACEDIM; ++d) {
        ivp[d] = BoxLib::BASISV(d);
        ivpp[d] = ivp[d] + BoxLib::BASISV(d);
        ivm[d] = -BoxLib::BASISV(d);
        ivmm[d] = ivm[d] - BoxLib::BASISV(d);
    }

    const PArray<Layout::MultiNodeFab>& nodes = layout.Nodes();
    const PArray<Layout::MultiNodeFab>& crseNodes = layout.CrseNodes();
    const Array<BoxArray>& bndryCells = layout.BndryCells();
    parallelInterpStencil.resize(nLevs);

    for (int lev=1; lev<nLevs; ++lev)
    {
        const Geometry& gl = layout.GeomArray()[lev];
        const IntVect& refRatio = layout.RefRatio()[lev-1];
        Array<IVSMap>& parInterpLev = parallelInterpStencil[lev];
        parInterpLev.resize(BL_SPACEDIM);
        for (MFIter mfi(nodes[lev]); mfi.isValid(); ++mfi)
        {
            const Layout::NodeFab& nodeFab = nodes[lev][mfi];
            const Box& cgbox = crseNodes[lev][mfi].box();
            for (int d=0; d<BL_SPACEDIM; ++d) {
                IVSMap& parInterpLevDir = parInterpLev[d];
                Array<int> dtan;
                for (int d0=0; d0<BL_SPACEDIM; ++d0) {
                    if (d!=d0) {
                        dtan.push_back(d0);
                    }
                }

                Box gdbox = Box(mfi.validbox()).grow(d,1) & gl.Domain();
                std::vector< std::pair<int,Box> > isects = bndryCells[lev].intersections(gdbox);
                for (int i=0; i<isects.size(); ++i) {
                    const Box& bndrySect = isects[i].second;
                    for (IntVect iv=bndrySect.smallEnd(), End=bndrySect.bigEnd(); iv<=End; bndrySect.next(iv)) {

                        Node nC = nodeFab(iv,0);
                        BL_ASSERT(cgbox.contains(nC.iv) && nC.type==Node::VALID);
                        std::map<int,Real> x;

                        Stencil& stencil = parInterpLevDir[iv];
                        stencil[nC] = 1;

                        for (int d0=0; d0<dtan.size(); ++d0) {
                            int itan = dtan[d0];
                            Stencil der, der2;
                            int r = refRatio[itan];

                            x[itan] = (iv[itan]%r - 0.5*(r-1))/r;

                            BL_ASSERT(crseNodes[lev][mfi].box().contains(nC.iv+ivp[itan]));
                            BL_ASSERT(crseNodes[lev][mfi].box().contains(nC.iv+ivm[itan]));
                            const Node& nR = crseNodes[lev][mfi](nC.iv+ivp[itan],0);
                            const Node& nL = crseNodes[lev][mfi](nC.iv+ivm[itan],0);

                            bool Rvalid = nR.type==Node::VALID;
                            bool Lvalid = nL.type==Node::VALID;

                            if (Rvalid && Lvalid) {
                                // Centered full
                                der[nL]  = -0.5;  der[nR] = +0.5;
                                der2[nL] = +1.0; der2[nC] = -2.0; der2[nR] = +1.0;
                            }
                            else if (Rvalid && crseNodes[lev][mfi].box().contains(nC.iv+ivpp[itan])) {
                                const Node& nRR = crseNodes[lev][mfi](nC.iv+ivpp[itan],0);
                                bool RRvalid = nRR.type==Node::VALID;
                                if (RRvalid) {
                                    // R-shifted full
                                    der[nC]  = -0.5;  der[nRR] = +0.5;
                                    der2[nC] = +1.0; der2[nR]  = -2.0; der2[nRR] = +1.0;
                                } else {
                                    // R-shifted linear
                                    der[nC] = -1.0; der[nR] = +1.0;
                                }
                            } else if (Lvalid && crseNodes[lev][mfi].box().contains(nC.iv+ivmm[itan])) {
                                const Node& nLL = crseNodes[lev][mfi](nC.iv+ivmm[itan],0);
                                bool LLvalid = nLL.type==Node::VALID;
                                if (LLvalid) {
                                    // L-shifted full
                                    der[nLL]  = -0.5;  der[nC] = +0.5;
                                    der2[nLL] = +1.0; der2[nL] = -2.0; der2[nC] = +1.0;
                                } else {
                                    // L-shifted linear
                                    der[nL] = -1.0; der[nC] = +1.0;
                                }
                            } else {
                                // piecewise constant (no derivatives)
                            }
                            stencil += x[itan]*der + 0.5*x[itan]*x[itan]*der2;
                        } // tangential direction

                        if (dtan.size()==2) {
                            const Node& npp = crseNodes[lev][mfi](nC.iv+ivp[dtan[0]]+ivp[dtan[1]],0);
                            const Node& nmm = crseNodes[lev][mfi](nC.iv+ivm[dtan[0]]+ivm[dtan[1]],0);
                            const Node& npm = crseNodes[lev][mfi](nC.iv+ivp[dtan[0]]+ivm[dtan[1]],0);
                            const Node& nmp = crseNodes[lev][mfi](nC.iv+ivm[dtan[0]]+ivp[dtan[1]],0);

                            bool PPvalid = npp.type==Node::VALID;
                            bool MMvalid = nmm.type==Node::VALID;
                            bool PMvalid = npm.type==Node::VALID;
                            bool MPvalid = nmp.type==Node::VALID;

                            if (PPvalid && MMvalid && PMvalid && MPvalid) {
                                Stencil crossterm;
                                crossterm[npp] = +0.25;
                                crossterm[nmm] = +0.25;
                                crossterm[npm] = -0.25;
                                crossterm[nmp] = -0.25;
                                stencil += x[dtan[0]]*x[dtan[1]]*crossterm;
                            }
                        }
                    } // bndry iv
                } // bndry box
            } // bc direction
        } // fine box
    } // lev
}

Layout&
MFTFillPatch::GetLayout()
{
    return layout;
}

MFTFillPatch::MyHash::MyHash(const BCRec& _bc, int _maxorder)
{
    maxorder = _maxorder;
    bc = BCRec(_bc.lo(),_bc.hi());
}

bool
operator==(const MFTFillPatch::MyHash& lhs, const MFTFillPatch::MyHash& rhs)
{
    if (lhs.maxorder!=rhs.maxorder) return false;

    for (int d=0; d<BL_SPACEDIM; ++d) {
        if (lhs.bc.lo(d) != rhs.bc.lo(d) ) return false;
        if (lhs.bc.hi(d) != rhs.bc.hi(d) ) return false;
    }
    return true;
}

bool
operator!=(const MFTFillPatch::MyHash& lhs, const MFTFillPatch::MyHash& rhs)
{
    return !operator==(lhs,rhs);
}

int // Return 0 if successful, 1 if bc not properly defined
MFTFillPatch::BuildStencil(const BCRec& bc,
                          int maxorder)
{
    if (MyHash(bc,maxorder)==myhash) {
        return 0; // we have already made what we need
    }

    int myproc = ParallelDescriptor::MyProc();

    BuildCFParallelInterpStencil();
    const PArray<Layout::MultiNodeFab>& nodes = layout.Nodes();

    // Precompute an often-used interp stencil
    Array<Real> iCoefsZero(maxorder); BuildInterpCoefs(0,iCoefsZero); // value at wall
    Array<Real> iCoefsFO(1); BuildInterpCoefs(0.5,iCoefsFO); // value at grow center (FOEXTRAP)
    Array<Real> iCoefsHO(maxorder); BuildInterpCoefs(0.5,iCoefsHO); // value at grow center (HOEXTRAP)
    Array<Real> iCoefsRE(1,1); // FIXME: Need to go to maxorder
    Array<Real> iCoefsRO(1,-1);// FIXME: Need to go to maxorder

    Array<Array<Real> > iCoefsCF(BL_SPACEDIM, Array<Real>(maxorder));
    const Array<Geometry>& geomArray = layout.GeomArray();
    const Array<IntVect>& refRatio = layout.RefRatio();
    const Array<BoxArray>& gridArray = layout.GridArray();

    perpInterpStencil.resize(nLevs,Array<IVSMap>(BL_SPACEDIM));
    for (int lev=0; lev<nLevs; ++lev) {
        Array<IVSMap>& perpInterpLev = perpInterpStencil[lev];
        const Box& dbox = geomArray[lev].Domain();
        Array<Array<Box> > bndry(2, Array<Box>(BL_SPACEDIM));
        for (int d=0; d<BL_SPACEDIM; ++d) {
            bndry[0][d] = BoxLib::adjCellLo(dbox,d);
            bndry[1][d] = BoxLib::adjCellHi(dbox,d);
        }

        // Precompute other often-used interp stencils
        if (lev>0) {
            for (int d=0; d<BL_SPACEDIM; ++d) {
                const IntVect& rat = refRatio[lev-1];
                BuildInterpCoefs(0.5*rat[d],iCoefsCF[d]); // value at wall
            }
        }

        const BoxArray& ba = gridArray[lev];
        const Layout::MultiNodeFab& fmn = nodes[lev];
        const Geometry& gl = geomArray[lev];
        for (MFIter mfi(fmn); mfi.isValid(); ++mfi) {
            const Box& vbox = mfi.validbox();
            const Layout::NodeFab& fn = fmn[mfi];
            for (int d=0; d<BL_SPACEDIM; ++d) {

                IVSMap& perpInterpLevDir = perpInterpLev[d];

                Array<Box> myBndry(2);
                myBndry[0] = BoxLib::adjCellLo(vbox,d);
                myBndry[1] = BoxLib::adjCellHi(vbox,d);

                for (int hilo=0; hilo<2; ++hilo)
                {
                    Box povlp = myBndry[hilo] & bndry[hilo][d];
                    int bc_flag = (hilo==0 ? bc.lo()[d] : bc.hi()[d]);
                    if (povlp.ok()) {

                        Array<Real>* interpCoef;
                        if (bc_flag==EXT_DIR) {
                            interpCoef = &iCoefsZero;
                        }
                        else if (bc_flag==FOEXTRAP) {
                            interpCoef = &iCoefsFO;
                        }
                        else if (bc_flag==HOEXTRAP) {
                            interpCoef = &iCoefsHO;
                        }
                        else if (bc_flag==REFLECT_EVEN) {
                            interpCoef = &iCoefsRE;
                        }
                        else if (bc_flag==REFLECT_ODD) {
                            interpCoef = &iCoefsRO;
                        }
                        else {
			    return 1;
                        }

                        int sgn = (hilo==0 ? +1  : -1); // Direction of interp stencil (inward)
                        for (IntVect iv=povlp.smallEnd(), End=povlp.bigEnd(); iv<=End; povlp.next(iv)) {
                            Stencil& stencil = perpInterpLevDir[iv];

                            int icnt = 0;
                            if (bc_flag == EXT_DIR) {
                                Node n = fn(iv,0); // This will have been an invalid node until now
                                n.level = lev; n.iv = iv; n.type=Node::VALID;
                                stencil[n] = (*interpCoef)[icnt++];
                            }
                            for (int k=0; icnt<interpCoef->size(); ++k, ++icnt) {
                                IntVect siv = iv + sgn*(k+1)*BoxLib::BASISV(d);
                                stencil[fn(siv,0)] = (*interpCoef)[icnt];
                            }
                        }
                    }
                    else if (lev>0) {
                        // Build c-f stencil
                        BoxArray sba = BoxLib::complementIn(myBndry[hilo],ba);
                        if (gl.isPeriodic(d)) {
                            BoxArray per_ba = BoxLib::intersect(sba,BoxArray(ba).shift(d,dbox.length(d)));
                            for (int j=0; j<per_ba.size(); ++j) {
                                sba = BoxLib::complementIn(per_ba[j],sba);
                            }
                            per_ba = BoxLib::intersect(sba,BoxArray(ba).shift(d,-dbox.length(d)));
                            for (int j=0; j<per_ba.size(); ++j) {
                                sba = BoxLib::complementIn(per_ba[j],sba);
                            }
                        } else {
                            sba = BoxLib::intersect(sba,dbox);
                        }

                        // Now, with coefs create stencil entries
                        const IVSMap& parStencil = parallelInterpStencil[lev][d];
                        for (int j=0; j<sba.size(); ++j) {
                            const Box& sbox = sba[j];
                            for (IntVect iv=sbox.smallEnd(), End=sbox.bigEnd(); iv<=End; sbox.next(iv)) {

                                IVScit it = parStencil.find(iv);
                                BL_ASSERT(it!=parStencil.end());
                                const Stencil& parallelStencil = it->second;

                                Stencil& stencil = perpInterpLevDir[iv];
                                stencil = iCoefsCF[d][0]*parallelStencil;
                                int sgn = (hilo==0 ? +1  : -1); // Direction of interp stencil (inward)
                                for (int k=1; k<iCoefsCF[d].size(); ++k) {
                                    IntVect siv = iv + sgn*k*BoxLib::BASISV(d);
                                    stencil[fn(siv,0)] = iCoefsCF[d][k];
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    myhash = MyHash(bc,maxorder);
    return 0;
}

void
MFTFillPatch::DoCoarseFineParallelInterp(MFTower& mft,
					 int      sComp,
					 int      nComp,
					 int      numLevs) const
{
    if (numLevs<0) numLevs = nLevs;
    BL_ASSERT(layout.IsCompatible(mft));
    BL_ASSERT(numLevs<=layout.NumLevels());
    BL_ASSERT(numLevs<=mft.NumLevels());
    const Array<BoxArray>& bndryCells = layout.BndryCells();

    const Array<Geometry>& geomArray = layout.GeomArray();
    const Array<IntVect>& refRatio = layout.RefRatio();
    const Array<BoxArray>& gridArray = layout.GridArray();
    for (int lev=1; lev<numLevs; ++lev) {

        MultiFab& mf = mft[lev];
        BL_ASSERT(mf.nGrow()>=1);
        BL_ASSERT(sComp+nComp<=mf.nComp());

        BoxArray bacgd = BoxArray(mf.boxArray()).coarsen(refRatio[lev-1]).grow(1);
        MultiFab crseMF(bacgd,nComp,0);
        crseMF.copy(mft[lev-1],sComp,0,nComp); // parallel copy

        BoxArray bnd = bndryCells[lev];
        const Geometry& gl = geomArray[lev];
        const Array<IVSMap>& parInterpLev = parallelInterpStencil[lev];

        for (MFIter mfi(mf); mfi.isValid(); ++mfi) {

            int gridIdx = mfi.index();

            FArrayBox& crseFab = crseMF[mfi];
            FArrayBox& fineFab = mf[mfi];

            for (int d=0; d<BL_SPACEDIM; ++d) {

                const IVSMap& parInterp = parInterpLev[d];
                Box boxgd = Box(mfi.validbox()).grow(d,1) & gl.Domain();
                std::vector< std::pair<int,Box> > isects = bnd.intersections(boxgd);
                for (int i=0; i<isects.size(); ++i) {
                    const Box& bndrySect = isects[i].second;
                    for (IntVect iv=bndrySect.smallEnd(), End=bndrySect.bigEnd(); iv<=End; bndrySect.next(iv)) {


                        IVScit it=parInterp.find(iv);
                        if (it!=parInterp.end()) {
                            const Stencil& s = it->second;

                            for (int n=0; n<nComp; ++n) {
                                Real res = 0;
                                for (Stencil::const_iterator it=s.begin(), End=s.end(); it!=End; ++it) {
                                    const IntVect& ivs=(it->first).iv;
                                    BL_ASSERT(crseFab.box().contains(ivs));
                                    BL_ASSERT((it->first).level==lev-1);
                                    res += crseFab(ivs,sComp+n) * it->second;
                                }
                                fineFab(iv,sComp+n) = res;
                            }
                        }
                    }
                }
            }
        }

        mf.FillBoundary(sComp,nComp);
    }
}

#include <iomanip>
void
MFTFillPatch::FillGrowCells(MFTower& mft,
                            int      sComp,
                            int      nComp,
                            bool     do_piecewise_constant,
			    int      numLevs) const
{
    if (numLevs<0) numLevs=nLevs;
    if (do_piecewise_constant) {
      FillGrowCellsSimple(mft,sComp,nComp,numLevs);
        return;
    }

    // Note: Assumes that Dirichlet values have been loaded into the MultiFab grow cells at
    // physical boundaries (and that these values are to be applied at the cell walls).
    BL_ASSERT(layout.IsCompatible(mft));
    BL_ASSERT(numLevs<=layout.NumLevels());

    BL_ASSERT(perpInterpStencil.size()>=numLevs);
    const Array<Geometry>& geomArray = layout.GeomArray();
    const Array<IntVect>& refRatio = layout.RefRatio();
    const Array<BoxArray>& gridArray = layout.GridArray();
    for (int lev=0; lev<numLevs; ++lev) {
        const Array<IVSMap>& perpInterpLev = perpInterpStencil[lev];
        MultiFab& mf = mft[lev];
        BL_ASSERT(mf.nGrow()>=1);
        BL_ASSERT(sComp+nComp<=mf.nComp());

        MultiFab crseMF;
        if (lev>0) {
            BoxArray bacgd = BoxArray(mf.boxArray()).coarsen(refRatio[lev-1]).grow(1);
            crseMF.define(bacgd,nComp,0,Fab_allocate);
            crseMF.copy(mft[lev-1],sComp,0,nComp); // parallel copy (maybe excessive data copied?)
        }

        const Geometry& gl = geomArray[lev];
        BoxArray bnd = Layout::GetBndryCells(gridArray[lev],IntVect::TheUnitVector(),gl); // Note: layout Bndrycells excluded phys

        for (MFIter mfi(mf); mfi.isValid(); ++mfi) {

            FArrayBox* crseFab;
            if (lev>0) {
                crseFab = &(crseMF[mfi]);
            }
            FArrayBox& fineFab = mf[mfi];

            // Fill grow cells, but only the ones for which we have a stencil
            for (int d=0; d<BL_SPACEDIM; ++d) {
                const IVSMap& perpInterpLevDir = perpInterpLev[d];
                Box boxdg = Box(mfi.validbox()).grow(d,1);

                std::vector< std::pair<int,Box> > isects = bnd.intersections(boxdg);
                for (int i=0; i<isects.size(); ++i) {
                    const Box& bndrySect = isects[i].second;
                    for (IntVect iv=bndrySect.smallEnd(), End=bndrySect.bigEnd(); iv<=End; bndrySect.next(iv)) {
                        IVScit it=perpInterpLevDir.find(iv);
                        if (it!=perpInterpLevDir.end()) {
			    const Stencil& s = it->second;
                            for (int n=0; n<nComp; ++n) {
                                Real res = 0;
                                for (Stencil::const_iterator it=s.begin(), End=s.end(); it!=End; ++it) {
                                    const IntVect& ivs=(it->first).iv;
                                    int slev = (it->first).level;
                                    if (slev==lev) {
                                        BL_ASSERT(fineFab.box().contains(ivs));
                                        res += fineFab(ivs,sComp+n) * it->second;
                                    }
                                    else if (slev==lev-1) {
                                        BL_ASSERT(crseFab->box().contains(ivs));
                                        res += (*crseFab)(ivs,sComp+n) * it->second;
                                    }
                                }
                                fineFab(iv,sComp+n) = res;
                            }
                        }
                    }
                }
            }
        }

        mf.FillBoundary(sComp,nComp);
    }
}

void
MFTFillPatch::FillGrowCellsSimple(MFTower& mft,
				  int      sComp,
				  int      nComp,
				  int      numLevs) const
{
    // Note: Assumes that Dirichlet values have been loaded into the MultiFab grow cells at
    // physical boundaries (and that these values are to be applied at the cell walls).
    //
    // The simple version effectively does the full stuff, but with a first order interpolant
    // both parallel and perpendicular....ie, piecewise constant.  Since this is trivial to
    // construct, we do not bother with the heavy guns
    if (numLevs<0) numLevs=nLevs;
    BL_ASSERT(layout.IsCompatible(mft));
    BL_ASSERT(layout.NumLevels()>=numLevs);
    BL_ASSERT(mft.NumLevels()>=numLevs);

    const PArray<Layout::MultiNodeFab>& nodes = layout.Nodes();
    const Array<Geometry>& geomArray = layout.GeomArray();
    const Array<IntVect>& refRatio = layout.RefRatio();
    const Array<BoxArray>& gridArray = layout.GridArray();
    for (int lev=0; lev<numLevs; ++lev) {

        MultiFab& mf = mft[lev];
        BL_ASSERT(mf.nGrow()>=1);
        BL_ASSERT(sComp+nComp<=mf.nComp());

        MultiFab crseMF;
        if (lev>0) {
            BoxArray bacgd = BoxArray(mf.boxArray()).coarsen(refRatio[lev-1]).grow(1);
            crseMF.define(bacgd,nComp,0,Fab_allocate);
            crseMF.copy(mft[lev-1],sComp,0,nComp); // parallel copy (maybe excessive data copied?)
        }

        const Geometry& gl = geomArray[lev];
        BoxArray bnd = Layout::GetBndryCells(gridArray[lev],IntVect::TheUnitVector(),gl); // Note: layout Bndrycells excluded phys
        const Layout::MultiNodeFab& nodesLev = nodes[lev];


        if (lev>0) {
            for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
                FArrayBox& crseFab = crseMF[mfi];
                FArrayBox& fineFab = mf[mfi];
                const Layout::NodeFab& nodesFab = nodesLev[mfi];

                // fill c-f grow cells (physbc values already assumed to be in place
                for (int d=0; d<BL_SPACEDIM; ++d) {
                    Box boxdg = Box(mfi.validbox()).grow(d,1)  &  gl.Domain();
                    std::vector< std::pair<int,Box> > isects = bnd.intersections(boxdg);
                    for (int i=0; i<isects.size(); ++i) {
                        const Box& bndrySect = isects[i].second;
                        for (IntVect iv=bndrySect.smallEnd(), End=bndrySect.bigEnd(); iv<=End; bndrySect.next(iv)) {
                            const Node& node = nodesFab(iv,0);
                            int slev = node.level;
                            BL_ASSERT(slev == lev-1);
                            const IntVect& ivs = node.iv;
                            BL_ASSERT(crseFab.box().contains(ivs));
                            for (int n=0; n<nComp; ++n) {
                                fineFab(iv,sComp+n) = crseFab(ivs,sComp+n);
                            }
                        }
                    }
                }
            }
        }

        mf.FillBoundary(sComp,nComp);
    }
}

void
MFTower::SetVal(Real     val,
		int      sComp,
		int      nComp,
		int      numLevs)
{
  BL_PROFILE("MFTower::SetVal()");
    if (numLevs<0) numLevs=nLevs;
    BL_ASSERT(numLevs <= nLevs);
    const Layout& layout = GetLayout();
    for (int lev=0; lev<numLevs; ++lev)
    {
        int nGrow = mft[lev].nGrow();
        mft[lev].setVal(val,sComp,nComp,nGrow);
    }
}

void
MFTower::Write(const std::string& fileName,
	       const std::string& varname,
	       Real               time) const
{
  BL_PROFILE("MFTower::Write()");
  Array<MFTower*> mfta(1,(MFTower*)this);
  Array<std::string> varnames(1,std::string(varname));
  if (NComp()>1) {
    varnames.resize(NComp());
    for (int i=0; i<NComp(); ++i) {
      varnames[i] = BoxLib::Concatenate(varname,i,2);
    }
  }
  MFTower::WriteSet(fileName,mfta,varnames,time);
}

void
MFTower::WriteSet(const std::string&          fileName,
		  const Array<MFTower*> mfta,
		  const Array<std::string>&   varnames,
		  Real                        time)
{
  BL_PROFILE("MFTower::WriteSet()");
  std::string pfversion = "MFTowerData-0.1";
  int ncomps = 0;
  for (int i=0; i<mfta.size(); ++i) {
    ncomps += mfta[i]->NComp();
  }
  BL_ASSERT(ncomps>0 && ncomps==varnames.size());
  const MFTower& mft = *(mfta[0]);
  for (int i=1; i<mfta.size(); ++i) {
    BL_ASSERT(mfta[i]->IsCompatible(mft));
  }

  const Layout& layout = mft.GetLayout();
  const Array<Geometry>& geomArray = layout.GeomArray();
  const Array<IntVect>& refRatio = layout.RefRatio();
  int num_levels = mft.NumLevels();

  // Currently, amrvis only supports integer refinement ratios
  BL_ASSERT(refRatio.size()==num_levels-1);
  Array<int> irat(num_levels-1);
  for (int i=0; i<num_levels-1; ++i) {
    irat[i] = refRatio[i][0];
    for (int j=1; j<BL_SPACEDIM; ++j) {
      BL_ASSERT(refRatio[i][j] == irat[i]);
    }
  }

  Array<Array<Real> > dxLevel(num_levels,Array<Real>(BL_SPACEDIM));
  Array<Box> probDomain(num_levels);
  for (int i=0; i<num_levels; ++i) {
    for (int j=0; j<BL_SPACEDIM; ++j) {
      dxLevel[i][j] = geomArray[i].CellSize()[j];
    }
    probDomain[i] = geomArray[i].Domain();
  }

  int coordSys = geomArray[0].IsCartesian();
  bool plt_verbose = false;
  bool isCartGrid = false;
  Real vfeps = 1.e-12;
  Array<int> levelSteps(num_levels,0);

  Array<Array<MultiFab*> > data(mfta.size(),Array<MultiFab*>(num_levels));
  for (int i=0; i<mfta.size(); ++i) {
    for (int lev=0; lev<num_levels; ++lev) {
      data[i][lev] = &( (*mfta[i])[lev] );
    }
  }
  WritePlotfile(pfversion,data,time,geomArray[0].ProbLo(),geomArray[0].ProbHi(),irat,probDomain,
		dxLevel,coordSys,fileName,varnames,plt_verbose,isCartGrid,&vfeps,levelSteps.dataPtr());

}
