/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <MFVector.H>

namespace Amanzi {

  MFVector::MFVector()
    : MultiFab() {}

  MFVector::MFVector(const MultiFab& mf, int scomp, int ncomp, int ngrow)
    : MultiFab(), srccomp(scomp) {
    numcomp = (ncomp < 0 ? mf.nComp() : std::min(ncomp,mf.nComp()));
    numgrow = (ngrow < 0 ? mf.nGrow() : std::min(ngrow,mf.nGrow()));
    srccomp = (scomp < 0 ? 0 : std::min(mf.nComp()-1,scomp));
    define(mf.boxArray(),numcomp,numgrow,mf.DistributionMap(),Fab_allocate);
    MultiFab::Copy(this->multiFab(),mf,srccomp,0,numcomp,numgrow);
  }

  MFVector::MFVector(const MFVector& vec, int scomp, int ncomp, int ngrow)
    : MultiFab(), srccomp(scomp) {
    numcomp = (ncomp < 0 ? vec.nComp() : std::min(ncomp,vec.nComp()));
    numgrow = (ngrow < 0 ? vec.nGrow() : std::min(ngrow,vec.nGrow()));
    srccomp = (scomp < 0 ? 0 : std::min(vec.nComp()-1,scomp));
    define(vec.boxArray(),numcomp,numgrow,vec.DistributionMap(),Fab_allocate);
    MultiFab::Copy(this->multiFab(),vec.multiFab(),srccomp,0,numcomp,numgrow);
  }

  void
  MFVector::AXPBY(const MFVector& Y, Real a, Real b) {// this = a * X  +  b * Y compoenentwise
    BL_ASSERT(nComp()==Y.nComp());
    BL_ASSERT(boxArray()==Y.boxArray());
    FArrayBox t;
    for (MFIter mfi(*this); mfi.isValid(); ++mfi) {
      FArrayBox& fabX = (*this)[mfi];
      const FArrayBox& fabY = Y[mfi];
      if (a!=1) {
        fabX.mult(a,0,nComp());
      }
      if (b!=0) {
        t.resize(mfi.validbox(),nComp());
        t.copy(fabY,0,0,nComp());
        if (b!=1) {
          t.mult(b,0,nComp());
        }
        fabX.plus(t,0,0,nComp());
      }
    }
  }

  void
  MFVector::AXPBYI(const MFVector& Y, Real a, Real b) {// this = a * X  +  b * (1/Y) compoenentwise
    BL_ASSERT(nComp()<=Y.nComp());
    BL_ASSERT(boxArray()==Y.boxArray());
    FArrayBox t;
    for (MFIter mfi(*this); mfi.isValid(); ++mfi) {
      FArrayBox& fabX = (*this)[mfi];
      const FArrayBox& fabY = Y[mfi];
      if (a!=1) {
        fabX.mult(a,0,nComp());
      }
      if (b!=0) {
        t.resize(mfi.validbox(),nComp());
        t.copy(fabY,0,0,nComp());
        t.invert(b,0,nComp());
        fabX.plus(t,0,0,nComp());
      }
    }
  }

  void
  MFVector::SCALE(Real a) { // this *= a  componentwise
    if (a!=1) {
      this->mult(a,0,nComp(),nGrow());
    }
  }

  void
  MFVector::MULTAY(const MFVector& Y, Real a) { // this *= a * Y componentwise
    if (a!=0) {
      BL_ASSERT(nComp()<=Y.nComp());
      BL_ASSERT(boxArray()==Y.boxArray());
      if (a==1) {
        MultiFab::Multiply(*this,Y,0,0,nComp(),nGrow());
      } else {
        FArrayBox t;
        for (MFIter mfi(*this); mfi.isValid(); ++mfi) {
          FArrayBox& fabX = (*this)[mfi];
          const FArrayBox& fabY = Y[mfi];
          t.resize(mfi.validbox(),nComp());
          t.copy(fabY,0,0,nComp());
          t.mult(a,0,nComp());
          fabX.mult(t,0,0,nComp());
        }
      }
    }
  }

  void
  MFVector::MULTAYI(const MFVector& Y, Real a) { // this *= a * (1/Y) componentwise
    if (a!=0) {
      FArrayBox t;
      for (MFIter mfi(*this); mfi.isValid(); ++mfi) {
        FArrayBox& fabX = (*this)[mfi];
        const FArrayBox& fabY = Y[mfi];
        t.resize(mfi.validbox(),nComp());
        t.copy(fabY,0,0,nComp());
        t.invert(a,0,nComp());
        fabX.mult(t,0,0,nComp());
      }
    }
  }

  void
  MFVector::COPY(const MFVector& rhs) {
    MultiFab::Copy(*this,rhs,0,0,nComp(),0);
  }

  void
  MFVector::INVERT(Real a) {
    this->invert(a,0,nGrow());
  }

} /* namespace Amanzi */
