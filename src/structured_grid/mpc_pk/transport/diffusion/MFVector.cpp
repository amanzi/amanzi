#include <MFVector.H>

namespace Amanzi {

  MFVector::MFVector(MultiFab& mf, int scomp, int ncomp, int ngrow)
    : MultiFab(), srccomp(scomp) {
    int num_comp = (ncomp < 0 ? mf.nComp() : std::min(ncomp,mf.nComp()));
    int num_grow = (ngrow < 0 ? mf.nGrow() : std::min(ngrow,mf.nGrow()));
    define(mf.boxArray(),num_comp,num_grow,mf.DistributionMap(),Fab_allocate);
  }

  MFVector*
  MFVector::Clone(int scomp, int ncomp, int ngrow) const {
    MultiFab *newmf = new MultiFab();
    int num_comp = (ncomp < 0 ? nComp() : std::max(ncomp,nComp()));
    int num_grow = (ngrow < 0 ? nGrow() : std::max(ngrow,nGrow()));
    newmf->define(boxArray(),num_comp,num_grow,DistributionMap(),Fab_allocate);
    return new MFVector(*newmf,scomp,ncomp,ngrow);
  }

  void
  MFVector::AXPBY(const MFVector& Y, Real a, Real b) {// this = a * X  +  b * Y compoenentwise
    BL_ASSERT(nComp()==Y.numComp());
    BL_ASSERT(boxArray()==Y.boxArray());
    FArrayBox t;
    for (MFIter mfi(*this); mfi.isValid(); ++mfi) {
      FArrayBox& fabX = (*this)[mfi];
      const FArrayBox& fabY = Y[mfi];
      if (a!=1) {
        fabX.mult(a,srcComp(),numComp());
      }
      if (b!=0) {
        t.resize(fabX.box(),numComp());
        t.copy(fabY,Y.srcComp(),0,numComp());
        if (b!=1) {
          t.mult(b,srcComp(),numComp());
        }
        fabX.plus(t,0,srcComp(),numComp());
      }
    }
  }

  void
  MFVector::AXPBYI(const MFVector& Y, Real a, Real b) {// this = a * X  +  b * (1/Y) compoenentwise
    BL_ASSERT(nComp()==Y.numComp());
    BL_ASSERT(boxArray()==Y.boxArray());
    FArrayBox t;
    for (MFIter mfi(*this); mfi.isValid(); ++mfi) {
      FArrayBox& fabX = (*this)[mfi];
      const FArrayBox& fabY = Y[mfi];
      if (a!=1) {
        fabX.mult(a,srcComp(),numComp());
      }
      if (b!=0) {
        t.resize(fabX.box(),numComp());
        t.copy(fabY,Y.srcComp(),0,numComp());
        t.invert(b,0,numComp());
        fabX.plus(t,0,srcComp(),numComp());
      }
    }
  }

  void
  MFVector::SCALE(Real a) { // this *= a  componentwise
    if (a!=1) {
      this->mult(a,srcComp(),numComp(),numGrow());
    }
  }

  void
  MFVector::MULTAY(const MFVector& Y, Real a) { // this *= a * Y componentwise
    if (a!=0) {
      if (a==1) {
        MultiFab::Multiply(*this,Y,Y.srcComp(),srcComp(),numComp(),numGrow());
      } else {
        FArrayBox t;
        for (MFIter mfi(*this); mfi.isValid(); ++mfi) {
          FArrayBox& fabX = (*this)[mfi];
          const FArrayBox& fabY = Y[mfi];
          t.resize(fabX.box(),numComp());
          t.copy(fabY,Y.srcComp(),0,numComp());
          if (a!=1) {
            t.mult(a,0,numComp());
          }
          fabX.plus(t,0,srcComp(),numComp());
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
        t.resize(fabX.box(),numComp());
        t.copy(fabY,Y.srcComp(),0,numComp());
        t.invert(a,0,numComp());
        fabX.mult(t,0,srcComp(),numComp());
      }
    }
  }

  void
  MFVector::COPY(const MFVector& rhs) {
    MultiFab::Copy(*this,rhs,rhs.srcComp(),srcComp(),numComp(),numGrow());
  }

  void
  MFVector::INVERT(Real a) {
    this->invert(a,srcComp(),numGrow());
  }

} /* namespace Amanzi */
