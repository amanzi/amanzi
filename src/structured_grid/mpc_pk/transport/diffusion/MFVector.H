#ifndef _MFVector_H_
#define _MFVector_H_

#include <MultiFab.H>

namespace Amanzi {

  struct MFVector
    : public MultiFab
  {
    MFVector(MultiFab& mf, int scomp=-1, int ncomp=-1, int ngrow=-1);
    const MultiFab& multiFab() const {return static_cast<const MultiFab&>(*this);}
    const int srcComp() const {return srccomp;}
    const int numComp() const {return numcomp;}
    const int numGrow() const {return numgrow;}
    MFVector* Clone(int scomp=-1, int ncomp=-1, int ngrow=-1) const;

    void AXPBY(const MFVector& Y, Real a=1, Real b=1); // this = a * X  +  b * Y compoenentwise
    void AXPBYI(const MFVector& Y, Real a, Real b); // this = a * X  +  b * (1/Y) compoenentwise
    void SCALE(Real a);  // this *= a  componentwise
    void MULTAY(const MFVector& Y, Real a=1); // this *= a * Y componentwise
    void MULTAYI(const MFVector& Y, Real a=1); // this *= a * (1/Y) componentwise
    void COPY(const MFVector& rhs); // this = Y componentwise
    void INVERT(Real a); // this = 1/(this) componentwise
  protected:
    int srccomp,numcomp,numgrow;
  };
  
} /* namespace Amanzi */
#endif /* _MFVector_H_ */
