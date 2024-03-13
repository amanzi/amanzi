/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#ifndef _MFVector_H_
#define _MFVector_H_

#include <MultiFab.H>

namespace Amanzi {

  struct MFVector
    : public MultiFab
  {
    MFVector();
    MFVector(const MFVector& vec, int scomp=-1, int ncomp=-1, int ngrow=-1);
    MFVector(const MultiFab& mf, int scomp=-1, int ncomp=-1, int ngrow=-1);
    MultiFab& multiFab() {return static_cast<MultiFab&>(*this);}
    const MultiFab& multiFab() const {return static_cast<const MultiFab&>(*this);}

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
