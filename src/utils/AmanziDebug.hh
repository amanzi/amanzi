#pragma once

#include "AmanziTypes.hh"
#include "AmanziVector.hh"

namespace Amanzi {
namespace Debug {

template<class Vector>
double get0(const Vector& v) { return -999.; }

template<> inline double get0(const Vector_type& v) {
  // auto vv = v.getLocalViewHost();
  // return vv(0,0);
  return v.normInf();
}

template<> inline double get0(const MultiVector_type& v) {
  // auto vv = v.getLocalViewHost();
  // return vv(0,0);
  return v.getVector(0)->normInf();
}


template<class Vector>
int isSame(const Vector& v1, const Vector& v2) { return -1; }

template<>
inline int isSame(const Vector_type& v1, const Vector_type& v2) {
  return int(v1.getLocalViewHost(Tpetra::Access::ReadOnly) == v2.getLocalViewHost(Tpetra::Access::ReadOnly));
}



} // namespace Debug
} // namespace Amanzi

