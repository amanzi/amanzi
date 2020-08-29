/*
  Transport PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/
//! Simple adapter that changes units.

#pragma once

#include "TransportBoundaryFunction_Alquimia.hh"
#include "TransportSourceFunction_Alquimia.hh"

namespace Amanzi {
namespace Transport {

template<class T>
class TransportDomainFunction_UnitConversion : public T {
 public:
  using T::T;

  void Compute(double t_old, double t_new) {
    T::Compute(t_old, t_new);

    AMANZI_ASSERT(vector_ != Teuchos::null);
    const Epetra_MultiVector& vector = *vector_;

    if (reciprocal_) {
      for (auto& it : *this) {
        for (auto& val : it.second) {
          val *= scalar_ / vector[0][it.first];
        }
      }
    } else {
      for (auto& it : *this) {
        for (auto& val : it.second) {
          val *= scalar_ * vector[0][it.first];
        }
      }
    }
  }

  void set_conversion(double scalar,
                      const Teuchos::RCP<const Epetra_MultiVector>& vector,
                      bool reciprocal) {
    scalar_ = scalar;
    vector_ = vector;
    reciprocal_ = reciprocal;
  }

protected:
  double scalar_;
  Teuchos::RCP<const Epetra_MultiVector> vector_;
  bool reciprocal_;
};

using TransportBoundaryFunction_Alquimia_Units =
  TransportDomainFunction_UnitConversion<TransportBoundaryFunction_Alquimia>;
using TransportSourceFunction_Alquimia_Units =
  TransportDomainFunction_UnitConversion<TransportSourceFunction_Alquimia>;

} // namespace Transport
} // namespace Amanzi
