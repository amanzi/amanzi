/*
  Transport PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/
//! Simple adapter that changes units.

/*

This class is a simple adapter for TransportDomainFunctions, or potentially
other domain functions, that need a unit conversion.

It uses the CRTP to override the DomainFunction's Compute method, multiplying
(dividing) by unit conversion factors provided through a set_conversion()
method.

Input and output units are dependent upon context.  As an example, Amanzi's
TransportSourceFunction_* methods provide source terms in units of
[mol / L / s].  ATS wishes to use the same quantity, but converted to
[mol-C / mol-H2O / s].  The difference is a factor of:

  1000 [L / m^3] / molar_density_water [mol-H2O / m^3].

To do this, assuming the molar density of water is a field
(e.g. temperature-dependent) and not a scalar:

  Teuchos::RCP<Epetra_MultiVector> mol_dens_water = ...;
  TransportDomainFunction_UnitConversion<TransportSourceFunction_Alquimia> source(...);
  source.set_conversion(1000, mol_dens_water, true);
  source.Compute(t_old, t_new);

Note the third (boolean) argument reciprocal indicates to divide by the
_vector_ coefficient, but not the scalar one.

*/

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
