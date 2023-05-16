/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  EOS

  EOS for tabulated density for temperature between 0.5 and 800 C and
  pressure between 634 Pa and 110 MPa
*/

#ifndef AMANZI_EOS_DENSITY_TABULAR_HH_
#define AMANZI_EOS_DENSITY_TABULAR_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "EOS_Density.hh"
#include "LookupTable.hh"

namespace Amanzi {
namespace AmanziEOS {

// Equation of State model
class DensityTabular : public EOS_Density {
 public:
  explicit DensityTabular(Teuchos::ParameterList& eos_plist);

  virtual double Density(double T, double p) override
  {
    double val = table_->Function(T, p, &ierr_);
    if (ierr_ != 0) error_msg_ = table_->ErrorMessage(T, p);
    return val;
  }

  virtual double DDensityDT(double T, double p) override
  {
    double val = table_->DFunctionDT(T, p, &ierr_);
    if (ierr_ != 0) error_msg_ = table_->ErrorMessage(T, p);
    return val;
  }

  virtual double DDensityDp(double T, double p) override
  {
    double val = table_->DFunctionDp(T, p, &ierr_);
    if (ierr_ != 0) error_msg_ = table_->ErrorMessage(T, p);
    return val;
  }

  virtual double MolarDensity(double T, double p) override { return Density(T, p) / M_; }
  virtual double DMolarDensityDT(double T, double p) override { return DDensityDT(T, p) / M_; }
  virtual double DMolarDensityDp(double T, double p) override { return DDensityDp(T, p) / M_; }

  static Utils::RegisteredFactory<EOS_Density, DensityTabular> reg_;

 private:
  Teuchos::RCP<LookupTable> table_;
};

} // namespace AmanziEOS
} // namespace Amanzi

#endif
