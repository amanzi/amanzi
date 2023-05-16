/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  EOS

  EOS for tabulated viscosity.
*/

#ifndef AMANZI_VISCOSITY_TABULAR_HH_
#define AMANZI_VISCOSITY_TABULAR_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "EOS_Viscosity.hh"
#include "Factory.hh"
#include "LookupTable.hh"

namespace Amanzi {
namespace AmanziEOS {

class ViscosityTabular : public EOS_Viscosity {
 public:
  explicit ViscosityTabular(Teuchos::ParameterList& eos_plist);

  virtual double Viscosity(double T, double p) override
  {
    double val = table_->Function(T, p, &ierr_);
    if (ierr_ != 0) error_msg_ = table_->ErrorMessage(T, p);
    return val;
  }

  virtual double DViscosityDT(double T, double p) override
  {
    double val = table_->DFunctionDT(T, p, &ierr_);
    if (ierr_ != 0) error_msg_ = table_->ErrorMessage(T, p);
    return val;
  }

  virtual double DViscosityDp(double T, double p) override
  {
    double val = table_->DFunctionDp(T, p, &ierr_);
    if (ierr_ != 0) error_msg_ = table_->ErrorMessage(T, p);
    return val;
  }

  static Utils::RegisteredFactory<EOS_Viscosity, ViscosityTabular> reg_;

 private:
  Teuchos::RCP<LookupTable> table_;
};

} // namespace AmanziEOS
} // namespace Amanzi

#endif
