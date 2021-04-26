/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  EOS for water viscosity for temperature between 0.5 and 800 C and 
  pressure between 634 Pa and 110 MPa
*/

#ifndef AMANZI_VISCOSITY_WATER_TABULAR_HH_
#define AMANZI_VISCOSITY_WATER_TABULAR_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "ViscosityBase.hh"
#include "LookupTable.hh"

namespace Amanzi {
namespace AmanziEOS {

class ViscosityWaterTabular : public ViscosityBase {
 public:
  explicit ViscosityWaterTabular(Teuchos::ParameterList& eos_plist);

  virtual double Viscosity(double T, double p) override { return table_->Function(T, p); }
  virtual double DViscosityDT(double T, double p) override { return table_->DFunctionDT(T, p); }
  virtual double DViscosityDp(double T, double p) override { return table_->DFunctionDp(T, p); }

  static Utils::RegisteredFactory<ViscosityBase, ViscosityWaterTabular> factory_;

 private:
  Teuchos::RCP<LookupTable> table_;
};

}  // namespace AmanziEOS
}  // namespace Amanzi

#endif
