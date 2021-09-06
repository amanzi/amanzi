/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  EOS for liquid water for temperature between 0.5 and 800 C and 
  pressure between 634 Pa and 110 MPa
*/

#ifndef AMANZI_EOS_H2O_DENSITY_TABULAR_HH_
#define AMANZI_EOS_H2O_DENSITY_TABULAR_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "EOS_Density.hh"
#include "LookupTable.hh"

namespace Amanzi {
namespace AmanziEOS {

// Equation of State model
class H2O_DensityTabular : public EOS_Density {
 public:
  explicit H2O_DensityTabular(Teuchos::ParameterList& eos_plist);

  virtual double Density(double T, double p) override { return table_->Function(T, p); }
  virtual double DDensityDT(double T, double p) override { return table_->DFunctionDT(T, p); }
  virtual double DDensityDp(double T, double p) override { return table_->DFunctionDp(T, p); }

  virtual double MolarDensity(double T, double p) override { return Density(T, p) / M_; }
  virtual double DMolarDensityDT(double T, double p) override { return DDensityDT(T, p) / M_; }
  virtual double DMolarDensityDp(double T, double p) override { return DDensityDp(T, p) / M_; }

  static Utils::RegisteredFactory<EOS_Density, H2O_DensityTabular> factory_;

 private:
  Teuchos::RCP<LookupTable> table_;
};

}  // namespace AmanziEOS
}  // namespace Amanzi

#endif
