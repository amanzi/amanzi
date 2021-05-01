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

#ifndef AMANZI_IEM_TABULAR_HH_
#define AMANZI_IEM_TABULAR_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "IEM.hh"
#include "Factory.hh"
#include "LookupTable.hh"

namespace Amanzi {
namespace Energy {

// Equation of State model
class IEM_Tabular : public IEM {
 public:
  IEM_Tabular(Teuchos::ParameterList& plist) {
    table_ = Teuchos::rcp(new AmanziEOS::LookupTable(plist));
  }

  virtual double InternalEnergy(double T, double p) override { return table_->Function(T, p); }
  virtual double DInternalEnergyDT(double T, double p) override { return table_->DFunctionDT(T, p); }
  virtual double DInternalEnergyDp(double T, double p) override { return table_->DFunctionDp(T, p); }

  static Utils::RegisteredFactory<IEM, IEM_Tabular> factory_;

 private:
  Teuchos::RCP<AmanziEOS::LookupTable> table_;
};

}  // namespace Energy
}  // namespace Amanzi

#endif
