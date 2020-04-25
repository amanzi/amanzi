/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  Constant prescribed viscosity.
*/

#ifndef AMANZI_EOS_VISCOSITY_CONSTANT_HH_
#define AMANZI_EOS_VISCOSITY_CONSTANT_HH_

#include "Teuchos_ParameterList.hpp"

#include "dbc.hh"
#include "Factory.hh"

#include "Viscosity_Base.hh"

namespace Amanzi {
namespace AmanziEOS {

class Viscosity_Constant : public Viscosity_Base {
 public:
  explicit
  Viscosity_Constant(Teuchos::ParameterList& visc_plist);

  virtual double Viscosity(double T) { return visc_; }
  virtual double DViscosityDT(double T) { return 0.; }

 protected:
  virtual void InitializeFromPlist_();

  Teuchos::ParameterList visc_plist_;
  double visc_;

 private:
  static Utils::RegisteredFactory<Viscosity_Base, Viscosity_Constant> factory_;
};

}  // namespace AmanziEOS
}  // namespace Amanzi

#endif
