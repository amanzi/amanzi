/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Self-registering factory for Viscosity models.
*/

#ifndef AMANZI_EOS_VISCOSITY_BASE_FACTORY_HH_
#define AMANZI_EOS_VISCOSITY_BASE_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "ViscosityBase.hh"
#include "Factory.hh"

namespace Amanzi {
namespace AmanziEOS {

class ViscosityBaseFactory : public Utils::Factory<ViscosityBase> {
 public:
  Teuchos::RCP<ViscosityBase> CreateViscosity(Teuchos::ParameterList& plist);
};

}  // namespace AmanziEOS
}  // namespace Amanzi

#endif
