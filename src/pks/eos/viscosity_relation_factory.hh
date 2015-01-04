/*
  This is the EOS component of the ATS and Amanzi codes.
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Self-registering factory for Viscosity implementations.
*/

#ifndef AMANZI_EOS_VISCOSITY_RELATION_FACTORY_HH_
#define AMANZI_EOS_VISCOSITY_RELATION_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "viscosity_relation.hh"
#include "factory.hh"

namespace Amanzi {
namespace Relations {

class ViscosityRelationFactory : public Utils::Factory<ViscosityRelation> {
 public:
  Teuchos::RCP<ViscosityRelation> createViscosity(Teuchos::ParameterList& plist);
};

}  // namespace Relations
}  // namespace Amanzi

#endif
