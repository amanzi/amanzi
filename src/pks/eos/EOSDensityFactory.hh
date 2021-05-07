/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Self-registering factory for EOS implementations.
*/

#ifndef AMANZI_PK_EOS_FACTORY_HH_
#define AMANZI_PK_EOS_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "EOS_Density.hh"
#include "Factory.hh"

namespace Amanzi {
namespace AmanziEOS {

class EOSDensityFactory : public Utils::Factory<EOS_Density> {
 public:
  Teuchos::RCP<EOS_Density> CreateEOS(Teuchos::ParameterList& plist);
};

}  // namespace AmanziEOS
}  // namespace Amanzi

#endif
