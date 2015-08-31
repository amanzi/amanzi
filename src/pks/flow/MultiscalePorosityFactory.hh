/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov

  Self-registering factory for mulscale porosity models.
*/

#ifndef PK_FLOW_MULTISCALE_POROSITY_FACTORY_HH_
#define PK_FLOW_MULTISCALE_POROSITY_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "MultiscalePorosity.hh"
#include "factory.hh"

namespace Amanzi {
namespace Flow {

class MultiscalePorosityFactory : public Utils::Factory<MultiscalePorosity> {
 public:
  Teuchos::RCP<MultiscalePorosity> Create(Teuchos::ParameterList& plist);
};

}  // namespace Flow
}  // namespace Amanzi

#endif
