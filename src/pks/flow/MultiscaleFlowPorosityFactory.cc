/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
*/

/*
  Flow PK

  Self-registering factory for multiscale porosity models.
*/

#include <string>
#include "MultiscaleFlowPorosityFactory.hh"

namespace Amanzi {
namespace Flow {

// method for instantiating a multiscale porosity model
Teuchos::RCP<MultiscaleFlowPorosity>
MultiscaleFlowPorosityFactory::Create(Teuchos::ParameterList& plist)
{
  std::string msp_typename = plist.get<std::string>("multiscale model");
  return Teuchos::rcp(CreateInstance(msp_typename, plist));
};

} // namespace Flow
} // namespace Amanzi
