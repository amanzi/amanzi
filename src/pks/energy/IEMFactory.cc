/*
  Energy

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Self-registering factory for IEM implementations.
*/

#include <string>
#include "IEMFactory.hh"

namespace Amanzi {
namespace Energy {

/* ******************************************************************
* method for instantiating IEM implementations
****************************************************************** */
Teuchos::RCP<IEM>
IEMFactory::CreateIEM(Teuchos::ParameterList& plist)
{
  std::string iem_typename = plist.get<std::string>("iem type");
  return Teuchos::rcp(CreateInstance(iem_typename, plist));
};

} // namespace Energy
} // namespace Amanzi
