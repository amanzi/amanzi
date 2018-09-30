/*
  Flow PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov

  A collection of multiscale porosity models along with a mesh partition.
*/

#include "dbc.hh"
#include "MultiscaleFlowPorosityFactory.hh"
#include "MultiscaleFlowPorosityPartition.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* Non-member factory.
****************************************************************** */
Teuchos::RCP<MultiscaleFlowPorosityPartition> CreateMultiscaleFlowPorosityPartition(
    Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
    Teuchos::RCP<Teuchos::ParameterList> plist)
{
  MultiscaleFlowPorosityFactory factory;
  std::vector<Teuchos::RCP<MultiscaleFlowPorosity> > msp_list;
  std::vector<std::vector<std::string> > region_list;

  for (auto lcv = plist->begin(); lcv != plist->end(); ++lcv) {
    std::string name = lcv->first;
    if (plist->isSublist(name)) {
      Teuchos::ParameterList sublist = plist->sublist(name);
      region_list.push_back(sublist.get<Teuchos::Array<std::string> >("regions").toVector());
      msp_list.push_back(factory.Create(sublist));
    } else {
      AMANZI_ASSERT(0);
    }
  }

  auto partition = Teuchos::rcp(new Functions::MeshPartition());
  partition->Initialize(mesh, AmanziMesh::CELL, region_list, -1);
  partition->Verify();

  return Teuchos::rcp(new MultiscaleFlowPorosityPartition(partition, msp_list));
}

}  // namespace Flow
}  // namespace Amanzi

