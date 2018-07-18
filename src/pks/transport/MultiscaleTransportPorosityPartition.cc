/*
  Transport PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov

  A collection of multiscale porosity models along with a mesh partition.
*/

#include "dbc.hh"
#include "MultiscaleTransportPorosityFactory.hh"
#include "MultiscaleTransportPorosityPartition.hh"

namespace Amanzi {
namespace Transport {

/* ******************************************************************
* Non-member factory.
****************************************************************** */
Teuchos::RCP<MultiscaleTransportPorosityPartition> CreateMultiscaleTransportPorosityPartition(
    Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
    Teuchos::RCP<Teuchos::ParameterList> plist)
{
  MultiscaleTransportPorosityFactory factory;
  std::vector<Teuchos::RCP<MultiscaleTransportPorosity> > msp;
  std::vector<std::vector<std::string> > regions;

  for (auto it = plist->begin(); it != plist->end(); ++it) {
    if (plist->isSublist(it->first)) {
      Teuchos::ParameterList sublist = plist->sublist(it->first);
      regions.push_back(sublist.get<Teuchos::Array<std::string> >("regions").toVector());
      msp.push_back(factory.Create(sublist));
    } else {
      AMANZI_ASSERT(0);
    }
  }

  Teuchos::RCP<Functions::MeshPartition> partition = Teuchos::rcp(new Functions::MeshPartition());
  partition->Initialize(mesh, AmanziMesh::CELL, regions, -1);
  partition->Verify();

  return Teuchos::rcp(new MultiscaleTransportPorosityPartition(partition, msp));
}

}  // namespace Transport
}  // namespace Amanzi

