/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Flow PK

  A collection of WRMs along with a Mesh Partition.
*/

#include "dbc.hh"
#include "WRMFactory.hh"
#include "WRMPartition.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* Non-member factory.
****************************************************************** */
Teuchos::RCP<WRMPartition>
CreateWRMPartition(Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                   Teuchos::RCP<Teuchos::ParameterList> plist)
{
  WRMFactory factory;
  std::vector<Teuchos::RCP<WRM>> wrm_list;
  std::vector<std::vector<std::string>> region_list;

  for (auto lcv = plist->begin(); lcv != plist->end(); ++lcv) {
    std::string name = lcv->first;
    if (plist->isSublist(name)) {
      Teuchos::ParameterList sublist = plist->sublist(name);
      region_list.push_back(sublist.get<Teuchos::Array<std::string>>("regions").toVector());
      wrm_list.push_back(factory.Create(sublist));
    } else {
      AMANZI_ASSERT(0);
    }
  }

  auto partition = Teuchos::rcp(new Functions::MeshPartition());
  partition->Initialize(mesh, AmanziMesh::CELL, region_list, -1);
  partition->Verify();

  return Teuchos::rcp(new WRMPartition(partition, wrm_list));
}

} // namespace Flow
} // namespace Amanzi
