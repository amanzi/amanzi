/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Mechanics PK

  A collection of SSMs along with a mesh partition.
*/

#ifndef AMANZI_MECHANICS_SSM_PARTITION_HH_
#define AMANZI_MECHANICS_SSM_PARTITION_HH_

#include "Mesh.hh"
#include "MeshPartition.hh"
#include "SSM_HardinDrnevich.hh"

namespace Amanzi {
namespace Mechanics {

typedef std::vector<Teuchos::RCP<SSM>> SSMList;
typedef std::pair<Teuchos::RCP<Functions::MeshPartition>, SSMList> SSMPartition;

// Non-member factory
inline
Teuchos::RCP<SSMPartition>
CreateSSMPartition(Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                   Teuchos::RCP<Teuchos::ParameterList> plist)
{
  std::vector<Teuchos::RCP<SSM>> ssm_list;
  std::vector<std::vector<std::string>> region_list;

  for (auto lcv = plist->begin(); lcv != plist->end(); ++lcv) {
    std::string name = lcv->first;
    if (plist->isSublist(name)) {
      Teuchos::ParameterList sublist = plist->sublist(name);
      region_list.push_back(sublist.get<Teuchos::Array<std::string>>("regions").toVector());

      std::string model = sublist.get<std::string>("model");
      if (model == "Hardin Drnevich") {
        ssm_list.push_back(Teuchos::rcp(new SSM_HardinDrnevich(sublist)));
      }
    } else {
      AMANZI_ASSERT(0);
    }
  }

  auto partition = Teuchos::rcp(new Functions::MeshPartition());
  partition->Initialize(mesh, AmanziMesh::Entity_kind::CELL, region_list, -1);
  partition->Verify();

  return Teuchos::rcp(new SSMPartition(partition, ssm_list));
}

} // namespace Mechanics
} // namespace Amanzi

#endif
