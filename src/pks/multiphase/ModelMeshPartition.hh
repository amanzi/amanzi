/*
  Multiphase PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  A collection of models along with a mesh partition.
*/

#ifndef AMANZI_MULTIPHASE_MODEL_MESH_PARTITION_HH_
#define AMANZI_MULTIPHASE_MODEL_MESH_PARTITION_HH_

#include <utility>

#include "dbc.hh"
#include "Factory.hh"
#include "Mesh.hh"
#include "MeshPartition.hh"

namespace Amanzi {
namespace Multiphase {

template <class Model>
using ModelPartition =
  std::pair<Teuchos::RCP<Functions::MeshPartition>, std::vector<Teuchos::RCP<Model>>>;

/* ******************************************************************
* Non-member factory.
****************************************************************** */
template <class Model>
Teuchos::RCP<ModelPartition<Model>>
CreateModelPartition(Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                     Teuchos::RCP<Teuchos::ParameterList> plist,
                     const std::string& key)
{
  Utils::Factory<Model> factory;
  std::vector<Teuchos::RCP<Model>> model_list;
  std::vector<std::vector<std::string>> region_list;

  for (auto lcv = plist->begin(); lcv != plist->end(); ++lcv) {
    std::string name = lcv->first;
    if (plist->isSublist(name)) {
      Teuchos::ParameterList sublist = plist->sublist(name);
      region_list.push_back(sublist.get<Teuchos::Array<std::string>>("regions").toVector());

      std::string model_name = sublist.get<std::string>(key);
      model_list.push_back(Teuchos::rcp(factory.CreateInstance(model_name, sublist)));
    } else {
      AMANZI_ASSERT(false);
    }
  }

  auto partition = Teuchos::rcp(new Functions::MeshPartition());
  partition->Initialize(mesh, AmanziMesh::CELL, region_list, -1);
  partition->Verify();

  return Teuchos::rcp(new ModelPartition<Model>(partition, model_list));
}

} // namespace Multiphase
} // namespace Amanzi

#endif
