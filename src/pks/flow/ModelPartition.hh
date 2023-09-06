/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Flow PK

  A collection of models along with a mesh partition.
*/

#ifndef AMANZI_FLOW_MODEL_PARTITION_HH_
#define AMANZI_FLOW_MODEL_PARTITION_HH_

#include "dbc.hh"
#include "Mesh.hh"
#include "MeshPartition.hh"

#include "ModelFactory.hh"

namespace Amanzi {
namespace Flow {

template <typename Model>
using ModelList = std::vector<Teuchos::RCP<Model>>;

template <typename Model>
using ModelPartition = std::pair<Teuchos::RCP<Functions::MeshPartition>, ModelList<Model>>;

// Non-member factory
template <typename Model>
Teuchos::RCP<ModelPartition<Model>>
CreateModelPartition(Teuchos::RCP<const AmanziMesh::Mesh>& mesh, Teuchos::ParameterList& plist)
{
  ModelFactory<Model> factory;
  ModelList<Model> model_list;
  std::vector<std::vector<std::string>> region_list;

  for (auto lcv = plist.begin(); lcv != plist.end(); ++lcv) {
    std::string name = lcv->first;
    if (plist.isSublist(name)) {
      Teuchos::ParameterList sublist = plist.sublist(name);
      region_list.push_back(sublist.get<Teuchos::Array<std::string>>("regions").toVector());
      auto model = factory.Create(sublist);
      if (!model.get()) {
        Errors::Message msg;
        msg << "Unknown model name in sublist \"" << name  << "\"";
        Exceptions::amanzi_throw(msg);
      }
      model_list.push_back(model);
    } else {
      AMANZI_ASSERT(0);
    }
  }

  auto partition = Teuchos::rcp(new Functions::MeshPartition());
  partition->Initialize(mesh, AmanziMesh::CELL, region_list, -1);
  partition->Verify();

  return Teuchos::rcp(new ModelPartition<Model>(partition, model_list));
}

} // namespace Flow
} // namespace Amanzi

#endif
