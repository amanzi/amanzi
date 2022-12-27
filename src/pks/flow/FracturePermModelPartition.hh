/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Flow PK

  A collection of fracture permeability models along with a mesh partition.
*/

#ifndef AMANZI_FLOW_FRACTURE_PERM_MODEL_PARTITION_HH_
#define AMANZI_FLOW_FRACTURE_PERM_MODEL_PARTITION_HH_

#include "Mesh.hh"
#include "MeshPartition.hh"
#include "FracturePermModel.hh"

namespace Amanzi {
namespace Flow {

typedef std::vector<Teuchos::RCP<FracturePermModel>> FracturePermModelList;
typedef std::pair<Teuchos::RCP<Functions::MeshPartition>, FracturePermModelList>
  FracturePermModelPartition;

// Non-member factory
Teuchos::RCP<FracturePermModelPartition>
CreateFracturePermModelPartition(Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                                 Teuchos::RCP<Teuchos::ParameterList> plist);

} // namespace Flow
} // namespace Amanzi

#endif
