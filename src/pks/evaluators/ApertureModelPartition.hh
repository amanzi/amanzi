/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Common evaluators

  A collection of aperture models along with a mesh partition.
*/

#ifndef AMANZI_EVALUATORS_APERTURE_MODEL_PARTITION_HH_
#define AMANZI_EVALUATORS_APERTURE_MODEL_PARTITION_HH_

#include "Mesh.hh"
#include "MeshPartition.hh"
#include "ApertureModel.hh"

namespace Amanzi {
namespace Evaluators {

typedef std::vector<Teuchos::RCP<ApertureModel>> ApertureModelList;
typedef std::pair<Teuchos::RCP<Functions::MeshPartition>, ApertureModelList> ApertureModelPartition;

// Non-member factory
Teuchos::RCP<ApertureModelPartition> CreateApertureModelPartition(
  Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
  Teuchos::RCP<Teuchos::ParameterList> plist);

} // namespace Evaluators
} // namespace Amanzi

#endif
