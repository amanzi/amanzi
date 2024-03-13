/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Transport PK

  A collection of MDMs along with a mesh partition and
  supporting parameters.
*/

#ifndef AMANZI_TRANSPORT_MDM_PARTITION_HH_
#define AMANZI_TRANSPORT_MDM_PARTITION_HH_

#include "Mesh.hh"
#include "MeshPartition.hh"
#include "MDM.hh"

namespace Amanzi {
namespace Transport {

typedef std::vector<Teuchos::RCP<MDM>> MDMList;
typedef std::pair<Teuchos::RCP<Functions::MeshPartition>, MDMList> MDMPartition;

// Non-member factory returns false if all models generate zero dispersion tensor.
Teuchos::RCP<MDMPartition>
CreateMDMPartition(Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                   Teuchos::RCP<Teuchos::ParameterList> plist,
                   bool& flag);

} // namespace Transport
} // namespace Amanzi

#endif
