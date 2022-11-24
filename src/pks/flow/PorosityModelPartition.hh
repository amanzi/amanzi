/*
  Flow PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)

  A collection of porosity models along with a mesh partition.
*/

#ifndef AMANZI_FLOW_POROSITY_MODEL_PARTITION_HH_
#define AMANZI_FLOW_POROSITY_MODEL_PARTITION_HH_

#include "Mesh.hh"
#include "MeshPartition.hh"
#include "PorosityModel.hh"

namespace Amanzi {
namespace Flow {

typedef std::vector<Teuchos::RCP<PorosityModel>> PorosityModelList;
typedef std::pair<Teuchos::RCP<Functions::MeshPartition>, PorosityModelList> PorosityModelPartition;

// Non-member factory
Teuchos::RCP<PorosityModelPartition>
CreatePorosityModelPartition(Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                             Teuchos::RCP<Teuchos::ParameterList> plist);

} // namespace Flow
} // namespace Amanzi

#endif
