/*
  Flow PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  A collection of multiscale porosity models along with a mesh partition.
*/

#ifndef MULTISCALE_FLOW_POROSITY_PARTITION_HH_
#define MULTISCALE_FLOW_POROSITY_PARTITION_HH_

#include "Mesh.hh"
#include "MeshPartition.hh"
#include "MultiscaleFlowPorosity.hh"

namespace Amanzi {
namespace Flow {

typedef std::vector<Teuchos::RCP<MultiscaleFlowPorosity>> MsFList;
typedef std::pair<Teuchos::RCP<Functions::MeshPartition>, MsFList> MultiscaleFlowPorosityPartition;

// Non-member factory
Teuchos::RCP<MultiscaleFlowPorosityPartition>
CreateMultiscaleFlowPorosityPartition(Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                                      Teuchos::RCP<Teuchos::ParameterList> plist);

// Non-member function quering partition
int
NumberMatrixNodes(Teuchos::RCP<MultiscaleFlowPorosityPartition>& msp);

} // namespace Flow
} // namespace Amanzi

#endif
