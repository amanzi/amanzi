/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)

  A collection of MultiscalePorosity models along with a mesh partition.
*/

#ifndef AMANZI_FLOW_MULTISCALE_POROSITY_PARTITION_HH_
#define AMANZI_FLOW_MULTISCALE_POROSITY_PARTITION_HH_

#include "Mesh.hh"
#include "MeshPartition.hh"
#include "MultiscalePorosity.hh"

namespace Amanzi {
namespace Flow {

typedef std::vector<Teuchos::RCP<MultiscalePorosity> > MsPList;
typedef std::pair<Teuchos::RCP<Functions::MeshPartition>, MsPList> MultiscalePorosityPartition;

// Non-member factory
Teuchos::RCP<MultiscalePorosityPartition> CreateMultiscalePorosityPartition(
    Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
    Teuchos::RCP<Teuchos::ParameterList> plist);

}  // namespace Flow
}  // namespace Amanzi

#endif
