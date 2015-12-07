/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  A collection of WRMs along with a Mesh Partition.
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOW_RELATIONS_WRM_PARTITION_
#define AMANZI_FLOW_RELATIONS_WRM_PARTITION_

#include "wrm.hh"
#include "wrm_permafrost_model.hh"
#include "MeshPartition.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

typedef std::vector<Teuchos::RCP<WRM> > WRMList;
typedef std::pair<Teuchos::RCP<Functions::MeshPartition>, WRMList> WRMPartition;

typedef std::vector<Teuchos::RCP<WRMPermafrostModel> > WRMPermafrostModelList;
typedef std::pair<Teuchos::RCP<Functions::MeshPartition>, WRMPermafrostModelList> WRMPermafrostModelPartition;

// Non-member factory
Teuchos::RCP<WRMPartition>
createWRMPartition(Teuchos::ParameterList& plist);

Teuchos::RCP<WRMPermafrostModelPartition>
createWRMPermafrostModelPartition(Teuchos::ParameterList& plist,
        Teuchos::RCP<WRMPartition>& wrms);

} // namespace
} // namespace
} // namespace

#endif
