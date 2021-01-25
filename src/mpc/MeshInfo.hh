/*
  Multi-Process Coordinator

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Checkpointing mesh info.
*/

#ifndef AMANZI_MESHINFO_CHECKPOINT_HH_
#define AMANZI_MESHINFO_CHECKPOINT_HH_

#include <map>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"

#include "State.hh"
#include "Checkpoint.hh"

#include "VerboseObject.hh"


namespace Amanzi {

class MeshInfo : public Checkpoint {

public:
  MeshInfo(Teuchos::ParameterList& plist, const State& S): Checkpoint (plist, S) {};
  MeshInfo() : Checkpoint(true) {};

  void WriteMeshCentroids(std::string domain, const AmanziMesh::Mesh& mesh);
};

}  // namespace Amanzi

#endif
