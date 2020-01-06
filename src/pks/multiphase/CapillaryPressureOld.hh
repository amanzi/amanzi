/*
This is the multiphase component of the Amanzi code. 

Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Authors: Quan Bui (mquanbui@math.umd.edu)
*/

#ifndef AMANZI_MULTIPHASE_CAPILLARY_PRESSURE_OLD_HH_
#define AMANZI_MULTIPHASE_CAPILLARY_PRESSURE_OLD_HH_

#include <string>
#include <vector>

#include "Epetra_IntVector.h"
#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Mesh.hh"
#include "CompositeVector.hh"
#include "ParallelCommunication.hh"
//#include "tensor.hh"
#include "VerboseObject.hh"

#include "WRMmp.hh"
#include "MultiphaseTypeDefs.hh"

namespace Amanzi {
namespace Multiphase {

class CapillaryPressureOld {
 public:
  CapillaryPressureOld(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) : mesh_(mesh) {};
  ~CapillaryPressureOld() {};

  // main methods
  void Init(Teuchos::ParameterList& plist);
  void Compute(const CompositeVector& saturation_water);

  // access methods
  Teuchos::RCP<WRMmpPartition>& wrm() { return wrm_; }
  
  Teuchos::RCP<const AmanziMesh::Mesh> Mesh() { return mesh_; }
  Teuchos::RCP<CompositeVector> dPc_dS() { return dPc_dS_; }
  Teuchos::RCP<CompositeVector> Pc() { return Pc_; }

 private:
  void ProcessParameterList_(Teuchos::ParameterList& plist);

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;

  Teuchos::RCP<WRMmpPartition> wrm_;

  Teuchos::RCP<CompositeVector> Pc_;  // realitive permeability 
  Teuchos::RCP<CompositeVector> dPc_dS_;  // derivative of realitive permeability wrt saturation
};

}  // namespace Flow
}  // namespace Amanzi

#endif

