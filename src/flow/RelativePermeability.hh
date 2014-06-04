/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_RELATIVE_PERMEABILITY_HH_
#define AMANZI_RELATIVE_PERMEABILITY_HH_

#include <string>
#include <vector>

#include "Epetra_IntVector.h"
#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Mesh.hh"
#include "CompositeVector.hh"
#include "ParallelCommunication.hh"
#include "tensor.hh"
#include "VerboseObject.hh"

#include "WaterRetentionModel.hh"
#include "FlowTypeDefs.hh"

namespace Amanzi {
namespace AmanziFlow {

class RelativePermeability {
 public:
  RelativePermeability(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) : mesh_(mesh) {};
  ~RelativePermeability() {};

  // main methods
  void Init(double p0, Teuchos::ParameterList& plist);
  void Compute(const CompositeVector& pressure);

  void DerivedSdP(const Epetra_MultiVector& p, Epetra_MultiVector& ds);
  void DerivedKdP(const Epetra_MultiVector& p, Epetra_MultiVector& dk);

  double Value(int c, double pc) const { return WRM_[(*map_c2mb_)[c]]->k_relative(atm_pressure - pc); } 

  void ComputeGravityFlux(const std::vector<WhetStone::Tensor>& K, const AmanziGeometry::Point& g,
                          Teuchos::RCP<CompositeVector> flux);

  void VerifyWRMparameters(double m, double alpha, double sr, double pc0);
  void VerifyStringMualemBurdine(const std::string name);
  void ProcessStringRelativePermeability(const std::string name);
  void PlotWRMcurves();

  // access methods
  std::vector<Teuchos::RCP<WaterRetentionModel> >& WRM() { return WRM_; }

  Teuchos::RCP<CompositeVector> dKdP() { return dKdP_; }
  Teuchos::RCP<CompositeVector> Krel() { return Krel_; }

  int method() { return method_; }
  const Epetra_IntVector& map_c2mb() { return *map_c2mb_; }

 private:
  void ProcessParameterList_(Teuchos::ParameterList& plist);
  void PopulateMapC2MB_();

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<ParallelCommunication> pp_;

  std::vector<Teuchos::RCP<WaterRetentionModel> > WRM_;
  double atm_pressure;

  int method_;  // method for calculating relative permeability
  Teuchos::RCP<CompositeVector> Krel_;  // realitive permeability 
  Teuchos::RCP<CompositeVector> dKdP_;  // derivative of realitive permeability 

  Teuchos::RCP<Epetra_IntVector> map_c2mb_;  // cell->model map

 protected:
  VerboseObject* vo_;
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif

