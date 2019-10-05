/*
This is the multiphase component of the Amanzi code. 

Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Authors: Quan Bui (mquanbui@math.umd.edu)
*/

#ifndef AMANZI_MULTIPHASE_RELATIVE_PERMEABILITY_HH_
#define AMANZI_MULTIPHASE_RELATIVE_PERMEABILITY_HH_

#include <string>
#include <vector>

#include "Epetra_IntVector.h"
#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Mesh.hh"
#include "CompositeVector.hh"
#include "ParallelCommunication.hh"
#include "Tensor.hh"
#include "VerboseObject.hh"

#include "WaterRetentionModel.hh"
#include "MultiphaseTypeDefs.hh"

namespace Amanzi {
namespace Multiphase {

class RelativePermeability {
 public:
  RelativePermeability(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) : mesh_(mesh) {};
  ~RelativePermeability() { delete vo_; }

  // main methods
  void Init(std::string phase_name, Teuchos::ParameterList& plist);
  void Compute(const CompositeVector& saturation_water);

  void DerivedSdP(const Epetra_MultiVector& p, Epetra_MultiVector& ds);
  void DerivedKdP(const Epetra_MultiVector& p, Epetra_MultiVector& dk);

  double Value(int c, double Sw) const { return WRM_[(*map_c2mb_)[c]]->k_relative(Sw, phase_); } 
  // hard-coded S2 as primary variable, multiply by -1 for derivative
  double Derivative(int c, double Sw) const { return - WRM_[(*map_c2mb_)[c]]->dKdS(Sw, phase_); } 

  double Value(int c, double x, double Sw) const { 
    return x * WRM_[(*map_c2mb_)[c]]->k_relative(Sw, phase_);
  }
  double Derivative(int c, double x, double Sw) const {
    return x * WRM_[(*map_c2mb_)[c]]->dKdS(Sw, phase_);
  }

  void ComputeGravityFlux(const std::vector<WhetStone::Tensor>& K, const AmanziGeometry::Point& g,
                          Teuchos::RCP<CompositeVector> flux);

  void VerifyWRMparameters(double m, double alpha, double sr, double pc0) {}
  void VerifyStringMualemBurdine(const std::string name) {}
  void ProcessStringRelativePermeability(const std::string name);
  void PlotWRMcurves(Teuchos::ParameterList& plist) {}

  // access methods
  std::vector<Teuchos::RCP<WaterRetentionModel> >& WRM() { return WRM_; }
  
  Teuchos::RCP<const AmanziMesh::Mesh> Mesh() { return mesh_; }
  Teuchos::RCP<CompositeVector> dKdS() { return dKdS_; }
  Teuchos::RCP<CompositeVector> Krel() { return Krel_; }

  //int method() { return method_; }
  const Epetra_IntVector& map_c2mb() { return *map_c2mb_; }

 private:
  void ProcessParameterList_(Teuchos::ParameterList& plist);
  void PopulateMapC2MB_();

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<ParallelCommunication> pp_;

  std::vector<Teuchos::RCP<WaterRetentionModel> > WRM_;
  std::string phase_;
  //double atm_pressure;

  //int method_;  // method for calculating relative permeability
  Teuchos::RCP<CompositeVector> Krel_;  // realitive permeability 
  Teuchos::RCP<CompositeVector> dKdS_;  // derivative of realitive permeability wrt saturation

  Teuchos::RCP<Epetra_IntVector> map_c2mb_;  // cell->model map

 protected:
  VerboseObject* vo_;
};

typedef double(RelativePermeability::*RelativePermeabilityUpwindFn)(int c, double s) const; 
typedef double(RelativePermeability::*FracKrelUpwindFn)(int c, double x, double s) const; 

}  // namespace Flow
}  // namespace Amanzi

#endif

