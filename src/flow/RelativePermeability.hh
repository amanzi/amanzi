/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
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

#include "VerboseObject.hh"
#include "Mesh.hh"
#include "tensor.hh"

#include "State.hh"
#include "WaterRetentionModel.hh"
#include "FlowTypeDefs.hh"

namespace Amanzi {
namespace AmanziFlow {

class RelativePermeability {
 public:
  RelativePermeability(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                       Teuchos::ParameterList& list)
      : mesh_(mesh), list_(list) {};
  ~RelativePermeability() {};

  // main methods
  void Init(double p0, Teuchos::RCP<State> S);
  void ProcessParameterList();

  void Compute(const CompositeVector& pressure, 
               const std::vector<int>& bc_model, const std::vector<bc_tuple>& bc_values);
  void ComputeInCells(const CompositeVector& pressure);
  void ComputeOnFaces(const CompositeVector& pressure,
                      const std::vector<int>& bc_model, const std::vector<bc_tuple>& bc_values);
  void ComputeDerivativeOnFaces(
      const CompositeVector& pressure,
      const std::vector<int>& bc_model, const std::vector<bc_tuple>& bc_values);

  void DerivedSdP(const Epetra_MultiVector& p, Epetra_MultiVector& ds);
  void DerivedKdP(const Epetra_MultiVector& p, Epetra_MultiVector& dk);

  void CalculateKVectorUnit(const std::vector<WhetStone::Tensor>& K, const AmanziGeometry::Point& g);
  void SetFullySaturated();
  void PopulateMapC2MB();

  void VerifyWRMparameters(double m, double alpha, double sr, double pc0);
  void VerifyStringMualemBurdine(const std::string name);
  void ProcessStringRelativePermeability(const std::string name);
  void PlotWRMcurves();

  // access methods
  std::vector<Teuchos::RCP<WaterRetentionModel> >& WRM() { return WRM_; }

  CompositeVector& dKdP() { return *dKdP_; }
  CompositeVector& Krel() { return *Krel_; }
  std::vector<std::vector<double> >& Krel_amanzi() { return Krel_amanzi_; }
  std::vector<AmanziGeometry::Point >& Kgravity_unit() {return Kgravity_unit_;}

  int method() { return method_; }
  Epetra_Vector& map_c2mb() { return *map_c2mb_; }
  void set_experimental_solver(int solver) { experimental_solver_ = solver; }

 private:
  void FaceArithmeticMean_(const CompositeVector& pressure);
  void FaceUpwindGravityInit_();
  void FaceUpwindGravityInit_(const AmanziGeometry::Point& g);
  void FaceUpwindGravity_(
      const CompositeVector& pressure,
      const std::vector<int>& bc_model, const std::vector<bc_tuple>& bc_values);
  void FaceUpwindGravityInSoil_(
      const CompositeVector& pressure,
      const std::vector<int>& bc_model, const std::vector<bc_tuple>& bc_values);
  void FaceUpwindFlux_(
      const CompositeVector& pressure, const Epetra_MultiVector& flux,
      const std::vector<int>& bc_model, const std::vector<bc_tuple>& bc_values);

  void DerivativeFaceUpwindGravity_(
      const CompositeVector& pressure,
      const std::vector<int>& bc_model, const std::vector<bc_tuple>& bc_values);
  void DerivativeFaceUpwindFlux_(
      const CompositeVector& pressure, const Epetra_MultiVector& flux,
      const std::vector<int>& bc_model, const std::vector<bc_tuple>& bc_values);

 protected:
  VerboseObject* vo_;

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::ParameterList list_;

  std::vector<Teuchos::RCP<WaterRetentionModel> > WRM_;
  double atm_pressure;

  int ncells_owned, ncells_wghost;
  int nfaces_owned, nfaces_wghost;

  int method_;  // method for calculating relative permeability
  Teuchos::RCP<CompositeVector> Krel_;  // realitive permeability 
  Teuchos::RCP<CompositeVector> dKdP_;  // derivative of realitive permeability 

  Teuchos::RCP<Epetra_IntVector> upwind_cell;
  Teuchos::RCP<Epetra_IntVector> downwind_cell;
  Teuchos::RCP<Epetra_IntVector> face_flag;
  std::vector<std::vector<double> > Krel_amanzi_;

  std::vector<AmanziGeometry::Point> Kgravity_unit_;  // normalized vector Kg

  // Miscallenous maps
  Teuchos::RCP<Epetra_Vector> map_c2mb_;

  // obsolete, must go away (lipnikov@lanl.gov)
  int experimental_solver_; 
  Teuchos::RCP<State> S_;
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif

