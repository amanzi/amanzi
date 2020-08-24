/*
 Shallow Water PK
 
 Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
 Amanzi is released under the three-clause BSD License.
 The terms of use and "as is" disclaimer for this license are
 provided in the top-level COPYRIGHT file.
 
 Author: Svetlana Tokareva (tokareva@lanl.gov)
 */

#ifndef AMANZI_SHALLOW_WATER_PK_HH_
#define AMANZI_SHALLOW_WATER_PK_HH_

// TPLs
#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_Import.h"
#include "Teuchos_RCP.hpp"

// Amanzi
#include "BCs.hh"
#include "CompositeVector.hh"
#include "DenseMatrix.hh"
#include "DenseVector.hh"
#include "Key.hh"
#include "LimiterCell.hh"
#include "PK.hh"
#include "PK_Explicit.hh"
#include "PK_Factory.hh"
#include "PK_Physical.hh"
#include "ReconstructionCell.hh"
#include "State.hh"
#include "Tensor.hh"
#include "Units.hh"
#include "VerboseObject.hh"
#include "WhetStoneMeshUtils.hh"

#include "ShallowWaterBoundaryFunction.hh"

namespace Amanzi {
namespace ShallowWater {
    
class ShallowWater_PK : public PK_Physical,
                        public PK_Explicit<Epetra_Vector> {
 public:
  ShallowWater_PK(Teuchos::ParameterList& pk_tree,
                  const Teuchos::RCP<Teuchos::ParameterList>& glist,
                  const Teuchos::RCP<State>& S,
                  const Teuchos::RCP<TreeVector>& soln);

  ShallowWater_PK(const Teuchos::RCP<Teuchos::ParameterList>& glist,
                  Teuchos::RCP<State> S,
                  const std::string& pk_list_name,
                  std::vector<std::string>& component_names);

  ~ShallowWater_PK() {};

  virtual void Setup(const Teuchos::Ptr<State>& S) override;
  virtual void Initialize(const Teuchos::Ptr<State>& S) override;

  virtual double get_dt() override;
  virtual void set_dt(double dt) override {};

  // Advance PK by step size dt.
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit=false) override;

  virtual void FunctionalTimeDerivative(double t, const Epetra_Vector& component,
                                        Epetra_Vector& f_component) override;

  // Commit any secondary (dependent) variables.
  virtual void CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S) override {};

  // Calculate any diagnostics prior to doing vis
  virtual void CalculateDiagnostics(const Teuchos::RCP<State>& S) override {};

  virtual std::string name() override { return "Shallow water PK"; }

  std::vector<double> PhysFlux_x(std::vector<double>);

  std::vector<double> PhysFlux_y(std::vector<double>);

  std::vector<double> NumFlux_x(std::vector<double>&, std::vector<double>&);

  std::vector<double> NumFlux_x_Rus(std::vector<double>&, std::vector<double>&);

  std::vector<double> NumFlux_x_central_upwind(std::vector<double>&, std::vector<double>&);

  std::vector<double> PhysSrc(std::vector<double>);

  std::vector<double> NumSrc(std::vector<double>,int);

 protected:
  Teuchos::RCP<Teuchos::ParameterList> glist_;
  Teuchos::RCP<Teuchos::ParameterList> sw_list_;
  Teuchos::RCP<TreeVector> soln_;
  Teuchos::RCP<State> S_;

  Key domain_;

  // names of state fields
  Key velocity_x_key_, velocity_y_key_;
  Key discharge_x_key_, discharge_y_key_;
  Key ponded_depth_key_;
  Key total_depth_key_;
  Key bathymetry_key_;
  Key discharge_y_grad_key_;

  std::string passwd_;

  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  int dim_;

 private:
  // boundary conditions
  std::vector<Teuchos::RCP<ShallowWaterBoundaryFunction> > bcs_; 
  std::vector<Teuchos::RCP<Operators::BCs> > op_bcs_;

  // gravity magnitude
  double g_;

  // limited reconstruction
  Teuchos::RCP<Operators::ReconstructionCell> total_depth_grad_, bathymetry_grad_;
  Teuchos::RCP<Operators::ReconstructionCell> velocity_x_grad_, velocity_y_grad_;
  Teuchos::RCP<Operators::ReconstructionCell> discharge_x_grad_, discharge_y_grad_;
  Teuchos::RCP<Operators::LimiterCell> limiter_;

 private:
  // factory registration
  static RegisteredPKFactory<ShallowWater_PK> reg_;
};
    
}  // namespace ShallowWater
}  // namespace Amanzi

#endif

