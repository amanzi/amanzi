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
#include "PK.hh"
#include "PK_Explicit.hh"
#include "PK_Factory.hh"
#include "PK_Physical.hh"
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

  std::vector<double> NumFlux_x(std::vector<double>,std::vector<double>);

  std::vector<double> NumFlux_x_Rus(std::vector<double>,std::vector<double>);

  std::vector<double> NumFlux_x_central_upwind(std::vector<double>,std::vector<double>);

  std::vector<double> PhysSrc(std::vector<double>);

  std::vector<double> NumSrc(std::vector<double>,int);

  void BJ_lim(const WhetStone::DenseMatrix&,WhetStone::DenseMatrix&,int,const Key&);

  void ComputeGradients(Key,Key,Key);

  double Reconstruction(double,double,int,const Key&,const Key&,const Key&);
  double Reconstruction(double,double,int,const Key&);

 protected:
  Teuchos::RCP<Teuchos::ParameterList> glist_;
  Teuchos::RCP<Teuchos::ParameterList> sw_list_;
  Teuchos::RCP<TreeVector> soln_;
  Teuchos::RCP<State> S_;

  double dummy_dt;
  int step_count;

  Key domain_;

  // names of state fields
  Key velocity_x_key_;
  Key velocity_y_key_;
  Key discharge_x_key_;
  Key discharge_y_key_;
  Key ponded_depth_key_;
  Key total_depth_key_;
  Key bathymetry_key_;
  Key velocity_x_dx_key_;
  Key velocity_x_dy_key_;
  Key velocity_y_dx_key_;
  Key velocity_y_dy_key_;
  Key discharge_x_dx_key_;
  Key discharge_x_dy_key_;
  Key discharge_y_dx_key_;
  Key discharge_y_dy_key_;
  Key ponded_depth_dx_key_;
  Key ponded_depth_dy_key_;
  Key total_depth_dx_key_;
  Key total_depth_dy_key_;
  Key bathymetry_dx_key_;
  Key bathymetry_dy_key_;
  Key myPID_;

  std::string passwd_;

  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  int dim_;

 private:
  // boundary conditions
  std::vector<Teuchos::RCP<ShallowWaterBoundaryFunction> > bcs_; 
  std::vector<Teuchos::RCP<Operators::BCs> > op_bcs_;

 private:
  // factory registration
  static RegisteredPKFactory<ShallowWater_PK> reg_;
};
    
}  // namespace ShallowWater
}  // namespace Amanzi

#endif

