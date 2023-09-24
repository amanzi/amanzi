/*
 Shallow Water PK

 Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
 Amanzi is released under the three-clause BSD License.
 The terms of use and "as is" disclaimer for this license are
 provided in the top-level COPYRIGHT file.

 Authors: Svetlana Tokareva (tokareva@lanl.gov)
          Giacomo Capodaglio (gcapodaglio@lanl.gov)
 */

#ifndef AMANZI_SHALLOW_WATER_PK_HH_
#define AMANZI_SHALLOW_WATER_PK_HH_

#include <memory>

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
#include "Explicit_TI_RK.hh"
#include "Key.hh"
#include "LimiterCell.hh"
#include "NumericalFlux.hh"
#include "PK.hh"
#include "PK_Explicit.hh"
#include "PK_Factory.hh"
#include "PK_Physical.hh"
#include "PK_DomainFunction.hh"
#include "ReconstructionCellLinear.hh"
#include "State.hh"
#include "Tensor.hh"
#include "Units.hh"
#include "VerboseObject.hh"
#include "WhetStoneMeshUtils.hh"

#include "ShallowWaterBoundaryFunction.hh"

namespace Amanzi {
namespace ShallowWater {

class ShallowWater_PK : public PK_Physical, public PK_Explicit<TreeVector> {
 public:
  ShallowWater_PK(Teuchos::ParameterList& pk_tree,
                  const Teuchos::RCP<Teuchos::ParameterList>& glist,
                  const Teuchos::RCP<State>& S,
                  const Teuchos::RCP<TreeVector>& soln);
  ~ShallowWater_PK(){};

  virtual void Setup() override;
  virtual void Initialize() override;

  virtual double get_dt() override;
  virtual void set_dt(double dt) override{};

  // Advance PK by step size dt.
  virtual bool
  AdvanceStep(double t_old, double t_new, bool reinit = false) override;

  virtual void FunctionalTimeDerivative(double t, const TreeVector& A, TreeVector& f) override;

  virtual void ModifySolution(double t, TreeVector& A) override { VerifySolution_(A); }

  virtual void SetupPrimaryVariableKeys();

  virtual void SetupExtraEvaluatorsKeys(){};

  virtual void ScatterMasterToGhostedExtraEvaluators(){};

  virtual void UpdateExtraEvaluators(){};

  virtual void SetPrimaryVariableBC(Teuchos::RCP<Teuchos::ParameterList> &bc_list);

  virtual void InitializeFields();

  virtual void ComputeCellArrays();

  virtual void ComputeExternalForcingOnCells(std::vector<double> &forcing); 

  // Commit any secondary (dependent) variables.
  virtual void CommitStep(double t_old, double t_new, const Tag& tag) override;

  // Calculate any diagnostics prior to doing vis
  virtual void CalculateDiagnostics(const Tag& tag) override{};

  virtual std::string name() override { return "shallow water"; }

  // Bathymetry reconstruction on cell edge midpoints
  double BathymetryEdgeValue(int e, const Epetra_MultiVector& Bn);

  // Recalculate total depth for positivity of ponded depth
  double TotalDepthEdgeValue(int c, int e, double htc,
                             double Bc, double Bmax, const Epetra_MultiVector& B_n);

  virtual std::vector<double> ComputeFieldsOnEdge(int c, int e, double htc, 
                                                  double Bc, double Bmax, const Epetra_MultiVector& B_n); 


  // due to rotational invariance of SW equations, we need flux in the
  // x-direction only.
  std::vector<double> PhysicalFlux_x(const std::vector<double>&);

  std::vector<double>
  NumericalFlux_x(std::vector<double>&, std::vector<double>&);
  std::vector<double> NumericalFlux_x_Rusanov(const std::vector<double>&, const std::vector<double>&);
  std::vector<double> NumericalFlux_x_CentralUpwind(const std::vector<double>&, const std::vector<double>&);

  std::vector<double> NumericalSourceBedSlope(int c, double hc);

  virtual std::vector<double> NumericalSourceBedSlope(int c, double htc, double Bc,
                                                      double Bmax, const Epetra_MultiVector& B_n,
                                                      std::vector<int> bc_model, std::vector<double> bc_value_h); 

  virtual double NumericalSourceFriction(double h, double qx, double WettedAngle){return 0.0;};

  virtual double ComputeWaterDepth(double WettedAngle){return 0.0;};

  virtual double ComputeWettedAngle(double WaterDepth){return -1.0;};

  virtual double ComputeWettedAngleNewton(double WettedArea){return -1.0;};

  virtual double ComputeWettedArea(double WettedAngle){return 0.0;};

  virtual double ComputeTotalDepth(double PrimaryVar, double Bathymetry, double WettedAngle){return PrimaryVar + Bathymetry;};

  virtual double ComputePressureHead(double WettedArea){return 0.0;};

  virtual void UpdateSecondaryFields();

  virtual void ProjectNormalOntoMeshDirection(int c, AmanziGeometry::Point &normal) {};

  virtual double ComputeHydrostaticPressureForce (std::vector<double> SolArray){return g_ * 0.5 * SolArray[0] * SolArray[0];};

  void PushBackBC(Teuchos::RCP<ShallowWaterBoundaryFunction> bc){bcs_.push_back(bc);};

  double inverse_with_tolerance(double h, double tol);

  // access
  double get_total_source() const { return total_source_; }

 private:
  void
  InitializeFieldFromField_(const std::string& field0,
                            const std::string& field1, bool call_evaluator);

  void VerifySolution_(TreeVector& A);
  int ErrorDiagnostics_(double t, int c, double h);

 protected:
  Teuchos::RCP<Teuchos::ParameterList> glist_;
  Teuchos::RCP<Teuchos::ParameterList> sw_list_;
  Teuchos::RCP<TreeVector> soln_;
  Teuchos::RCP<State> S_;

  Key domain_;

  // numerical flux
  std::shared_ptr<NumericalFlux> numerical_flux_;

  // names of state fields
  Key velocity_key_, discharge_key_;
  // the primary variable is:
  // ponded depth for shallow water
  // wetted area for pipe flow
  Key primary_variable_key_, prev_primary_variable_key_; 
  Key total_depth_key_, bathymetry_key_;
  Key hydrostatic_pressure_key_;
  Key riemann_flux_key_;
  Key wetted_angle_key_;
  Key source_key_;

  std::string passwd_;

  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  int dim_;
  
  // source terms
  std::vector<Teuchos::RCP<PK_DomainFunction>> srcs_;
  double total_source_;

  // gravity magnitude
  double g_;

  int shallow_water_model_; 

  double pipe_diameter_;

  double celerity_;

  double velocity_desingularization_eps_;

  double Pi = 3.14159265359;
  double TwoPi = 6.28318530718;

  std::vector<int> model_cells_owned_;
  std::vector<int> junction_cells_owned_;
  std::vector<int> model_cells_wghost_;
  std::vector<int> junction_cells_wghost_;
  bool cellArraysInitDone_ = false;

 private:
  // boundary conditions
  std::vector<Teuchos::RCP<ShallowWaterBoundaryFunction>> bcs_;
  std::vector<Teuchos::RCP<Operators::BCs>> op_bcs_;

  // limited reconstruction
  bool use_limiter_;
  Teuchos::RCP<Operators::ReconstructionCellLinear> total_depth_grad_, bathymetry_grad_;
  Teuchos::RCP<Operators::ReconstructionCellLinear> discharge_x_grad_, discharge_y_grad_;
  Teuchos::RCP<Operators::LimiterCell> limiter_;

  // advanced cfl control
  double cfl_, cfl_positivity_;
  int iters_, max_iters_;

  // control of PK
  int temporal_disc_order_;
  double cell_area2_max_;

 private:
  // factory registration
  static RegisteredPKFactory<ShallowWater_PK> reg_;
};

} // namespace ShallowWater
} // namespace Amanzi

#endif
