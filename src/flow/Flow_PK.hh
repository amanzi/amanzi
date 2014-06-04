/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Neil Carlson (version 1) 
           Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)

  This is a base virtual class.
*/

#ifndef AMANZI_FLOW_PK_HH_
#define AMANZI_FLOW_PK_HH_

#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_FECrsMatrix.h"
#include "Teuchos_RCP.hpp"

#include "BDFFnBase.hh"
#include "CompositeVectorSpace.hh"
#include "TI_Specs.hh"
#include "VerboseObject.hh"

#include "checkpoint.hh"
#include "primary_variable_field_evaluator.hh"
#include "tensor.hh"

#include "FlowDefs.hh"
#include "FlowTypeDefs.hh"
#include "FlowBoundaryFunction.hh"
#include "FlowDomainFunction.hh"
#include "Flow_BC_Factory.hh"
#include "Flow_SourceFactory.hh"

namespace Amanzi {
namespace AmanziFlow {

double bestLSfit(const std::vector<double>& h, const std::vector<double>& error);

class Flow_PK : public Amanzi::BDFFnBase<CompositeVector> {
 public:
  Flow_PK();
  virtual ~Flow_PK() {};

  // main methods
  void Init();
  virtual void InitPK() = 0;
  virtual void InitPicard(double T0) = 0;
  virtual void InitSteadyState(double T0, double dT0) = 0;
  virtual void InitTransient(double T0, double dT0) = 0;

  virtual double CalculateFlowDt() = 0;
  virtual int Advance(double dT, double& dT_actual) = 0; 
  virtual int AdvanceToSteadyState(double T0, double dT0) = 0;
  virtual void InitializeAuxiliaryData() = 0;
  virtual void InitializeSteadySaturated() = 0;

  virtual void CommitState(Teuchos::RCP<State> S) = 0;

  void UpdateAuxilliaryData();  // auxilliary data management
  void InitializeFields();

  // boundary and source teerms
  void ProcessBCs();
  void ComputeBCs(const CompositeVector& pressure);

  void AddSourceTerms(CompositeVector& rhs);

  // absolute permeability
  void SetAbsolutePermeabilityTensor();
  void CalculatePermeabilityFactorInWell();

  void ProcessShiftWaterTableList(const Teuchos::ParameterList& list);
  void CalculateShiftWaterTable(const std::string region);

  // miscallenous members
  void ResetPKtimes(double T0, double dT0) { T_physics = T0; dT = dT0; }
  void DeriveFaceValuesFromCellValues(const Epetra_MultiVector& ucells, Epetra_MultiVector& ufaces);
  int FindPosition(int f, AmanziMesh::Entity_ID_List faces);

  // io members
  void ProcessParameterList(Teuchos::ParameterList& list);
  void ProcessSublistTimeIntegration(Teuchos::ParameterList& list, const std::string name, TI_Specs& ti_specs);
  void ProcessStringSourceDistribution(const std::string name, int* method);
  void ProcessStringMFD3D(const std::string name, int* method);
  void ProcessStringTimeIntegration(const std::string name, int* method);
  void ProcessStringLinearSolver(const std::string& name, LinearSolver_Specs* ls_specs);
  void ProcessStringPreconditioner(const std::string& name, int* preconditioner);
  void ProcessStringErrorOptions(Teuchos::ParameterList& list, int* control);

  std::string FindStringLinearSolver(const Teuchos::ParameterList& plist);
  std::string FindStringPreconditioner(const Teuchos::ParameterList& list);

  void OutputTimeHistory(const Teuchos::ParameterList& plist, std::vector<dt_tuple>& dT_history);
  void WriteGMVfile(Teuchos::RCP<State> S) const;

  // utilities
  double WaterVolumeChangePerSecond(std::vector<int>& bc_model, Epetra_MultiVector& darcy_flux);

  void CalculateDarcyVelocity(std::vector<AmanziGeometry::Point>& xyz, 
                              std::vector<AmanziGeometry::Point>& velocity);
  void CalculatePoreVelocity(std::vector<AmanziGeometry::Point>& xyz, 
                             std::vector<AmanziGeometry::Point>& velocity,
                             std::vector<double>& porosity, std::vector<double>& saturation);
  void WriteWalkabout(const Teuchos::Ptr<Checkpoint>& wlk);

  // V&V
  void VV_ValidateBCs() const;
  void VV_PrintHeadExtrema(const CompositeVector& pressure) const;

  // extensions 
  int BoundaryFaceGetCell(int f);  // of AmanziMesh

  void set_intersection(const std::vector<AmanziMesh::Entity_ID>& v1,  // of std
                        const std::vector<AmanziMesh::Entity_ID>& v2, 
                        std::vector<AmanziMesh::Entity_ID>* vv);

  // support of unit tests
  double rho() { return rho_; }
  double mu() { return mu_; }
  const AmanziGeometry::Point& gravity() { return gravity_; }
  const TI_Specs& ti_specs_sss() { return ti_specs_sss_; }

 public:
  int ncells_owned, ncells_wghost;
  int nfaces_owned, nfaces_wghost;

  double T_physics, dT, dTnext;

  int MyPID;  // parallel information: will be moved to private
  int missed_bc_faces_;
  int ti_phase_counter;

 public:
  Teuchos::ParameterList linear_operator_list_;
  Teuchos::ParameterList preconditioner_list_;

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  int dim;

  Teuchos::RCP<State> S_;
  std::string passwd_;

  // Stationary physical quantatities
  std::vector<WhetStone::Tensor> K; 
  AmanziGeometry::Point gravity_;
  double g_, rho_, mu_, atm_pressure_;

  Teuchos::RCP<Epetra_Vector> Kxy;

  // boundary conditons
  Functions::FlowBoundaryFunction* bc_pressure; 
  Functions::FlowBoundaryFunction* bc_head;
  Functions::FlowBoundaryFunction* bc_flux;
  Functions::FlowBoundaryFunction* bc_seepage;
  int nseepage_prev;

  std::vector<int> bc_model, bc_submodel; 
  std::vector<double> bc_value, bc_coef;

  std::vector<double> rainfall_factor;
  Teuchos::RCP<Epetra_Vector> shift_water_table_;

  // source and sink terms
  Functions::FlowDomainFunction* src_sink;
  int src_sink_distribution; 

  // discretization and solvers
  int mfd3d_method_, mfd3d_method_preconditioner_;

  // time integration phases
  TI_Specs ti_specs_igs_;
  TI_Specs ti_specs_sss_;
  TI_Specs ti_specs_trs_;
  TI_Specs* ti_specs;

  // field evaluators (MUST GO AWAY lipnikov@lanl.gov)
  Teuchos::RCP<PrimaryVariableFieldEvaluator> darcy_flux_eval;

 protected:
  VerboseObject* vo_;
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
