/*
This is the Flow component of the Amanzi code. 
License: BSD
Authors: Neil Carlson (version 1) 
         Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
Usage: 
  Flow_PK TPK(Teuchos::ParameterList& list, Teuchos::RCP<Flow_State> FS);
  double time_step = FPK.calculateFlowDt();
  FPK.advance(any_dT);
*/

#ifndef __FLOW_PK_HH__
#define __FLOW_PK_HH__

#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_FECrsMatrix.h"
#include "Teuchos_RCP.hpp"

#include "BoundaryFunction.hh"
#include "DomainFunction.hh"
#include "mfd3d.hh"
#include "BDF_FnBase.hh"

#include "Flow_State.hh"
#include "Flow_constants.hh"
#include "Flow_typedefs.hh"
#include "Matrix_MFD.hh"
#include "Matrix_MFD_PLambda.hh"
#include "TI_Specs.hh"


namespace Amanzi {
namespace AmanziFlow {

double bestLSfit(const std::vector<double>& h, const std::vector<double>& error);

class Flow_PK : public BDF2::fnBase {
 public:
  Flow_PK() {};
  virtual ~Flow_PK() {};

  // main methods
  void Init(Teuchos::RCP<Flow_State> FS_MPC);
  virtual void InitPK() = 0;
  virtual void InitSteadyState(double T0, double dT0) = 0;
  virtual void InitTransient(double T0, double dT0) = 0;
  virtual void InitPicard(double T0) = 0;

  virtual double CalculateFlowDt() = 0;
  virtual int Advance(double dT) = 0; 
  virtual int AdvanceToSteadyState(double T0, double dT0) = 0;
  virtual void InitializeAuxiliaryData() = 0;
  virtual void InitializeSteadySaturated() = 0;

  virtual void CommitState(Teuchos::RCP<Flow_State> FS) = 0;

  // boundary condition members
  void ProcessStaticBCsubmodels(const std::vector<int>& bc_submodel,
                                std::vector<double>& rainfall_factor);
  void ProcessBoundaryConditions(
      BoundaryFunction* bc_pressure, BoundaryFunction* bc_head,
      BoundaryFunction* bc_flux, BoundaryFunction* bc_seepage,
      const Epetra_Vector& pressure_cells, 
      const Epetra_Vector& pressure_faces, const double atm_pressure,
      const std::vector<double>& rainfall_factor,
      const std::vector<int>& bc_submodel,
      std::vector<int>& bc_model, std::vector<bc_tuple>& bc_values);

  void CalculatePermeabilityFactorInWell(const std::vector<WhetStone::Tensor>& K, Epetra_Vector& Kxy);
  void AddSourceTerms(DomainFunction* src_sink, Epetra_Vector& rhs);

  void ProcessShiftWaterTableList(
      const Teuchos::ParameterList& list, BoundaryFunction* bc_head,
      Teuchos::RCP<Epetra_Vector> shift_water_table_);
  void CalculateShiftWaterTable(
      const std::string region, Teuchos::RCP<Epetra_Vector> shift_water_table_);

  double WaterVolumeChangePerSecond(std::vector<int>& bc_model, Epetra_Vector& darcy_flux);

  // gravity members
  void AddGravityFluxes_MFD(
      std::vector<WhetStone::Tensor>& K,
      const Epetra_Vector& Krel_cells, const Epetra_Vector& Krel_faces, int method,
      Matrix_MFD* matrix);
  void AddGravityFluxes_DarcyFlux(
      std::vector<WhetStone::Tensor>& K,
      const Epetra_Vector& Krel_cells, const Epetra_Vector& Krel_faces, int method,
      Epetra_Vector& darcy_mass_flux);

  // Picard-Newton members
  void AddNewtonFluxes_MFD(const Epetra_Vector& dKdP_faces, const Epetra_Vector& Krel_faces,
                           const Epetra_Vector& pressure_faces, const Epetra_Vector& flux,
                           Epetra_Vector& rhs, Matrix_MFD_PLambda* matrix);

  // access members 
  Teuchos::RCP<Flow_State> flow_state() { return FS; }
  int flow_status() { return flow_status_; }

  // control members
  void ValidateBoundaryConditions(
      BoundaryFunction *bc_pressure, BoundaryFunction *bc_head, BoundaryFunction *bc_flux) const;
  void WriteGMVfile(Teuchos::RCP<Flow_State> FS) const;
 
  void ResetPKtimes(double T0, double dT0) { T_physics = T0; dT = dT0; }
  
  // miscallenous members
  Epetra_Map* CreateSuperMap();
  void DeriveFaceValuesFromCellValues(const Epetra_Vector& ucells, Epetra_Vector& ufaces);
  void IdentifyUpwindCells(Epetra_IntVector& upwind_cell, Epetra_IntVector& downwind_cell);
  int FindPosition(int f, AmanziMesh::Entity_ID_List faces);

  // io members
  void ProcessSublistTimeIntegration(Teuchos::ParameterList& list, const std::string name, TI_Specs& ti_specs);
  void ProcessStringSourceDistribution(const std::string name, int* method);
  void ProcessStringMFD3D(const std::string name, int* method);
  void ProcessStringTimeIntegration(const std::string name, int* method);
  void ProcessStringVerbosity(const std::string name, int* verbosity);
  void ProcessStringPreconditioner(const std::string name, int* preconditioner);
  std::string FindStringLinearSolver(const Teuchos::ParameterList& list, const Teuchos::ParameterList& solver_list);
  void OutputTimeHistory(std::vector<dt_tuple>& dT_history);

  // extension of mesh API
  int BoundaryFaceGetCell(int f);

  // extension of STL
  void set_intersection(const std::vector<AmanziMesh::Entity_ID>& v1, 
                        const std::vector<AmanziMesh::Entity_ID>& v2, 
                        std::vector<AmanziMesh::Entity_ID>* vv);
 public:
  int ncells_owned, ncells_wghost;
  int nfaces_owned, nfaces_wghost;

  int MyPID;  // parallel information: will be moved to private
  int verbosity, verbosity_AztecOO;  // output information
  int missed_bc_faces_;
 
  Teuchos::RCP<Flow_State> FS;
  Teuchos::RCP<Flow_State> FS_aux;  // adds ghosts to selected state variables 
  
  int ti_phase_counter;  
  double T_physics, dT, dTnext;
  int flow_status_;
  int dim;

 private:
  int nseepage_prev;

 private:
  Teuchos::RCP<AmanziMesh::Mesh> mesh_;
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
