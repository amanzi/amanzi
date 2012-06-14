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

#ifndef __FLOW_PK_HPP__
#define __FLOW_PK_HPP__

#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_FECrsMatrix.h"
#include "Teuchos_RCP.hpp"

#include "boundary_function.hh"
#include "domain_function.hh"
#include "mfd3d.hpp"
#include "BDF2_fnBase.hpp"

#include "Flow_State.hpp"
#include "Matrix_MFD.hpp"

namespace Amanzi {
namespace AmanziFlow {

const int FLOW_STATUS_NULL = 0;  // used for internal debuging
const int FLOW_STATUS_INIT = 2;
const int FLOW_STATUS_STEADY_STATE_INIT = 4;
const int FLOW_STATUS_STEADY_STATE_COMPLETE = 6;
const int FLOW_STATUS_TRANSIENT_STATE_INIT = 8;
const int FLOW_STATUS_TRANSIENT_STATE_COMPLETE = 10;

const int FLOW_BC_FACE_NULL = 0; 
const int FLOW_BC_FACE_PRESSURE = 1; 
const int FLOW_BC_FACE_HEAD = 2; 
const int FLOW_BC_FACE_FLUX = 4; 

const int FLOW_TIME_INTEGRATION_PICARD = 1;
const int FLOW_TIME_INTEGRATION_BACKWARD_EULER = 2;  // Only for testing.
const int FLOW_TIME_INTEGRATION_BDF1 = 3;
const int FLOW_TIME_INTEGRATION_BDF2 = 4;
const double FLOW_INITIAL_DT = 1e-8;
const double FLOW_MAXIMAL_DT = 3.15e+10;  // 1000 years

const int FLOW_RELATIVE_PERM_CENTERED = 1; 
const int FLOW_RELATIVE_PERM_UPWIND_GRAVITY = 2; 
const int FLOW_RELATIVE_PERM_UPWIND_DARCY_FLUX = 3;
const int FLOW_RELATIVE_PERM_ARITHMETIC_MEAN = 4;
const double FLOW_RELATIVE_PERM_TOLERANCE = 1e-10;

const int FLOW_MFD3D_POLYHEDRA = 1;
const int FLOW_MFD3D_POLYHEDRA_MONOTONE = 2;  // under development
const int FLOW_MFD3D_HEXAHEDRA_MONOTONE = 3;
const int FLOW_MFD3D_TWO_POINT_FLUX = 4;  // without consistency
const int FLOW_MFD3D_SUPPORT_OPERATOR = 5;  // rc1 compatibility
const int FLOW_MFD3D_OPTIMIZED = 6;

const int FLOW_TI_ERROR_CONTROL_PRESSURE = 1;  // binary mask for error control
const int FLOW_TI_ERROR_CONTROL_SATURATION = 2;
const int FLOW_TI_ERROR_CONTROL_RESIDUAL = 4;

const double FLOW_TI_ABSOLUTE_TOLERANCE = 1e-4;  // defaults for time integrations
const double FLOW_TI_RELATIVE_TOLERANCE = 0.0;
const double FLOW_TI_NONLINEAR_RESIDUAL_TOLERANCE = 1e-6;
const int FLOW_TI_MAX_ITERATIONS = 400;

const int FLOW_HEX_FACES = 6;  // Hexahedron is the common element
const int FLOW_HEX_NODES = 8;
const int FLOW_HEX_EDGES = 12;

const int FLOW_QUAD_FACES = 4;  // Quadrilateral is the common element
const int FLOW_QUAD_NODES = 4;
const int FLOW_QUAD_EDGES = 4;

const int FLOW_MAX_FACES = 14;  // Kelvin's tetrakaidecahedron
const int FLOW_MAX_NODES = 47;  // These polyhedron parameters must
const int FLOW_MAX_EDGES = 60;  // be calculated in Init().

const int FLOW_INTERNAL_ERROR = 911;  // contact (lipnikov@lanl.gov)

const int FLOW_VERBOSITY_NONE = 0;
const int FLOW_VERBOSITY_LOW = 1;
const int FLOW_VERBOSITY_MEDIUM = 2;
const int FLOW_VERBOSITY_HIGH = 3;
const int FLOW_VERBOSITY_EXTREME = 4;

const int FLOW_AMANZI_VERSION = 2;  


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
  virtual int AdvanceToSteadyState() = 0;
  virtual void InitializeAuxiliaryData() = 0;
  virtual void InitializeSteadySaturated() = 0;

  virtual void CommitState(Teuchos::RCP<Flow_State> FS) = 0;

  // boundary condition members
  void ProcessBoundaryConditions(
      BoundaryFunction* bc_pressure, BoundaryFunction* bc_head,
      BoundaryFunction* bc_flux, BoundaryFunction* bc_seepage,
      const Epetra_Vector& pressure_faces, const double atm_pressure,
      std::vector<int>& bc_markers, std::vector<double>& bc_values);

  void ApplyBoundaryConditions(std::vector<int>& bc_markers,
                               std::vector<double>& bc_values,
                               Epetra_Vector& pressure_faces);

  void AddSourceTerms(DomainFunction* src_sink, Epetra_Vector& rhs);

  // gravity members
  void AddGravityFluxes_MFD(std::vector<WhetStone::Tensor>& K,
                            const Epetra_Vector& Krel_cells,
                            const Epetra_Vector& Krel_faces, 
                            Matrix_MFD* matrix);
  void AddGravityFluxes_DarcyFlux(std::vector<WhetStone::Tensor>& K,
                                  const Epetra_Vector& Krel_cells,
                                  const Epetra_Vector& Krel_faces,
                                  Epetra_Vector& darcy_mass_flux);

  // access members  
  Teuchos::RCP<Flow_State> flow_state() { return FS; }
  int flow_status() { return flow_status_; }

  // control members
  void ValidateBoundaryConditions(
      BoundaryFunction *bc_pressure, BoundaryFunction *bc_head, BoundaryFunction *bc_flux) const;
  void WriteGMVfile(Teuchos::RCP<Flow_State> FS) const;
 
  void set_time(double T0, double dT0) { T_physics = T0; dT = dT0; }
  void set_verbosity(int level) { verbosity = level; }
  
  // miscallenous members
  Epetra_Map* CreateSuperMap();
  void IdentifyUpwindCells(Epetra_IntVector& upwind_cell, Epetra_IntVector& downwind_cell);

  // io members
  void ProcessStringMFD3D(const std::string name, int* method);
  void ProcessStringVerbosity(const std::string name, int* verbosity);

 public:
  int ncells_owned, ncells_wghost;
  int nfaces_owned, nfaces_wghost;

  int MyPID;  // parallel information: will be moved to private
  int verbosity, internal_tests;  // output information
 
  Teuchos::RCP<Flow_State> FS;
  
  double T_physics, dT, dTnext;
  int flow_status_;
 
  int dim;

 private:
  Teuchos::RCP<AmanziMesh::Mesh> mesh_;
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
