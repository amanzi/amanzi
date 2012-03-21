/*
This is the flow component of the Amanzi code. 
License: BSD
Authors: Neil Carlson (version 1) 
         Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
Usage: 
  Flow_PK TPK(Teuchos::ParameterList& list, Teuchos::RCP<Flow_State> FS);
  double time_step = FPK.calculate_flow_dT();
  FPK.advance(time_step);
*/

#ifndef __FLOW_PK_HPP__
#define __FLOW_PK_HPP__

#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_FECrsMatrix.h"
#include "Teuchos_RCP.hpp"

#include "boundary-function.hh"
#include "mfd3d.hpp"
#include "BDF2_fnBase.hpp"

#include "Flow_State.hpp"
#include "Matrix_MFD.hpp"

namespace Amanzi {
namespace AmanziFlow {

const int FLOW_NULL = 0;
const int FLOW_INIT_COMPLETE = 1;
const int FLOW_STEADY_STATE_COMPLETE = 2;
const int FLOW_NEXT_STATE_BEGIN = 3;
const int FLOW_NEXT_STATE_COMPLETE = 4;

const int FLOW_BC_FACE_NULL = 0; 
const int FLOW_BC_FACE_PRESSURE = 1; 
const int FLOW_BC_FACE_HEAD = 2; 
const int FLOW_BC_FACE_FLUX = 4; 

const int FLOW_STEADY_STATE_PICARD = 1;
const int FLOW_STEADY_STATE_BACKWARD_EULER = 2;
const int FLOW_STEADY_STATE_BDF2 = 3;
const int FLOW_STEADY_STATE_MAX_ITERATIONS = 100;
const double FLOW_STEADY_STATE_TOLERANCE = 1e-6;
const double FLOW_STEADY_STATE_INITIAL_DT = 1e-8;

const int FLOW_RELATIVE_PERM_CENTERED = 1; 
const int FLOW_RELATIVE_PERM_UPWIND_GRAVITY = 2; 
const int FLOW_RELATIVE_PERM_UPWIND_DARCY_FLUX = 3;
const int FLOW_RELATIVE_PERM_ARITHMETIC_MEAN = 4; 

const int FLOW_MFD3D_POLYHEDRA = 1;
const int FLOW_MFD3D_POLYHEDRA_MONOTONE = 2;  // under development
const int FLOW_MFD3D_HEXAHEDRA_MONOTONE = 3;

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

const int FLOW_VERBOSITY_NULL = 0;
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
  virtual void InitPK(Matrix_MFD* matrix_ = NULL, Matrix_MFD* preconditioner_ = NULL) = 0;
  virtual void InitSteadyState(double T0, double dT0) = 0;
  virtual void InitTransient(double T0, double dT0) = 0;

  virtual int advance(double dT) = 0; 
  virtual int advanceToSteadyState() = 0;
  virtual void commitState(Teuchos::RCP<Flow_State> FS) = 0;
  virtual void commitStateForTransport(Teuchos::RCP<Flow_State> FS) = 0;

  double calculate_flow_dT() { return dT; }

  virtual void deriveDarcyVelocity(const Epetra_Vector& flux, Epetra_MultiVector& velocity) = 0;

  // boundary condition members
  void updateBoundaryConditions(
      BoundaryFunction* bc_pressure, BoundaryFunction* bc_head,
      BoundaryFunction* bc_flux, BoundaryFunction* bc_seepage,
      const Epetra_Vector& pressure_cells, const double atm_pressure,
      std::vector<int>& bc_markers, std::vector<double>& bc_values);

  void applyBoundaryConditions(std::vector<int>& bc_markers,
                               std::vector<double>& bc_values,
                               Epetra_Vector& pressure_faces);

  // gravity members
  void addGravityFluxes_MFD(std::vector<WhetStone::Tensor>& K, 
                            const Epetra_Vector& Krel_faces, 
                            Matrix_MFD* matrix);
  void addGravityFluxes_DarcyFlux(std::vector<WhetStone::Tensor>& K, 
                                  const Epetra_Vector& Krel_faces,
                                  Epetra_Vector& darcy_mass_flux);

  // access members  
  Teuchos::RCP<Flow_State> flow_state() { return FS; }
  int flow_status() { return flow_status_; }

  // control members
  void validate_boundary_conditions(
      BoundaryFunction *bc_pressure, BoundaryFunction *bc_head, BoundaryFunction *bc_flux) const;
  inline void set_standalone_mode(bool mode) { standalone_mode = mode; }
  void writeGMVfile(Teuchos::RCP<Flow_State> FS) const;
 
  void set_time(double T0, double dT0) { T_internal = T0; dT = dT0; }
  void set_verbosity(int level) { verbosity = level; }

  // miscallenous members
  Epetra_Map* createSuperMap();
  void identifyUpwindCells(Epetra_IntVector& upwind_cell, Epetra_IntVector& downwind_cell);

 public:
  int ncells_owned, ncells_wghost;
  int nfaces_owned, nfaces_wghost;

  int MyPID;  // parallel information: will be moved to private
  int verbosity, internal_tests;  // output information
 
  Teuchos::RCP<Flow_State> FS;
  
  double T_internal, dT, dT0;
  int flow_status_;
  int standalone_mode;
 
  Teuchos::RCP<AmanziMesh::Mesh> mesh_;
  int dim;
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
