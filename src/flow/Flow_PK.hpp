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
#define __FLOw_PK_HPP__

#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_FECrsMatrix.h"
#include "Teuchos_RCP.hpp"

#include "boundary-function.hh"
#include "mfd3d.hpp"
#include "Flow_State.hpp"

namespace Amanzi {
namespace AmanziFlow {

const int FLOW_NULL = 0;
const int FLOW_STATE_BEGIN = 2;
const int FLOW_STATE_INIT = 3;
const int FLOW_STATE_COMPLETE = 4;

const int FLOW_BC_FACE_NULL = 0; 
const int FLOW_BC_FACE_PRESSURE = 1; 
const int FLOW_BC_FACE_HEAD = 2; 
const int FLOW_BC_FACE_FLUX = 4; 

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

const int FLOW_AMANZI_VERSION = 2;  


class Flow_PK {
 public:
  Flow_PK() {};
  ~Flow_PK() {};

  // main methods
  void Init(Teuchos::RCP<Flow_State> FS_MPC);

  virtual int advance(double dT) = 0; 
  virtual int advance_to_steady_state() = 0;
  virtual void commit_state() = 0;

  void update_data_boundary_faces(BoundaryFunction *bc_pressure, 
                                  BoundaryFunction *bc_head,
                                  BoundaryFunction *bc_flux,
                                  std::vector<int>& bc_markers,
                                  std::vector<double>& bc_values);

  void calculateGravityFluxes(int cell, WhetStone::Tensor& K, std::vector<double>& gravity_flux);

  // access members  
  Teuchos::RCP<Flow_State> get_flow_state() { return FS; }
  Teuchos::RCP<Flow_State> get_flow_state_next() { return FS_nextMPC; }
  Flow_State& ref_flow_state_next() { return *FS_nextBIG; }

  Teuchos::RCP<AmanziMesh::Mesh> get_mesh() { return mesh_; }

  inline double get_flow_dT() { return dT; }
  inline int get_flow_status() { return status; }

  // control members
  void validate_boundary_conditions(
      BoundaryFunction *bc_pressure, BoundaryFunction *bc_head, BoundaryFunction *bc_flux) const;
  inline bool set_standalone_mode(bool mode) { standalone_mode = mode; } 

  // miscallenous members
  Epetra_Map* create_super_map();

 public:
  int cmin, cmax_owned, cmax, number_owned_cells, number_wghost_cells;
  int fmin, fmax_owned, fmax, number_owned_faces, number_wghost_faces;
  int vmin, vmax;

  int MyPID;  // parallel information: will be moved to private
  int verbosity_level, internal_tests;  // output information
 
  Teuchos::RCP<Flow_State> FS;
  Teuchos::RCP<Flow_State> FS_nextBIG;  // involves both owned and ghost values
  Teuchos::RCP<Flow_State> FS_nextMPC;  // uses physical memory of FS_nextBIG
  
  double T_internal, T_physical, dT;
  int status;
  int standalone_mode;
 
  Teuchos::RCP<AmanziMesh::Mesh> mesh_;
  int dim;
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
