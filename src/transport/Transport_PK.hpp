/*
This is the transport component of the Amanzi code. 
License: BSD
Author: Konstantin Lipnikov (lipnikov@lanl.gov)
Usage: 
  Transport_PK TPK(Teuchos::ParameterList& list, Teuchos::RCP<Transport_State> TS);
  double time_step = TPK.calculate_transport_dT();
  TPK.advance( time_step);
*/

#ifndef __Transport_PK_hpp__
#define __Transport_PK_hpp__

#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_Import.h"
#include "Teuchos_RCP.hpp"

#include "tensor.hpp"
#include "Explicit_TI_fnBase.hpp"

#include "State.hpp"
#include "Transport_State.hpp"
#include "Transport_BCs.hpp"
#include "Reconstruction.hpp"

/*
This is Amanzi Transport Process Kernel (PK), release Beta.

The transport PK receives a reduced (optional) copy of 
a physical state at time n and returns a different state 
at time n+1. 

Unmodified physical quantaties in the returned state are
the smart pointers to the original variables.
*/

namespace Amanzi {
namespace AmanziTransport {

const int TRANSPORT_NULL = 0;
const int TRANSPORT_FLOW_AVAILABLE = 1;
const int TRANSPORT_STATE_BEGIN = 2;
const int TRANSPORT_STATE_COMPLETE = 3;

const int TRANSPORT_INTERNAL_ERROR = 911;  // contact (lipnikov@lanl.gov)

const double TRANSPORT_LARGE_TIME_STEP = 1e+99;
const double TRANSPORT_SMALL_TIME_STEP = 1e-12;

const int TRANSPORT_BC_CONSTANT_TCC = 1;
const int TRANSPORT_BC_DISPERSION_FLUX = 2;
const int TRANSPORT_BC_NULL = 3;

const double TRANSPORT_CONCENTRATION_OVERSHOOT = 1e-6;

const int TRANSPORT_MAX_FACES = 14;  // Kelvin's tetrakaidecahedron
const int TRANSPORT_MAX_NODES = 47;  // These olyhedron parameters must
const int TRANSPORT_MAX_EDGES = 60;  // be calculated in Init().

const int TRANSPORT_DISPERSIVITY_MODEL_NULL = 1;
const int TRANSPORT_DISPERSIVITY_MODEL_ISOTROPIC = 2;
const int TRANSPORT_DISPERSIVITY_MODEL_BEAR = 3;
const int TRANSPORT_DISPERSIVITY_MODEL_LICHTNER = 4;

const int TRANSPORT_LIMITER_BARTH_JESPERSEN = 1; 
const int TRANSPORT_LIMITER_TENSORIAL = 2;

const int TRANSPORT_AMANZI_VERSION = 2;  


class Transport_PK : public Explicit_TI::fnBase {
 public:
  Transport_PK();
  Transport_PK(Teuchos::ParameterList& parameter_list_MPC,
               Teuchos::RCP<Transport_State> TS_MPC);
  ~Transport_PK() {};

  // primary members
  double calculate_transport_dT();
  void advance(double dT);
  void commit_state(Teuchos::RCP<Transport_State> TS) {};  // pointer to state is known

  void check_divergence_property();
  void check_GEDproperty(Epetra_MultiVector& tracer) const; 
  void check_tracer_bounds(Epetra_MultiVector& tracer, 
                           int component,
                           double lower_bound,
                           double upper_bound,
                           double tol = 0.0) const;

  // access members  
  Teuchos::RCP<Transport_State> get_transport_state() { return TS; }
  Teuchos::RCP<Transport_State> get_transport_state_next() { return TS_nextMPC; }
  Transport_State& ref_transport_state_next() { return *TS_nextBIG; }

  inline double get_transport_dT() { return dT; }
  inline double get_cfl() { return cfl; }
  inline int get_transport_status() { return status; }

  void print_statistics() const;
 
 private:
  // advection routines
  void advance_donor_upwind(double dT);
  void advance_second_order_upwind(double dT);
  void advance_second_order_upwind_testing(double dT);
  void advance_arbitrary_order_upwind(double dT);
  void fun(const double t, const Epetra_Vector& component, Epetra_Vector& f_component);

  void limiterBarthJespersen(const int component,
                             Teuchos::RCP<Epetra_Vector> scalar_field, 
                             Teuchos::RCP<Epetra_MultiVector> gradient, 
                             std::vector<double>& field_local_min,
                             std::vector<double>& field_local_max,
                             Teuchos::RCP<Epetra_Vector> limiter);

  void limiterTensorial(const int component,
                        Teuchos::RCP<Epetra_Vector> scalar_field, 
                        Teuchos::RCP<Epetra_MultiVector> gradient);

  void apply_directional_limiter(AmanziGeometry::Point& normal, 
                                 AmanziGeometry::Point& p,
                                 AmanziGeometry::Point& direction, 
                                 AmanziGeometry::Point& gradient);

  void process_parameter_list();
  void identify_upwind_cells();

  const Teuchos::RCP<Epetra_IntVector>& get_upwind_cell() { return upwind_cell_; }
  const Teuchos::RCP<Epetra_IntVector>& get_downwind_cell() { return downwind_cell_; }  

  // dispersion routines
  void calculate_dispersion_tensor();
  void extract_boundary_conditions(const int component,
                                   std::vector<int>& bc_face_id,
                                   std::vector<double>& bc_face_value);
  void populate_harmonic_points_values(int component,
                                       Teuchos::RCP<Epetra_MultiVector> tcc,
                                       std::vector<int>& bc_face_id,
                                       std::vector<double>& bc_face_values);
  void add_dispersive_fluxes(int component,
                             Teuchos::RCP<Epetra_MultiVector> tcc,
                             std::vector<int>& bc_face_id,
                             std::vector<double>& bc_face_values,
                             Teuchos::RCP<Epetra_MultiVector> tcc_next);

 public:
  std::vector<double> calculate_accumulated_influx();
  std::vector<double> calculate_accumulated_outflux();

  int MyPID;  // parallel information: will be moved to private
  int spatial_disc_order, temporal_disc_order, limiter_model;

  int verbosity_level, internal_tests;  // output information
  double tests_tolerance;

 private:
  int Init();  // move here all code from constructor 

  Teuchos::RCP<Transport_State> TS;
  Teuchos::RCP<Transport_State> TS_nextBIG;  // involves both owned and ghost values
  Teuchos::RCP<Transport_State> TS_nextMPC;  // uses physical memory of TS_nextBIG
  
  Teuchos::ParameterList parameter_list;

  Teuchos::RCP<Epetra_IntVector> upwind_cell_;
  Teuchos::RCP<Epetra_IntVector> downwind_cell_;

  int current_component_;
  Teuchos::RCP<Epetra_Vector> component_, component_next_;
  Teuchos::RCP<Epetra_Vector> limiter_;
  Reconstruction lifting;

  Teuchos::RCP<Epetra_Import> cell_importer;  // parallel communicators
  Teuchos::RCP<Epetra_Import> face_importer;

  int dispersivity_model;  // data for dispersion 
  double dispersivity_longitudinal, dispersivity_transverse;

  std::vector<AmanziGeometry::Point> harmonic_points;
  std::vector<double> harmonic_points_weight;
  std::vector<double> harmonic_points_value;
  std::vector<WhetStone::Tensor> dispersion_tensor;

  int advection_limiter;  // data for limiters

  double cfl, dT, dT_debug, T_internal;  
  int number_components; 
  int status;

  std::vector<Transport_BCs> bcs;  // BCs for each components and each side set

  int cmin, cmax_owned, cmax, number_owned_cells, number_wghost_cells;
  int fmin, fmax_owned, fmax, number_owned_faces, number_wghost_faces;
  int vmin, vmax;
 
  Teuchos::RCP<AmanziMesh::Mesh> mesh_;
  int dim;
};

}  // namespace AmanziTransport
}  // namespace Amanzi

#endif

