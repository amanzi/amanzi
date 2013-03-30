/*
This is the transport component of the Amanzi code. 

Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Author: Konstantin Lipnikov (lipnikov@lanl.gov)
Usage: 
  Transport_PK TPK(Teuchos::ParameterList& list, Teuchos::RCP<Transport_State> TS);
  double time_step = TPK.calculate_transport_dT();
  TPK.advance(time_step);
*/

#ifndef __Transport_PK_hpp__
#define __Transport_PK_hpp__

#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_Import.h"
#include "Teuchos_RCP.hpp"

#include "tensor.hh"
#include "Explicit_TI_FnBase.hh"
#include "BoundaryFunction.hh"

#include "State.hh"
#include "Transport_State.hh"
#include "Transport_SourceFactory.hh"
#include "Transport_constants.hh"
#include "Reconstruction.hh"

/*
This is Amanzi Transport Process Kernel (PK).

The transport PK receives a reduced (optional) copy of 
a physical state at time n and returns a different state 
at time n+1. 

Unmodified physical quantaties in the returned state are
the smart pointers to the original variables.
*/

namespace Amanzi {
namespace AmanziTransport {

double bestLSfit(const std::vector<double>& h, const std::vector<double>& error);

class Transport_PK : public Explicit_TI::fnBase {
 public:
  Transport_PK();
  Transport_PK(Teuchos::ParameterList& parameter_list_MPC,
               Teuchos::RCP<Transport_State> TS_MPC);
  ~Transport_PK() { for (int i=0; i<bcs.size(); i++) delete bcs[i]; }

  // primary members
  int InitPK();
  double EstimateTransportDt();
  double CalculateTransportDt();
  void Advance(double dT);
  void CommitState(Teuchos::RCP<Transport_State> TS) {};  // pointer to state is known

  void CheckDivergenceProperty();
  void CheckGEDproperty(Epetra_MultiVector& tracer) const; 
  void CheckTracerBounds(Epetra_MultiVector& tracer, int component,
                         double lower_bound, double upper_bound, double tol = 0.0) const;
  void CheckInfluxBC() const;

  // access members  
  Teuchos::RCP<Transport_State> transport_state() { return TS; }
  Teuchos::RCP<Transport_State> transport_state_next() { return TS_nextMPC; }
  Transport_State& ref_transport_state_next() { return *TS_nextBIG; }

  inline double cfl() { return cfl_; }
  inline int get_transport_status() { return status; }

  // control members
  void PrintStatistics() const;
  void WriteGMVfile(Teuchos::RCP<Transport_State> TS) const;

  // limiters
  void LimiterBarthJespersen(const int component,
                             Teuchos::RCP<Epetra_Vector> scalar_field, 
                             Teuchos::RCP<Epetra_MultiVector> gradient, 
                             Teuchos::RCP<Epetra_Vector> limiter);

  // sources and sinks
  void ComputeAddSourceTerms(double Tp, double dTp, 
                             DomainFunction* src_sink, Epetra_MultiVector& tcc);

 private:
  // advection members
  void AdvanceDonorUpwind(double dT);
  void AdvanceSecondOrderUpwindGeneric(double dT);
  void AdvanceSecondOrderUpwindRK1(double dT);
  void AdvanceSecondOrderUpwindRK2(double dT);

  // time integration members
  void fun(const double t, const Epetra_Vector& component, Epetra_Vector& f_component);

  // limiters 
  void LimiterTensorial(const int component,
                        Teuchos::RCP<Epetra_Vector> scalar_field, 
                        Teuchos::RCP<Epetra_MultiVector> gradient);

  void LimiterKuzmin(const int component,
                     Teuchos::RCP<Epetra_Vector> scalar_field, 
                     Teuchos::RCP<Epetra_MultiVector> gradient);

  void CalculateDescentDirection(std::vector<AmanziGeometry::Point>& normals,
                                 AmanziGeometry::Point& normal_new,
                                 double& L22normal_new, 
                                 AmanziGeometry::Point& direction);

  void ApplyDirectionalLimiter(AmanziGeometry::Point& normal, 
                               AmanziGeometry::Point& p,
                               AmanziGeometry::Point& direction, 
                               AmanziGeometry::Point& gradient);

  void IdentifyUpwindCells();

  const Teuchos::RCP<Epetra_IntVector>& get_upwind_cell() { return upwind_cell_; }
  const Teuchos::RCP<Epetra_IntVector>& get_downwind_cell() { return downwind_cell_; }  

  // dispersion routines
  void CalculateDispersionTensor();
  void ExtractBoundaryConditions(const int component,
                                 std::vector<int>& bc_face_id,
                                 std::vector<double>& bc_face_value);
  void PopulateHarmonicPointsValues(int component,
                                    Teuchos::RCP<Epetra_MultiVector> tcc,
                                    std::vector<int>& bc_face_id,
                                    std::vector<double>& bc_face_values);
  void AddDispersiveFluxes(int component,
                           Teuchos::RCP<Epetra_MultiVector> tcc,
                           std::vector<int>& bc_face_id,
                           std::vector<double>& bc_face_values,
                           Teuchos::RCP<Epetra_MultiVector> tcc_next);

  // I/O methods
  void ProcessParameterList();
  void ProcessStringFlowMode(const std::string name, int* method);
  void ProcessStringDispersionModel(const std::string name, int* method);
  void ProcessStringAdvectionLimiter(const std::string name, int* method);
  void ProcessStringVerbosity(const std::string name, int* verbosity);

  // obsolete methods
  void CreateConcentration(Teuchos::ParameterList& bcs_list);

  // miscaleneous methods
  double TracerVolumeChangePerSecond(int idx_tracer);

 public:
  int MyPID;  // parallel information: will be moved to private
  int spatial_disc_order, temporal_disc_order, limiter_model;

  int verbosity, internal_tests;  // output information
  double tests_tolerance;

 private:
  Teuchos::RCP<Transport_State> TS;
  Teuchos::RCP<Transport_State> TS_nextBIG;  // involves both owned and ghost values
  Teuchos::RCP<Transport_State> TS_nextMPC;  // uses physical memory of TS_nextBIG
  
  Teuchos::ParameterList parameter_list;

  Teuchos::RCP<Epetra_IntVector> upwind_cell_;
  Teuchos::RCP<Epetra_IntVector> downwind_cell_;

  Teuchos::RCP<Epetra_Vector> water_saturation_start;  // data for subcycling 
  Teuchos::RCP<Epetra_Vector> water_saturation_end;
  Teuchos::RCP<Epetra_Vector> ws_subcycle_start;  // ws = water saturation 
  Teuchos::RCP<Epetra_Vector> ws_subcycle_end;

  int advection_limiter;  // data for limiters
  int current_component_;
  Teuchos::RCP<Epetra_Vector> component_, component_next_;
  Teuchos::RCP<Epetra_Vector> limiter_;
  Reconstruction lifting;
  std::vector<double> component_local_min_;
  std::vector<double> component_local_max_;

  DomainFunction* src_sink;  // Source and sink terms
  int src_sink_distribution; 
  Teuchos::RCP<Epetra_Vector> Kxy;  // absolute permeability in plane xy

  Teuchos::RCP<Epetra_Import> cell_importer;  // parallel communicators
  Teuchos::RCP<Epetra_Import> face_importer;

  int dispersivity_model;  // data for dispersion 
  double dispersivity_longitudinal, dispersivity_transverse;

  std::vector<AmanziGeometry::Point> harmonic_points;
  std::vector<double> harmonic_points_weight;
  std::vector<double> harmonic_points_value;
  std::vector<WhetStone::Tensor> dispersion_tensor;

  double cfl_, dT, dT_debug, T_physics;  
  int number_components; 
  int status;
  int flow_mode;  // steady-sate or transient

  std::vector<BoundaryFunction*> bcs;  // influx BCs for each components
  std::vector<int> bcs_tcc_index; 
  double bc_scaling;

  double mass_tracer_exact;  // statistics for tracer

  int ncells_owned, ncells_wghost;
  int nfaces_owned, nfaces_wghost;
  int nnodes_wghost;
 
  Teuchos::RCP<AmanziMesh::Mesh> mesh_;
  int dim;
};

}  // namespace AmanziTransport
}  // namespace Amanzi

#endif

