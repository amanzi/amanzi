/*
  This is the transport component of Amanzi. 

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

#ifndef AMANZI_TRANSPORT_PK_HH_
#define AMANZI_TRANSPORT_PK_HH_

#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_Import.h"
#include "Teuchos_RCP.hpp"

#include "CompositeVector.hh"
#include "DiffusionPhase.hh"
#include "Explicit_TI_FnBase.hh"
#include "MaterialProperties.hh"
#include "PK.hh"
#include "ReconstructionCell.hh"
#include "State.hh"
#include "tensor.hh"
#include "TransportBoundaryFunction.hh"
#include "TransportDomainFunction.hh"
#include "TransportDefs.hh"
#include "TransportSourceFactory.hh"
#include "VerboseObject.hh"

#ifdef ALQUIMIA_ENABLED
#include "Chemistry_State.hh"
#include "ChemistryEngine.hh"
#endif

/*
  This is Amanzi Transport Process Kernel (PK).

  The transport PK receives a reduced (optional) copy of 
  a physical state at time n and returns a different state 
  at time n+1. 

  Unmodified physical quantaties in the returned state are
  the smart pointers to the original variables.
*/

namespace Amanzi {
namespace Transport {

double bestLSfit(const std::vector<double>& h, const std::vector<double>& error);
typedef double AnalyticFunction(const AmanziGeometry::Point&, const double);

class Transport_PK : public Explicit_TI::fnBase<Epetra_Vector> {
 public:
  Transport_PK();
#ifdef ALQUIMIA_ENABLED
  Transport_PK(Teuchos::ParameterList& glist, Teuchos::RCP<State> S,
               Teuchos::RCP<AmanziChemistry::Chemistry_State> chem_state = Teuchos::null,
               Teuchos::RCP<AmanziChemistry::ChemistryEngine> chem_engine = Teuchos::null);
#endif
  Transport_PK(Teuchos::ParameterList& glist, Teuchos::RCP<State> S,
               std::vector<std::string>& component_names);
  ~Transport_PK();

  // main PK members
  void Initialize(const Teuchos::Ptr<State>& S);
  int Advance(double dT, double &dT_actual); 
  void CommitState(double dummy_dT, const Teuchos::Ptr<State>& S);

  void SetState(const Teuchos::RCP<State>& S) { S_ = S; }
  double get_dt();
  std::string name() { return passwd_; }

  // main transport members
  void InitializeFields();
  double CalculateTransportDt();

  // access members  
  Teuchos::RCP<const State> state() { return S_; }
  Teuchos::RCP<CompositeVector> total_component_concentration() { return tcc_tmp; }

  double TimeStep() { return dT; }
  void TimeStep(double dT_) { dT = dT_; }
  inline double cfl() { return cfl_; }

  // control members
  void Policy(Teuchos::Ptr<State> S);
  void PrintStatistics() const;
  void WriteGMVfile(Teuchos::RCP<State> S) const;

  void CheckDivergenceProperty();
  void CheckGEDproperty(Epetra_MultiVector& tracer) const; 
  void CheckTracerBounds(Epetra_MultiVector& tracer, int component,
                         double lower_bound, double upper_bound, double tol = 0.0) const;
  void CheckInfluxBC() const;

  void CreateDefaultState(Teuchos::RCP<const AmanziMesh::Mesh>& mesh, int ncomponents);

  void CalculateLpErrors(AnalyticFunction f, double t, Epetra_Vector* sol, double* L1, double* L2);

  // sources and sinks for components from n0 to n1 including
  void ComputeAddSourceTerms(double Tp, double dTp, 
                             std::vector<TransportDomainFunction*>& src_sink, 
                             Epetra_MultiVector& tcc, int n0, int n1);

  bool PopulateBoundaryData(std::vector<int>& bc_model,
                            std::vector<double>& bc_value, int component);

  // limiters 
  void LimiterBarthJespersen(const int component,
                             Teuchos::RCP<const Epetra_Vector> scalar_field, 
                             Teuchos::RCP<CompositeVector>& gradient, 
                             Teuchos::RCP<Epetra_Vector>& limiter);

 private:
  // Helper for constructors.
  void Construct_(Teuchos::ParameterList& glist, 
                  Teuchos::RCP<State> S,
                  std::vector<std::string>& component_names);

  // advection members
  void AdvanceDonorUpwind(double dT);
  void AdvanceSecondOrderUpwindGeneric(double dT);
  void AdvanceSecondOrderUpwindRK1(double dT);
  void AdvanceSecondOrderUpwindRK2(double dT);

  // time integration members
  void Functional(const double t, const Epetra_Vector& component, Epetra_Vector& f_component);

  // limiters 
  void LimiterTensorial(const int component,
                        Teuchos::RCP<const Epetra_Vector> scalar_field, 
                        Teuchos::RCP<CompositeVector>& gradient);

  void LimiterKuzmin(const int component,
                     Teuchos::RCP<const Epetra_Vector>& scalar_field, 
                     Teuchos::RCP<CompositeVector>& gradient);

  void CalculateDescentDirection(std::vector<AmanziGeometry::Point>& normals,
                                 AmanziGeometry::Point& normal_new,
                                 double& L22normal_new, 
                                 AmanziGeometry::Point& direction);

  void ApplyDirectionalLimiter(AmanziGeometry::Point& normal, 
                               AmanziGeometry::Point& p,
                               AmanziGeometry::Point& direction, 
                               AmanziGeometry::Point& gradient);

  void IdentifyUpwindCells();

  void InterpolateCellVector(
      const Epetra_MultiVector& v0, const Epetra_MultiVector& v1, 
      double dT_int, double dT, Epetra_MultiVector& v_int);

  const Teuchos::RCP<Epetra_IntVector>& upwind_cell() { return upwind_cell_; }
  const Teuchos::RCP<Epetra_IntVector>& downwind_cell() { return downwind_cell_; }  

  // dispersion and diffusion
  void CalculateDispersionTensor_(
      const Epetra_MultiVector& darcy_flux, 
      const Epetra_MultiVector& porosity, const Epetra_MultiVector& saturation);

  void CalculateDiffusionTensor_(
      double md, int phase,
      const Epetra_MultiVector& porosity, const Epetra_MultiVector& saturation);

  int FindDiffusionValue(const std::string& tcc_name, double* md, int* phase);

  // I/O methods
  void ProcessParameterList();
  void ProcessStringDispersionModel(const std::string name, int* model);
  void ProcessStringDispersionMethod(const std::string name, int* method);
  void ProcessStringAdvectionLimiter(const std::string name, int* method);

  // miscaleneous methods
  int FindComponentNumber(const std::string component_name);
  double TracerVolumeChangePerSecond(int idx_tracer);

 public:
  Teuchos::ParameterList parameter_list;
  Teuchos::ParameterList preconditioners_list;
  Teuchos::ParameterList solvers_list;
  Teuchos::ParameterList nonlin_solvers_list;

  int MyPID;  // parallel information: will be moved to private
  int spatial_disc_order, temporal_disc_order, limiter_model;

  int nsubcycles;  // output information
  int internal_tests;
  double tests_tolerance;

 private:
  VerboseObject* vo_;

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  int dim;

  Teuchos::RCP<State> S_;  // state info
  std::string passwd_;

  Teuchos::RCP<CompositeVector> tcc_tmp;  // next tcc
  Teuchos::RCP<CompositeVector> tcc;  // smart mirrow of tcc 
  Teuchos::RCP<Epetra_MultiVector> darcy_flux;
  Teuchos::RCP<const Epetra_MultiVector> ws, ws_prev, phi;
  
#ifdef ALQUIMIA_ENABLED
  Teuchos::RCP<AmanziChemistry::Chemistry_State> chem_state_;
  Teuchos::RCP<AmanziChemistry::ChemistryEngine> chem_engine_;
#endif

  Teuchos::RCP<Epetra_IntVector> upwind_cell_;
  Teuchos::RCP<Epetra_IntVector> downwind_cell_;

  Teuchos::RCP<const Epetra_MultiVector> ws_start, ws_end;  // data for subcycling 
  Teuchos::RCP<Epetra_MultiVector> ws_subcycle_start, ws_subcycle_end;

  int advection_limiter;  // data for limiters
  int current_component_;
  Teuchos::RCP<Epetra_Vector> limiter_;
  Teuchos::RCP<Operators::ReconstructionCell> lifting_;
  std::vector<double> component_local_min_;
  std::vector<double> component_local_max_;

  std::vector<TransportDomainFunction*> srcs;  // Source or sink for components
  std::vector<TransportBoundaryFunction*> bcs;  // influx BC for components
  double bc_scaling;
  Teuchos::RCP<Epetra_Vector> Kxy;  // absolute permeability in plane xy

  Teuchos::RCP<Epetra_Import> cell_importer;  // parallel communicators
  Teuchos::RCP<Epetra_Import> face_importer;

  std::vector<Teuchos::RCP<MaterialProperties> > material_properties_;  // vector of materials
  std::vector<Teuchos::RCP<DiffusionPhase> >  diffusion_phase_;   // vector of phases

  std::vector<WhetStone::Tensor> D;
  int dispersion_models_;
  std::string dispersion_preconditioner;
  std::string dispersion_solver;

  double cfl_, dT, dT_debug, T_physics;  

  double mass_tracer_exact, mass_tracer_source;  // statistics for tracer

  int ncells_owned, ncells_wghost;
  int nfaces_owned, nfaces_wghost;
  int nnodes_wghost;
 
  std::vector<std::string> component_names_;  // details of components
  int num_aqueous, num_gaseous;

  // Forbidden.
  Transport_PK(const Transport_PK&);
  Transport_PK& operator=(const Transport_PK&);
};

}  // namespace Transport
}  // namespace Amanzi

#endif

