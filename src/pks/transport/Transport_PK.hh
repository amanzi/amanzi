/*
  Transport PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
          Daniil Svyatsky (dasvyat@lanl.gov)
*/

#ifndef AMANZI_TRANSPORT_PK_HH_
#define AMANZI_TRANSPORT_PK_HH_

// TPLs
#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_Import.h"
#include "Teuchos_RCP.hpp"

// Amanzi
#include "CompositeVector.hh"
#include "DenseVector.hh"
#include "FCT.hh"
#include "Key.hh"
#include "LimiterCell.hh"
#include "PK.hh"
#include "PK_Explicit.hh"
#include "PK_Factory.hh"
#include "PK_Physical.hh"
#include "ReconstructionCellLinear.hh"
#include "State.hh"
#include "Tensor.hh"
#include "Units.hh"
#include "VerboseObject.hh"

#include "Chemistry_PK.hh"
#ifdef ALQUIMIA_ENABLED
#include "Alquimia_PK.hh"
#include "ChemistryEngine.hh"
#endif

// Amanzi::Transport
#include "DiffusionPhase.hh"
#include "MaterialProperties.hh"
#include "MDMPartition.hh"
#include "MultiscaleTransportPorosityPartition.hh"
#include "TransportDefs.hh"
#include "TransportDomainFunction.hh"

namespace Amanzi {
namespace Transport {

typedef double AnalyticFunction(const AmanziGeometry::Point&, const double);

class Transport_PK : public PK_Physical {  
  public:
    Transport_PK(Teuchos::ParameterList& pk_tree,
                 const Teuchos::RCP<Teuchos::ParameterList>& glist,
                 const Teuchos::RCP<State>& S,
                 const Teuchos::RCP<TreeVector>& soln);

    Transport_PK(const Teuchos::RCP<Teuchos::ParameterList>& glist,
                 Teuchos::RCP<State> S,
                 const std::string& pk_list_name,
                 std::vector<std::string>& component_names);

    virtual ~Transport_PK() {};

  // members required by PK interface
  virtual void Setup() override;
  virtual void Initialize() override;

  virtual double get_dt() override;
  virtual void set_dt(double dt) override {};

  virtual void CommitStep(double t_old, double t_new, const Tag& tag) override;
  virtual void CalculateDiagnostics(const Tag& tag) override {};

  virtual std::string name() override { return "transport"; }

  // main transport members
  // -- calculation of a stable time step needs saturations and darcy flux
  double StableTimeStep();

  // -- coupling with chemistry
  void SetupChemistry(Teuchos::RCP<AmanziChemistry::Chemistry_PK> chem_pk) { chem_pk_ = chem_pk; }
  void SetupAlquimia();

  // -- access members  
  inline double cfl() { return cfl_; }
  Teuchos::RCP<const State> state() { return S_; }
  Teuchos::RCP<CompositeVector> total_component_concentration() { return tcc_tmp; }

  // -- control members
  void CreateDefaultState(Teuchos::RCP<const AmanziMesh::Mesh>& mesh, int ncomponents);
  void Policy(Teuchos::Ptr<State> S);

  void VV_CheckGEDproperty(Epetra_MultiVector& tracer) const; 
  void VV_CheckTracerBounds(Epetra_MultiVector& tracer, int component,
                            double lower_bound, double upper_bound, double tol = 0.0) const;
  void VV_CheckInfluxBC() const;
  void VV_PrintSoluteExtrema(const Epetra_MultiVector& tcc_next, double dT_MPC, const std::string& mesh_id);
  double VV_SoluteVolumeChangePerSecond(int idx_solute);
  void VV_PrintLimiterStatistics();

  void CalculateLpErrors(AnalyticFunction f, double t, Epetra_Vector* sol, double* L1, double* L2);

 protected:
  void InitializeFields_();

  // sources and sinks for components from n0 to n1 including
  void ComputeSources_(double tp, double dtp, Epetra_MultiVector& tcc,
                       const Epetra_MultiVector& tcc_prev, int n0, int n1);
  bool ComputeBCs_(std::vector<int>& bc_model, std::vector<double>& bc_value, int component);

  // tools
  void IdentifyUpwindCells();

  void InterpolateCellVector(
      const Epetra_MultiVector& v0, const Epetra_MultiVector& v1, 
      double dT_int, double dT, Epetra_MultiVector& v_int);

  // physical models
  // -- dispersion and diffusion
  void CalculateDispersionTensor_(
      const Epetra_MultiVector& darcy_flux, 
      const Epetra_MultiVector& porosity, const Epetra_MultiVector& saturation);

  void CalculateDiffusionTensor_(
      double md, int phase,
      const Epetra_MultiVector& porosity, const Epetra_MultiVector& saturation);

  int FindDiffusionValue(const std::string& tcc_name, double* md, int* phase);

  void CalculateAxiSymmetryDirection();

  void DispersionSolver(const Epetra_MultiVector& tcc_prev,
                        const Epetra_MultiVector& tcc_next,
                        double t_old, double t_new);

  // -- effective diffusion
  void CalculateDiffusionTensorEffective_(
      double mdl, double mdg, double kH,
      const Epetra_MultiVector& porosity, const Epetra_MultiVector& saturation);

  void DiffusionSolverEffective(const Epetra_MultiVector& tcc_next,
                                double t_old, double t_new);

  // -- air-water partitioning using Henry's law. This is a temporary
  //    solution to get things moving.
  void PrepareAirWaterPartitioning_();
  void MakeAirWaterPartitioning_();

  // -- multiscale methods
  void AddMultiscalePorosity_(double t_old, double t_new, double t_int1, double t_int2);

  // initialization methods
  void InitializeAll_();
  void InitializeFieldFromField_(const std::string& field0, const std::string& field1, bool call_evaluator);

  // miscaleneous methods
  int FindComponentNumber(const std::string component_name);

 public:
  Teuchos::RCP<Teuchos::ParameterList> tp_list_;
  Teuchos::RCP<const Teuchos::ParameterList> preconditioner_list_;
  Teuchos::RCP<const Teuchos::ParameterList> linear_solver_list_;
  Teuchos::RCP<const Teuchos::ParameterList> nonlinear_solver_list_;

  int MyPID;  // parallel information: will be moved to private
  int spatial_disc_order, temporal_disc_order, limiter_model;

  int nsubcycles;  // output information
  bool internal_tests_, genericRK_;
  double internal_tests_tol_;

 protected:
  Teuchos::RCP<TreeVector> soln_;

  // names of state fields 
  Key tcc_key_;
  Key darcy_flux_key_;
  Key porosity_key_, transport_porosity_key_, permeability_key_;
  Key saturation_liquid_key_, prev_saturation_liquid_key_;
  Key water_content_key_, prev_water_content_key_;

  Key porosity_msp_key_;
  Key water_content_msp_key_, prev_water_content_msp_key_;

  int ncells_owned, ncells_wghost;
  int nfaces_owned, nfaces_wghost;
  int nnodes_wghost;

  Teuchos::RCP<CompositeVector> tcc_tmp;  // next tcc
  Teuchos::RCP<CompositeVector> tcc;  // smart mirrow of tcc 
  Teuchos::RCP<const Epetra_MultiVector> ws, ws_prev, phi, transport_phi;
    
  Teuchos::RCP<const Epetra_MultiVector> ws_start, ws_end;  // data for subcycling 
  Teuchos::RCP<Epetra_MultiVector> ws_subcycle_start, ws_subcycle_end;

  std::vector<Teuchos::RCP<TransportDomainFunction> > srcs_;  // Sources and sinks
  std::vector<Teuchos::RCP<TransportDomainFunction> > bcs_;
  Teuchos::RCP<Epetra_Vector> Kxy;  // absolute permeability in plane xy

  double cfl_, dt_, dt_debug_, t_physics_;  

  std::string passwd_;
  Method_t method_;

  bool subcycling_, use_transport_porosity_, use_effective_diffusion_;
  int dim;

  Teuchos::RCP<AmanziChemistry::Chemistry_PK> chem_pk_;
#ifdef ALQUIMIA_ENABLED
  Teuchos::RCP<AmanziChemistry::Alquimia_PK> alquimia_pk_;
  Teuchos::RCP<AmanziChemistry::ChemistryEngine> chem_engine_;
#endif

  std::vector<std::vector<int> > upwind_cells_;  // fracture friendly 
  std::vector<std::vector<int> > downwind_cells_;
  std::vector<std::vector<double> > upwind_flux_, downwind_flux_;

  int current_component_;  // data for lifting
  Teuchos::RCP<Operators::ReconstructionCellLinear> lifting_;
  Teuchos::RCP<Operators::LimiterCell> limiter_;
  Teuchos::RCP<Operators::FCT> fct_;

  Teuchos::RCP<Epetra_Import> cell_importer;  // parallel communicators
  Teuchos::RCP<Epetra_Import> face_importer;

  // mechanical dispersion and molecual diffusion
  Teuchos::RCP<MDMPartition> mdm_;
  std::vector<WhetStone::Tensor> D_;

  bool flag_dispersion_;
  std::vector<int> axi_symmetry_;  // axi-symmetry direction of permeability tensor
  std::string dispersion_preconditioner, dispersion_solver;

  std::vector<Teuchos::RCP<MaterialProperties> > mat_properties_;  // vector of materials
  std::vector<Teuchos::RCP<DiffusionPhase> > diffusion_phase_;   // vector of phases

  // Hosting temporarily Henry law 
  bool henry_law_;
  std::vector<double> kH_;
  std::vector<int> air_water_map_;

  // multiscale models
  bool multiscale_porosity_;
  Teuchos::RCP<MultiscaleTransportPorosityPartition> msp_;
 
  std::vector<double> mass_solutes_exact_, mass_solutes_source_;  // mass for all solutes
  std::vector<std::string> runtime_solutes_;  // names of trached solutes
  std::vector<std::string> runtime_regions_;

  std::vector<std::string> component_names_;  // details of components
  int num_aqueous, num_gaseous;

  // io
  Utils::Units units_;

  // Forbidden.
  Transport_PK(const Transport_PK&);
  Transport_PK& operator=(const Transport_PK&);
};

}  // namespace Transport
}  // namespace Amanzi

#endif

