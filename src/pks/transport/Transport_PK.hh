/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Daniil Svyatsky (dasvyat@lanl.gov)
*/

/*!

Mathematical models
...................
A few PDE models can be instantiated using the parameters described below.


Single-phase transport
``````````````````````
The conceptual PDE model for the transport in partially saturated media is

.. math::
  \frac{\partial (\phi s_l C_l)}{\partial t} 
  =
  - \boldsymbol{\nabla} \cdot (\boldsymbol{q}_l C_l) 
  + \boldsymbol{\nabla} \cdot (\phi_e s_l\, (\boldsymbol{D}_l + \tau \boldsymbol{M}_l) \boldsymbol{\nabla} C_l) + Q,

where 
:math:`\phi` is total porosity [-],
:math:`\phi_e` is effective transport porosity [-],
:math:`s_l` is liquid saturation [-], 
:math:`Q` is source or sink term,
:math:`\boldsymbol{q}_l` is the Darcy velocity [m/s],
:math:`\boldsymbol{D}_l` is dispersion tensor,
:math:`\boldsymbol{M}_l` is diffusion coefficient,
and :math:`\tau` is tortuosity [-].
For an isotropic medium with no preferred axis of symmetry the dispersion 
tensor has the following form:

.. math::
  \boldsymbol{D}_l 
  = \alpha_t \|\boldsymbol{v}\| \boldsymbol{I} 
  + \left(\alpha_l-\alpha_t \right) 
    \frac{\boldsymbol{v} \boldsymbol{v}}{\|\boldsymbol{v}\|}, \qquad
  \boldsymbol{v} = \frac{\boldsymbol{q}}{\phi_e}

where
:math:`\alpha_l` is longitudinal dispersivity [m],
:math:`\alpha_t` is  transverse dispersivity [m],
and :math:`\boldsymbol{v}` is average pore velocity [m/s].
Amanzi supports two additional models for dispersivity with 3 and 4 parameters.


Single-phase transport with dual porosity model
```````````````````````````````````````````````
The dual porosity formulation of the solute transport consists of two equations
for the fracture and matrix regions. 
In the fracture region, we have \citep{simunek-vangenuchten_2008}

.. math::
  \frac{\partial (\phi_f\, s_{lf}\, C_{lf})}{\partial t} 
  =
  - \boldsymbol{\nabla} \cdot (\boldsymbol{q}_l C_{lf}) 
  + \boldsymbol{\nabla} \cdot (\phi_f\, s_{lf}\, (\boldsymbol{D}_l + \tau_f M) \boldsymbol{\nabla} C_{lf}) 
  - \frac{\phi_m\,\tau_m}{L_m}\, M \nabla C_m - \Sigma_w C^* + Q_f,

where 
:math:`\phi_f` is fracture porosity [-],
:math:`\phi_m` is matrix porosity [-],
:math:`s_{lf}` is liquid saturation in fracture [-], 
:math:`\boldsymbol{q}_l` is the Darcy velocity [m/s],
:math:`\boldsymbol{D}_l` is dispersion tensor,
:math:`\tau_f` is fracture tortuosity [-],
:math:`\tau_m` is matrix tortuosity [-],
:math:`M` is molecular diffusion coefficient [:math:`m^2/s`], and
:math:`L_m` is the characteristic matrix depth defined typically as the ratio of a matrix block [m],
:math:`\Sigma_w` is transfer rate due to flow from the matrix to the fracture, 
:math:`C^*` is equal to :math:`C_{lf}` if :math:`\Sigma_w > 0` and :math:`C_{lm}` is :math:`\Sigma_w < 0`,
and :math:`Q_f` is source or sink term.
In the matrix region, we have

.. math::
  \frac{\partial (\phi_m\, s_{lm}\, C_{lm})}{\partial t}
  = \nabla\cdot (\phi_m\, \tau_m\, M_m \nabla C_{lm}) + \Sigma_w C^* + Q_m,

where 
:math:`\phi_m` is matrix porosity [-],
:math:`s_{lm}` is liquid saturation in matrix [-], 
:math:`Q_m` is source or sink term.
The simplified one-node dual porosity model uses a finite difference approximation of the 
solute gradient:

.. math::
  \nabla C_{lm} \approx WR \, \frac{C_{lf} - C_{lm}}{L_m},

where 
:math:`WR` is the Warren-Root coefficient that estimates the poro-space geometry, [-]


Physical models and assumptions
...............................
This list is used to summarize physical models and assumptions, such as
coupling with other PKs.
This list is often generated or extended by a high-level MPC PK.

* `"gas diffusion`" [bool] indicates that air-water partitioning coefficients
  are used to distribute components between liquid and as phases. Default is *false*.

* `"permeability field is required`" [bool] indicates if some transport features
  require absolute permeability. Default is *false*.

* `"multiscale model`" [string] specifies a multiscale model.
  Available options are `"single porosity`" (default) and `"dual porosity`".

* `"effective transport porosity`" [bool] If *true*, effective transport porosity
  will be used by dispersive-diffusive fluxes instead of total porosity. 
  Default is *false*.

* `"eos lookup table`" [string] provides the name for optional EOS lookup table.

* `"use dispersion solver`" [bool] instructs PK to instantiate a solver but do
  not call it. It is used now by MPC to form a global solver. Default is *false*.

.. code-block:: xml

  <ParameterList name="_TRANSPORT">  <!-- parent list -->
  <ParameterList name="physical models and assumptions">
    <Parameter name="gas diffusion" type="bool" value="false"/>
    <Parameter name="permeability field is required" type="bool" value="false"/>
    <Parameter name="multiscale model" type="string" value="single porosity"/>
    <Parameter name="effective transport porosity" type="bool" value="false"/>
    <Parameter name="eos lookup table" type="string" value="h2o.eos"/>
    <Parameter name="use dispersion solver" type="bool" value="false"/>
  </ParameterList>
  </ParameterList>


Global parameters
.................
This list is used to summarize physical models and assumptions, such as
The transport component of Amanzi performs advection of aqueous and gaseous
components and their dispersion and diffusion. 
The main parameters control temporal stability, spatial 
and temporal accuracy, and verbosity:


* `"domain name`" [string] specifies mesh name that defined domain of this PK.
  Default is `"domain`".

* `"cfl`" [double] Time step limiter, a number less than 1. Default value is 1.
   
* `"method`" [string] defines flux method. Available options are `"muscl`" (default) and `"fct`".
   
* `"spatial discretization order`" [int] defines accuracy of spatial discretization.
  It permits values 1 or 2. Default value is 1. 
  
* `"temporal discretization order`" [int] defines accuracy of temporal discretization.
  It permits values 1 or 2 and values 3 or 4. Note that RK3 is not monotone.
  Default value is 1.

* `"reconstruction`" [list] collects reconstruction parameters. The available options are
  describe in the separate section below.

* `"solver`" [string] Specifies the dispersion/diffusion solver.

* `"preconditioner`" [string] specifies preconditioner for dispersion solver.

* `"number of aqueous components`" [int] The total number of aqueous components. 
  Default value is the total number of components.

* `"number of gaseous components`" [int] The total number of gaseous components. 
  Default value is 0.

.. code-block:: xml

  <ParameterList>  <!-- parent list -->
  <ParameterList name="_TRANSPORT">
    <Parameter name="domain name" type="string" value="domain"/>
    <Parameter name="cfl" type="double" value="1.0"/>
    <Parameter name="method" type="string" value="muscl"/>
    <Parameter name="spatial discretization order" type="int" value="1"/>
    <Parameter name="temporal discretization order" type="int" value="1"/>
    <Parameter name="solver" type="string" value="_PCG_SOLVER"/>

    <ParameterList name="reconstruction">
      <Parameter name="method" type="string" value="cell-based"/>
      <Parameter name="polynomial order" type="int" value="1"/>
      <Parameter name="limiter" type="string" value="tensorial"/>
      <Parameter name="limiter extension for transport" type="bool" value="true"/>
    </ParameterList>

    <ParameterList name="verbose object">
      <Parameter name="verbosity level" type="string" value="high"/>
    </ParameterList>
  </ParameterList>  
  </ParameterList>  

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
#  include "Alquimia_PK.hh"
#  include "ChemistryEngine.hh"
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

typedef double
AnalyticFunction(const AmanziGeometry::Point&, const double);

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

  virtual ~Transport_PK(){};

  // members required by PK interface
  virtual void Setup() override;
  virtual void Initialize() override;

  virtual double get_dt() override;
  virtual void set_dt(double dt) override{};

  virtual void CommitStep(double t_old, double t_new, const Tag& tag) override;
  virtual void CalculateDiagnostics(const Tag& tag) override{};

  virtual std::string name() override { return "transport"; }

  // main transport members
  // -- calculation of a stable time step needs saturations and darcy flux
  double StableTimeStep(int n);

  // -- coupling with chemistry
  void SetupChemistry(Teuchos::RCP<AmanziChemistry::Chemistry_PK> chem_pk) { chem_pk_ = chem_pk; }
  void SetupAlquimia();

  // -- access members
  double cfl() { return cfl_; }
  bool get_flag_dispersion() { return flag_dispersion_ || flag_diffusion_; }
  Teuchos::RCP<CompositeVector> total_component_concentration() { return tcc_tmp; }
  void get_discretization_order(int* spatial, int* temporal)
  {
    *spatial = spatial_disc_order;
    *temporal = temporal_disc_order;
  }
  // -- molecualr diffusion coefficient for n-th solute
  double getDiffusion(int n) { return (diffusion_phase_[0]->values())[n]; }

  // -- modifiers
  void set_current_component(int i) { current_component_ = i; }

  // -- control members
  void VV_CheckGEDproperty(Epetra_MultiVector& tracer) const;
  void VV_CheckTracerBounds(Epetra_MultiVector& tracer,
                            int component,
                            double lower_bound,
                            double upper_bound,
                            double tol = 0.0) const;
  void VV_CheckInfluxBC() const;
  void VV_PrintSoluteExtrema(const Epetra_MultiVector& tcc_next,
                             double dT_MPC,
                             const std::string& mesh_id);
  double VV_SoluteVolumeChangePerSecond(int idx_solute);
  void VV_PrintLimiterStatistics();

  void CalculateLpErrors(AnalyticFunction f, double t, Epetra_Vector* sol, double* L1, double* L2);

  // multi-purpose wrapper for dispersion-diffusion solvers
  // -- solves for all components is comp0 < 0.
  // -- otherwise, it returns global operator for component comp0.
  Teuchos::RCP<Operators::Operator> DispersionSolver(const Epetra_MultiVector& tcc_prev,
                                                     Epetra_MultiVector& tcc_next,
                                                     double t_old,
                                                     double t_new,
                                                     int comp0 = -1);

 protected:
  void InitializeFields_();

  void FunctionalTimeDerivative_MUSCL_(double t,
                                       const CompositeVector& component,
                                       CompositeVector& f,
                                       bool scale);
  void
  FunctionalTimeDerivative_FCT_(double t, const CompositeVector& component, CompositeVector& f);

  // sources and sinks for components from n0 to n1 including
  void ComputeSources_(double tp,
                       double dtp,
                       Epetra_MultiVector& tcc,
                       const Epetra_MultiVector& tcc_prev,
                       int n0,
                       int n1);
  bool ComputeBCs_(std::vector<int>& bc_model, std::vector<double>& bc_value, int component);

  // tools
  void IdentifyUpwindCells();

  void InterpolateCellVector(const Epetra_MultiVector& v0,
                             const Epetra_MultiVector& v1,
                             double dT_int,
                             double dT,
                             Epetra_MultiVector& v_int);

  // physical models
  // -- dispersion and diffusion
  void CalculateDispersionTensor_(const Epetra_MultiVector& porosity,
                                  const Epetra_MultiVector& water_content);

  void CalculateDiffusionTensor_(double md,
                                 int phase,
                                 const Epetra_MultiVector& porosity,
                                 const Epetra_MultiVector& saturation,
                                 const Epetra_MultiVector& water_content);

  int FindDiffusionValue(const std::string& tcc_name, double* md, int* phase);

  void CalculateAxiSymmetryDirection();

  // -- effective diffusion
  void CalculateDiffusionTensorEffective_(double mdl,
                                          double mdg,
                                          double kH,
                                          const Epetra_MultiVector& porosity,
                                          const Epetra_MultiVector& saturation);

  void DiffusionSolverEffective(Epetra_MultiVector& tcc_next, double t_old, double t_new);

  // -- air-water partitioning using Henry's law. This is a temporary
  //    solution to get things moving.
  void PrepareAirWaterPartitioning_();
  void MakeAirWaterPartitioning_();

  // -- multiscale methods
  void AddMultiscalePorosity_(double t_old, double t_new, double t_int1, double t_int2);

  // initialization methods
  void InitializeAll_();

  // miscaleneous methods
  int FindComponentNumber(const std::string component_name);

 public:
  Teuchos::RCP<Teuchos::ParameterList> tp_list_;
  Teuchos::RCP<const Teuchos::ParameterList> preconditioner_list_;
  Teuchos::RCP<const Teuchos::ParameterList> linear_solver_list_;
  Teuchos::RCP<const Teuchos::ParameterList> nonlinear_solver_list_;

  int MyPID; // parallel information: will be moved to private
  int spatial_disc_order, temporal_disc_order, limiter_model;

  int nsubcycles; // output information
  bool internal_tests_;
  double internal_tests_tol_;

 protected:
  Teuchos::RCP<TreeVector> soln_;

  // names of state fields
  Key tcc_key_;
  Key vol_flowrate_key_, aperture_key_;
  Key porosity_key_, transport_porosity_key_, permeability_key_;
  Key saturation_liquid_key_, tortuosity_key_;
  Key wc_key_, prev_wc_key_;

  Key porosity_msp_key_;
  Key water_content_msp_key_, prev_water_content_msp_key_;

  int ncells_owned, ncells_wghost;
  int nfaces_owned, nfaces_wghost;
  int nnodes_wghost;

  Teuchos::RCP<CompositeVector> tcc_tmp; // next tcc
  Teuchos::RCP<CompositeVector> tcc;     // smart mirrow of tcc
  Teuchos::RCP<const Epetra_MultiVector> phi, transport_phi;

  Teuchos::RCP<const Epetra_MultiVector> wc_start, wc_end; // data for subcycling
  Teuchos::RCP<Epetra_MultiVector> wc_subcycle_start, wc_subcycle_end;

  std::vector<Teuchos::RCP<TransportDomainFunction>> srcs_; // sources and sinks
  std::vector<Teuchos::RCP<TransportDomainFunction>> bcs_;
  Teuchos::RCP<Epetra_MultiVector> Kxy; // absolute permeability in plane xy

  double cfl_, dt_, dt_debug_, t_physics_, dt_prev_;

  std::string passwd_;
  Method_t method_;

  bool transport_on_manifold_;
  bool subcycling_, use_transport_porosity_, use_effective_diffusion_;
  int dim;

  Teuchos::RCP<AmanziChemistry::Chemistry_PK> chem_pk_;
#ifdef ALQUIMIA_ENABLED
  Teuchos::RCP<AmanziChemistry::Alquimia_PK> alquimia_pk_;
  Teuchos::RCP<AmanziChemistry::ChemistryEngine> chem_engine_;
#endif

  std::vector<std::vector<int>> upwind_cells_; // fracture friendly
  std::vector<std::vector<int>> downwind_cells_;
  std::vector<std::vector<double>> upwind_flux_, downwind_flux_;

  int current_component_; // data for lifting
  Teuchos::RCP<Operators::ReconstructionCellLinear> lifting_;
  Teuchos::RCP<Operators::LimiterCell> limiter_;
  Teuchos::RCP<Operators::FCT> fct_;

  double limiter_mean_;

  Teuchos::RCP<Epetra_Import> cell_importer; // parallel communicators
  Teuchos::RCP<Epetra_Import> face_importer;

  // mechanical dispersion and molecular diffusion
  Teuchos::RCP<MDMPartition> mdm_;
  std::vector<WhetStone::Tensor> D_;

  bool flag_dispersion_, flag_diffusion_, use_dispersion_;
  std::vector<int> axi_symmetry_; // axi-symmetry direction of permeability tensor
  std::string dispersion_preconditioner, dispersion_solver;

  std::vector<Teuchos::RCP<MaterialProperties>> mat_properties_; // vector of materials
  std::vector<Teuchos::RCP<DiffusionPhase>> diffusion_phase_;    // vector of phases

  // Hosting temporarily Henry law
  bool henry_law_;
  std::vector<double> kH_;
  std::vector<int> air_water_map_;

  // multiscale models
  bool multiscale_porosity_;
  Teuchos::RCP<MultiscaleTransportPorosityPartition> msp_;

  std::vector<double> mass_solutes_exact_, mass_solutes_source_; // mass for all solutes
  std::vector<std::string> runtime_solutes_;                     // names of trached solutes
  std::vector<std::string> runtime_regions_;

  std::vector<std::string> component_names_; // details of components
  int num_aqueous, num_gaseous;

  // io
  Utils::Units units_;
};

} // namespace Transport
} // namespace Amanzi

#endif
