/*
  Flow PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Neil Carlson (nnc@lanl.gov), 
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <vector>

// TPLs
#include "Epetra_IntVector.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

// Amanzi
#include "BoundaryFlux.hh"
#include "CommonDefs.hh"
#include "dbc.hh"
#include "exceptions.hh"
#include "independent_variable_field_evaluator_fromfunction.hh"
#include "Mesh.hh"
#include "OperatorDefs.hh"
#include "PDE_DiffusionFactory.hh"
#include "PK_Utils.hh"
#include "Point.hh"
#include "primary_variable_field_evaluator.hh"
#include "UpwindFactory.hh"
#include "XMLParameterListWriter.hh"
#include "LinearOperatorFactory.hh"

// Amanzi::Flow
#include "DarcyVelocityEvaluator.hh"
#include "PorosityModelEvaluator.hh"
#include "RelPermEvaluator.hh"
#include "Richards_PK.hh"
#include "VWContentEvaluator.hh"
#include "VWContentEvaluatorFactory.hh"
#include "WRMEvaluator.hh"
#include "WRM.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* Simplest possible constructor: extracts lists and requires fields.
****************************************************************** */
Richards_PK::Richards_PK(Teuchos::ParameterList& pk_tree,
                         const Teuchos::RCP<Teuchos::ParameterList>& glist,
                         const Teuchos::RCP<State>& S,
                         const Teuchos::RCP<TreeVector>& soln) :
  Flow_PK(pk_tree, glist, S, soln),
  glist_(glist),
  soln_(soln)
{
  S_ = S;

  std::string pk_name = pk_tree.name();
  auto found = pk_name.rfind("->");
  if (found != std::string::npos) pk_name.erase(0, found + 2);

  // We need the flow list
  Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(glist, "PKs", true);
  fp_list_ = Teuchos::sublist(pk_list, pk_name, true);
  
  // We also need miscaleneous sublists
  preconditioner_list_ = Teuchos::sublist(glist, "preconditioners", true);
  linear_operator_list_ = Teuchos::sublist(glist, "solvers", true);
  ti_list_ = Teuchos::sublist(fp_list_, "time integrator");

  // domain name
  domain_ = fp_list_->template get<std::string>("domain name", "domain");

  vo_ = Teuchos::null;
}


/* ******************************************************************
* Old constructor for unit tests.
****************************************************************** */
Richards_PK::Richards_PK(const Teuchos::RCP<Teuchos::ParameterList>& glist,
                         const std::string& pk_list_name,
                         Teuchos::RCP<State> S,
                         const Teuchos::RCP<TreeVector>& soln) :
    Flow_PK(),
    glist_(glist),
    soln_(soln)
{
  S_ = S;

  // We need the flow list
  Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(glist, "PKs", true);
  fp_list_ = Teuchos::sublist(pk_list, pk_list_name, true);
 
  // We also need miscaleneous sublists
  preconditioner_list_ = Teuchos::sublist(glist, "preconditioners", true);
  linear_operator_list_ = Teuchos::sublist(glist, "solvers", true);
  ti_list_ = Teuchos::sublist(fp_list_, "time integrator");

  // domain name
  domain_ = fp_list_->template get<std::string>("domain name", "domain");

  ms_itrs_ = 0;
  ms_calls_ = 0;

  vo_ = Teuchos::null;
}


/* ******************************************************************
* Clean memory.
****************************************************************** */
Richards_PK::~Richards_PK()
{
  if (vo_ != Teuchos::null) vo_ = Teuchos::null;
}


/* ******************************************************************
* Define structure of this PK. We request physical fields and their
* evaluators. Selection of a few models is available and driven by
* model factories, evaluator factories, and parameters of the list
* "physical models and assumptions".
****************************************************************** */
void Richards_PK::Setup(const Teuchos::Ptr<State>& S)
{
  dt_ = 0.0;
  mesh_ = S->GetMesh(domain_);
  dim = mesh_->space_dimension();

  // generate keys here to be available for setup of the base class
  pressure_key_ = Keys::getKey(domain_, "pressure"); 
  hydraulic_head_key_ = Keys::getKey(domain_, "hydraulic_head"); 

  darcy_flux_key_ = Keys::getKey(domain_, "darcy_flux"); 
  darcy_velocity_key_ = Keys::getKey(domain_, "darcy_velocity"); 

  permeability_key_ = Keys::getKey(domain_, "permeability"); 
  porosity_key_ = Keys::getKey(domain_, "porosity"); 

  saturation_liquid_key_ = Keys::getKey(domain_, "saturation_liquid"); 
  prev_saturation_liquid_key_ = Keys::getKey(domain_, "prev_saturation_liquid"); 

  // set up the base class 
  Flow_PK::Setup(S);

  // Our decision can be affected by the list of models
  Teuchos::RCP<Teuchos::ParameterList> physical_models =
      Teuchos::sublist(fp_list_, "physical models and assumptions");
  std::string vwc_model = physical_models->get<std::string>("water content model", "constant density");
  std::string multiscale_model = physical_models->get<std::string>("multiscale model", "single continuum");

  // Require primary field for this PK, which is pressure
  std::vector<std::string> names;
  std::vector<AmanziMesh::Entity_kind> locations;
  std::vector<int> ndofs;

  Teuchos::RCP<Teuchos::ParameterList> list1 = Teuchos::sublist(fp_list_, "operators", true);
  Teuchos::RCP<Teuchos::ParameterList> list2 = Teuchos::sublist(list1, "diffusion operator", true);
  Teuchos::RCP<Teuchos::ParameterList> list3 = Teuchos::sublist(list2, "matrix", true);
  std::string name = list3->get<std::string>("discretization primary");

  names.push_back("cell");
  locations.push_back(AmanziMesh::CELL);
  ndofs.push_back(1);
  if (name != "fv: default" && name != "nlfv: default") {
    names.push_back("face");
    locations.push_back(AmanziMesh::FACE);
    ndofs.push_back(1);
  }

  if (!S->HasField(pressure_key_)) {
    S->RequireField(pressure_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponents(names, locations, ndofs);

    Teuchos::ParameterList elist;
    elist.set<std::string>("evaluator name", pressure_key_);
    pressure_eval_ = Teuchos::rcp(new PrimaryVariableFieldEvaluator(elist));
    S->SetFieldEvaluator(pressure_key_, pressure_eval_);
  }

  // Require conserved quantity.
  // -- water content
  if (!S->HasField("water_content")) {
    S->RequireField("water_content", "water_content")->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);

    Teuchos::ParameterList vwc_list;
    vwc_list.set<std::string>("pressure key", pressure_key_)
            .set<std::string>("saturation key", saturation_liquid_key_)
            .set<std::string>("porosity key", porosity_key_);
    VWContentEvaluatorFactory fac;
    Teuchos::RCP<VWContentEvaluator> eval = fac.Create(vwc_model, vwc_list);
    S->SetFieldEvaluator("water_content", eval);
  }

  // -- water content from the previous time step
  if (!S->HasField("prev_water_content")) {
    S->RequireField("prev_water_content", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    S->GetField("prev_water_content", passwd_)->set_io_vis(false);
  }

  // -- multiscale extension: secondary (immobile water content)
  if (multiscale_model == "dual continuum discontinuous matrix") {
    if (!S->HasField("pressure_matrix")) {
      S->RequireField("pressure_matrix", passwd_)->SetMesh(mesh_)->SetGhosted(false)
        ->SetComponent("cell", AmanziMesh::CELL, 1);

      Teuchos::ParameterList elist;
      elist.set<std::string>("evaluator name", "pressure_matrix");
      pressure_matrix_eval_ = Teuchos::rcp(new PrimaryVariableFieldEvaluator(elist));
      S->SetFieldEvaluator("pressure_matrix", pressure_matrix_eval_);
    }

    Teuchos::RCP<Teuchos::ParameterList> msp_list = Teuchos::sublist(fp_list_, "multiscale models", true);
    msp_ = CreateMultiscaleFlowPorosityPartition(mesh_, msp_list);

    if (!S->HasField("water_content_matrix")) {
      S->RequireField("water_content_matrix", passwd_)->SetMesh(mesh_)->SetGhosted(true)
        ->SetComponent("cell", AmanziMesh::CELL, 1);
    }
    if (!S->HasField("prev_water_content_matrix")) {
      S->RequireField("prev_water_content_matrix", passwd_)->SetMesh(mesh_)->SetGhosted(true)
        ->SetComponent("cell", AmanziMesh::CELL, 1);
      S->GetField("prev_water_content_matrix", passwd_)->set_io_vis(false);
    }

    S->RequireField("porosity_matrix", "porosity_matrix")->SetMesh(mesh_)->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, 1);

    Teuchos::ParameterList elist;
    elist.set<std::string>("porosity key", "porosity_matrix");
    elist.set<std::string>("pressure key", "pressure_matrix");
    Teuchos::RCP<PorosityModelPartition> pom = CreatePorosityModelPartition(mesh_, msp_list);
    Teuchos::RCP<PorosityModelEvaluator> eval = Teuchos::rcp(new PorosityModelEvaluator(elist, pom));
    S->SetFieldEvaluator("porosity_matrix", eval);

    // Secondary matrix nodes are collected here.
    int nnodes, nnodes_tmp = NumberMatrixNodes(msp_);
    mesh_->get_comm()->MaxAll(&nnodes_tmp, &nnodes, 1);
    if (nnodes > 1) {
      S->RequireField("water_content_matrix_aux", passwd_)
        ->SetMesh(mesh_)->SetGhosted(false)
        ->SetComponent("cell", AmanziMesh::CELL, nnodes - 1);
      S->GetField("water_content_matrix_aux", passwd_)->set_io_vis(false);

      S->RequireField("pressure_matrix_aux", passwd_)->SetMesh(mesh_)->SetGhosted(false)
        ->SetComponent("cell", AmanziMesh::CELL, nnodes - 1);
      S->GetField("pressure_matrix_aux", passwd_)->set_io_vis(false);

      S->RequireField("porosity_matrix_aux", "porosity_matrix_aux")->SetMesh(mesh_)->SetGhosted(false)
        ->SetComponent("cell", AmanziMesh::CELL, nnodes - 1);
      S->SetFieldEvaluator("porosity_matrix_aux", eval);
      S->GetField("porosity_matrix_aux", passwd_)->set_io_vis(false);
    }
  }

  // Require additional fields and evaluators for this PK.
  // -- absolute permeability
  if (!S->HasField(permeability_key_)) {
    S->RequireField(permeability_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, dim);
  }

  // -- darcy flux
  if (!S->HasField(darcy_flux_key_)) {
    S->RequireField(darcy_flux_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("face", AmanziMesh::FACE, 1);

    Teuchos::ParameterList elist;
    elist.set<std::string>("evaluator name", darcy_flux_key_);
    darcy_flux_eval_ = Teuchos::rcp(new PrimaryVariableFieldEvaluator(elist));
    S->SetFieldEvaluator(darcy_flux_key_, darcy_flux_eval_);
  }

  // -- porosity
  if (!S->HasField(porosity_key_)) {
    S->RequireField(porosity_key_, porosity_key_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);

    std::string pom_name = physical_models->get<std::string>("porosity model", "constant porosity");

    if (pom_name == "compressible: pressure function") {
      Teuchos::RCP<Teuchos::ParameterList>
          pom_list = Teuchos::sublist(fp_list_, "porosity models", true);
      Teuchos::RCP<PorosityModelPartition> pom = CreatePorosityModelPartition(mesh_, pom_list);

      Teuchos::ParameterList elist;
      elist.set<std::string>("porosity key", porosity_key_)
           .set<std::string>("pressure key", pressure_key_);
      // elist.sublist("verbose object").set<std::string>("verbosity level", "extreme");
      Teuchos::RCP<PorosityModelEvaluator> eval = Teuchos::rcp(new PorosityModelEvaluator(elist, pom));
      S->SetFieldEvaluator(porosity_key_, eval);
    } else {
      S->RequireFieldEvaluator(porosity_key_);
    }
  }

  // -- viscosity: if not requested by any PK, we request its constant value.
  if (!S->HasField("viscosity_liquid")) {
    if (!S->HasField("fluid_viscosity")) {
      S->RequireScalar("fluid_viscosity", passwd_);
    }
    S->RequireField("viscosity_liquid", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    S->GetField("viscosity_liquid", passwd_)->set_io_vis(false);
  }

  // -- model for liquid density is constant density unless specified otherwise
  //    in high-level PKs.
  if (!S->HasField("molar_density_liquid")) {
    S->RequireField("molar_density_liquid", "molar_density_liquid")->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);

    double rho = glist_->sublist("state").sublist("initial conditions")
                        .sublist("fluid_density").get<double>("value", 1000.0);
    double n_l = rho / CommonDefs::MOLAR_MASS_H2O;

    Teuchos::ParameterList& wc_eval = S->FEList().sublist("molar_density_liquid");
    wc_eval.sublist("function").sublist("DOMAIN")
           .set<std::string>("region", "All")
           .set<std::string>("component","cell")
           .sublist("function").sublist("function-constant")
           .set<double>("value", n_l);
    wc_eval.set<std::string>("field evaluator type", "independent variable");

    S->RequireFieldEvaluator("molar_density_liquid");
  }
  
  // -- saturation
  Teuchos::RCP<Teuchos::ParameterList>
      wrm_list = Teuchos::sublist(fp_list_, "water retention models", true);
  wrm_ = CreateWRMPartition(mesh_, wrm_list);

  if (!S->HasField(saturation_liquid_key_)) {
    S->RequireField(saturation_liquid_key_, saturation_liquid_key_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);

    Teuchos::ParameterList elist;
    elist.set<std::string>("saturation key", saturation_liquid_key_);
    // elist.sublist("verbose object").set<std::string>("verbosity level", "extreme");
    Teuchos::RCP<WRMEvaluator> eval = Teuchos::rcp(new WRMEvaluator(elist, wrm_));
    S->SetFieldEvaluator(saturation_liquid_key_, eval);
  }

  if (!S->HasField(prev_saturation_liquid_key_)) {
    S->RequireField(prev_saturation_liquid_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    S->GetField(prev_saturation_liquid_key_, passwd_)->set_io_vis(false);
  }

  // Local fields and evaluators.
  // -- hydraulic head
  if (!S->HasField(hydraulic_head_key_)) {
    S->RequireField(hydraulic_head_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  // -- Darcy velocity vector
  if (!S->HasField(darcy_velocity_key_)) {
    S->RequireField(darcy_velocity_key_, darcy_velocity_key_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, dim);

    Teuchos::ParameterList elist;
    elist.set<std::string>("domain name", domain_);
    elist.set<std::string>("darcy velocity key", darcy_velocity_key_)
         .set<std::string>("darcy flux key", darcy_flux_key_);
    Teuchos::RCP<DarcyVelocityEvaluator> eval = Teuchos::rcp(new DarcyVelocityEvaluator(elist));
    S->SetFieldEvaluator(darcy_velocity_key_, eval);
  }

  // Require additional components for the existing fields
  Teuchos::ParameterList abs_perm = fp_list_->sublist("absolute permeability");
  coordinate_system_ = abs_perm.get<std::string>("coordinate system", "cartesian");
  int noff = abs_perm.get<int>("off-diagonal components", 0);
 
  if (noff > 0) {
    CompositeVectorSpace& cvs = *S->RequireField(permeability_key_, passwd_);
    cvs.SetOwned(false);
    cvs.AddComponent("offd", AmanziMesh::CELL, noff)->SetOwned(true);
  }
}


/* ******************************************************************
* This is a long but simple routine. It goes through flow parameter
* list and initializes various objects including those created during 
* the setup step.
****************************************************************** */
void Richards_PK::Initialize(const Teuchos::Ptr<State>& S)
{
  // Initialize miscalleneous defaults.
  // -- times
  double t_ini = S->time(); 
  dt_desirable_ = dt_;
  dt_next_ = dt_;

  // -- others
  error_control_ = FLOW_TI_ERROR_CONTROL_PRESSURE;

  mass_bc = 0.0;
  seepage_mass_ = 0.0;
  mass_initial = 0.0;
  initialize_with_darcy_ = true;
  num_itrs_ = 0;

  // Create verbosity object to print out initialiation statistics.
  Teuchos::ParameterList vlist;
  vlist.sublist("verbose object") = fp_list_->sublist("verbose object");
  vo_ = Teuchos::rcp(new VerboseObject("FlowPK::Richards", vlist)); 

  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os()<< "\nPK initialization started...\n";
  }

  // Initilize various base class data.
  Flow_PK::Initialize(S);

  // Initialize local fields and evaluators.
  InitializeFields_();
  UpdateLocalFields_(S);

  // Create BCs and source terms.
  InitializeBCsSources_(*fp_list_);

  // relative permeability
  // -- create basic fields, factories and control variables
  Teuchos::RCP<Teuchos::ParameterList> upw_list = Teuchos::sublist(fp_list_, "relative permeability", true);
  relperm_ = Teuchos::rcp(new RelPerm(*upw_list, mesh_, atm_pressure_, wrm_));

  Operators::UpwindFactory<RelPerm> upwind_factory;
  upwind_ = upwind_factory.Create(mesh_, relperm_, *upw_list);

  std::string upw_upd = upw_list->get<std::string>("upwind frequency", "every timestep");
  if (upw_upd == "every nonlinear iteration") upwind_frequency_ = FLOW_UPWIND_UPDATE_ITERATION;
  else upwind_frequency_ = FLOW_UPWIND_UPDATE_TIMESTEP;  

  // relative permeability and related stractures
  // -- create vectors using estimate of the space size
  Teuchos::RCP<CompositeVectorSpace> upw_cvs = upwind_->Map();
  krel_ = Teuchos::rcp(new CompositeVector(*upw_cvs));
  dKdP_ = Teuchos::rcp(new CompositeVector(*upw_cvs));

  // -- populate fields with default values
  krel_->PutScalarMasterAndGhosted(1.0);
  dKdP_->PutScalarMasterAndGhosted(0.0);

  // Process models and assumptions.
  flux_units_ = molar_rho_ / rho_;

  // -- coupling with other physical PKs
  Teuchos::RCP<Teuchos::ParameterList> physical_models = 
      Teuchos::sublist(fp_list_, "physical models and assumptions");
  vapor_diffusion_ = physical_models->get<bool>("vapor diffusion", false);
  multiscale_porosity_ = (physical_models->get<std::string>(
      "multiscale model", "single continuum") != "single continuum");

  // Process other fundamental structures.
  SetAbsolutePermeabilityTensor();

  // Select a proper matrix class. 
  const Teuchos::ParameterList& tmp_list = fp_list_->sublist("operators")
                                                    .sublist("diffusion operator");
  Teuchos::ParameterList oplist_matrix = tmp_list.sublist("matrix");
  Teuchos::ParameterList oplist_pc = tmp_list.sublist("preconditioner");

  std::string name = fp_list_->sublist("relative permeability").get<std::string>("upwind method");
  std::string nonlinear_coef("standard: cell");
  if (name == "upwind: darcy velocity") {
    nonlinear_coef = "upwind: face";
  } else if (name == "upwind: gravity") {
    nonlinear_coef = "upwind: face";
  } else if (name == "upwind: amanzi" || name == "upwind: amanzi new") {
    nonlinear_coef = "divk: cell-face";
    // nonlinear_coef = "divk: face";
  } else if (name == "other: arithmetic average") {
    nonlinear_coef = "upwind: face";
  }
  oplist_matrix.set<std::string>("nonlinear coefficient", nonlinear_coef);
  oplist_pc.set<std::string>("nonlinear coefficient", nonlinear_coef);

  Operators::PDE_DiffusionFactory opfactory;
  op_matrix_diff_ = opfactory.Create(oplist_matrix, mesh_, op_bc_, rho_, gravity_);
  op_matrix_ = op_matrix_diff_->global_operator();
  op_preconditioner_diff_ = opfactory.Create(oplist_pc, mesh_, op_bc_, rho_, gravity_);
  op_preconditioner_ = op_preconditioner_diff_->global_operator();
  op_acc_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, op_preconditioner_));

  if (vapor_diffusion_) {
    Teuchos::ParameterList oplist_vapor = tmp_list.sublist("vapor matrix");
    op_vapor_diff_ = opfactory.Create(oplist_vapor, mesh_, op_bc_);
    op_vapor_ = op_vapor_diff_->global_operator();
    op_preconditioner_->OpPushBack(op_vapor_diff_->local_op(),
                                   Operators::OPERATOR_PROPERTY_DATA_READ_ONLY);
  }

  // Create pointers to the primary flow field pressure.
  solution = S->GetFieldData(pressure_key_, passwd_);
  soln_->SetData(solution); 
  
  // Create auxiliary vectors for time history and error estimates.
  const Epetra_BlockMap& cmap_owned = mesh_->cell_map(false);
  pdot_cells_prev = Teuchos::rcp(new Epetra_Vector(cmap_owned));
  pdot_cells = Teuchos::rcp(new Epetra_Vector(cmap_owned));

  // Initialize flux copy for the upwind operator.
  darcy_flux_copy = Teuchos::rcp(new CompositeVector(*S->GetFieldData(darcy_flux_key_, passwd_)));

  // Conditional initialization of lambdas from pressures.
  CompositeVector& pressure = *S->GetFieldData(pressure_key_, passwd_);

  if (ti_list_->isSublist("pressure-lambda constraints") && pressure.HasComponent("face")) {
    DeriveFaceValuesFromCellValues(*pressure.ViewComponent("cell"),
                                   *pressure.ViewComponent("face"));
  }

  // error control options
  AMANZI_ASSERT(ti_list_->isParameter("error control options"));

  error_control_ = 0;
  std::vector<std::string> options;
  options = ti_list_->get<Teuchos::Array<std::string> >("error control options").toVector();

  for (int i=0; i < options.size(); i++) {
    if (options[i] == "pressure") {
      error_control_ += FLOW_TI_ERROR_CONTROL_PRESSURE;
    } else if (options[i] == "saturation") {
      error_control_ += FLOW_TI_ERROR_CONTROL_SATURATION;
    } else if (options[i] == "residual") {
      error_control_ += FLOW_TI_ERROR_CONTROL_RESIDUAL;
    }
  }

  if (error_control_ == 0) {
    error_control_ = FLOW_TI_ERROR_CONTROL_PRESSURE +  // usually 1 [Pa]
                     FLOW_TI_ERROR_CONTROL_SATURATION;  // usually 1e-4;
  }

  // initialize time integrator
  std::string ti_method_name = ti_list_->get<std::string>("time integration method", "none");
  AMANZI_ASSERT(ti_method_name == "BDF1");
  Teuchos::ParameterList& bdf1_list = ti_list_->sublist("BDF1");

  if (! bdf1_list.isSublist("verbose object"))
      bdf1_list.sublist("verbose object") = fp_list_->sublist("verbose object");

  bdf1_dae = Teuchos::rcp(new BDF1_TI<TreeVector, TreeVectorSpace>(*this, bdf1_list, soln_));

  // Initialize boundary conditions and source terms.
  UpdateSourceBoundaryData(t_ini, t_ini, pressure);

  // Initialize matrix and preconditioner operators.
  // -- setup phase
  // -- molar density requires to rescale gravity later.
  op_matrix_->Init();
  Teuchos::RCP<std::vector<WhetStone::Tensor> > Kptr = Teuchos::rcpFromRef(K);
  op_matrix_diff_->SetBCs(op_bc_, op_bc_);
  op_matrix_diff_->Setup(Kptr, krel_, dKdP_);

  op_preconditioner_->Init();
  op_preconditioner_diff_->SetBCs(op_bc_, op_bc_);
  op_preconditioner_diff_->Setup(Kptr, krel_, dKdP_);

  // -- assemble phase
  op_matrix_diff_->UpdateMatrices(Teuchos::null, solution.ptr());
  op_matrix_diff_->ApplyBCs(true, true, true);

  op_preconditioner_diff_->UpdateMatrices(darcy_flux_copy.ptr(), solution.ptr());
  op_preconditioner_diff_->UpdateMatricesNewtonCorrection(darcy_flux_copy.ptr(), solution.ptr(), molar_rho_);
  op_preconditioner_diff_->ApplyBCs(true, true, true);
  op_preconditioner_->SymbolicAssembleMatrix();

  if (vapor_diffusion_) {
    // op_vapor_diff_->SetBCs(op_bc_);
    op_vapor_diff_->SetScalarCoefficient(Teuchos::null, Teuchos::null);
  }

  // -- generic linear solver for most cases
  solver_name_ = ti_list_->get<std::string>("linear solver");

  // -- preconditioner or encapsulated preconditioner
  std::string pc_name = ti_list_->get<std::string>("preconditioner");
  Teuchos::ParameterList pc_list = preconditioner_list_->sublist(pc_name);
  op_preconditioner_->InitializePreconditioner(pc_list);
  
  op_pc_solver_ = op_preconditioner_;

  if (ti_list_->isParameter("preconditioner enhancement")) {
    std::string tmp_solver = ti_list_->get<std::string>("preconditioner enhancement");
    if (tmp_solver != "none") {
      AMANZI_ASSERT(linear_operator_list_->isSublist(tmp_solver));

      AmanziSolvers::LinearOperatorFactory<Operators::Operator, CompositeVector, CompositeVectorSpace> sfactory;
      op_pc_solver_ = sfactory.Create(tmp_solver, *linear_operator_list_, op_preconditioner_);
    }
  }
  
  // Optional step: calculate hydrostatic solution consistent with BCs
  // and clip it as requested. We have to do it only once at the beginning
  // of time period.
  if (ti_list_->isSublist("initialization") && initialize_with_darcy_ 
      && S->position() == Amanzi::TIME_PERIOD_START) {
    initialize_with_darcy_ = false;
    Teuchos::ParameterList& ini_list = ti_list_->sublist("initialization");
 
    std::string ini_method_name = ini_list.get<std::string>("method", "none");
    if (ini_method_name == "saturated solver") {
      name = ini_list.get<std::string>("linear solver");
      SolveFullySaturatedProblem(t_ini, *solution, name);

      bool clip(false);
      double clip_saturation = ini_list.get<double>("clipping saturation value", -1.0);
      if (clip_saturation > 0.0) {
        double pmin = atm_pressure_;
        Epetra_MultiVector& p = *solution->ViewComponent("cell");
        ClipHydrostaticPressure(pmin, clip_saturation, p);
        clip = true;
      }

      double clip_pressure = ini_list.get<double>("clipping pressure value", -1e+10);
      if (clip_pressure > -5 * atm_pressure_) {
        Epetra_MultiVector& p = *solution->ViewComponent("cell");
        ClipHydrostaticPressure(clip_pressure, p);
        clip = true;
      }

      if (clip && vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
        Teuchos::OSTab tab = vo_->getOSTab();
        *vo_->os() << "\nClipped pressure field.\n" << std::endl;
      }

      if (clip && solution->HasComponent("face")) {
        Epetra_MultiVector& p = *solution->ViewComponent("cell");
        Epetra_MultiVector& lambda = *solution->ViewComponent("face", true);
        DeriveFaceValuesFromCellValues(p, lambda);
      }
    }
    else if (ini_method_name == "picard") {
      AdvanceToSteadyState_Picard(ti_list_->sublist("initialization"));
    }
    pressure_eval_->SetFieldAsChanged(S.ptr());

    // initialization is usually done at time 0, so we need to update other
    // fields such as prev_saturation_liquid
    S->GetFieldEvaluator(saturation_liquid_key_)->HasFieldChanged(S.ptr(), "flow");
    CompositeVector& s_l = *S->GetFieldData(saturation_liquid_key_, saturation_liquid_key_);
    CompositeVector& s_l_prev = *S->GetFieldData(prev_saturation_liquid_key_, passwd_);
    s_l_prev = s_l;

    S->GetFieldEvaluator("water_content")->HasFieldChanged(S.ptr(), "flow");
    CompositeVector& wc = *S->GetFieldData("water_content", "water_content");
    CompositeVector& wc_prev = *S->GetFieldData("prev_water_content", passwd_);
    wc_prev = wc;

    // We start with pressure equilibrium
    if (multiscale_porosity_) {
      *S->GetFieldData("pressure_matrix", passwd_)->ViewComponent("cell") =
          *S->GetFieldData(pressure_key_)->ViewComponent("cell");
      pressure_matrix_eval_->SetFieldAsChanged(S.ptr());
    }
  }

  // Trigger update of secondary fields depending on the primary pressure.
  pressure_eval_->SetFieldAsChanged(S.ptr());

  // Derive mass flux (state may not have it at time 0)
  double tmp;
  darcy_flux_copy->Norm2(&tmp);
  if (tmp == 0.0) {
    op_matrix_diff_->UpdateFlux(solution.ptr(), darcy_flux_copy.ptr());

    Epetra_MultiVector& flux = *darcy_flux_copy->ViewComponent("face", true);
    for (int f = 0; f < nfaces_owned; f++) flux[0][f] /= molar_rho_;
  }

  // Subspace entering: re-initialize lambdas.
  if (ti_list_->isSublist("pressure-lambda constraints") && solution->HasComponent("face")) {
    solver_name_constraint_ = ti_list_->sublist("pressure-lambda constraints").get<std::string>("linear solver");

    if (S->position() == Amanzi::TIME_PERIOD_START) {
      EnforceConstraints(t_ini, solution);
      pressure_eval_->SetFieldAsChanged(S.ptr());

      // update mass flux
      op_matrix_->Init();
      op_matrix_diff_->UpdateMatrices(Teuchos::null, solution.ptr());
      op_matrix_diff_->UpdateFlux(solution.ptr(), darcy_flux_copy.ptr());

      // normalize to Darcy flux, m/s
      Epetra_MultiVector& flux = *darcy_flux_copy->ViewComponent("face", true);
      for (int f = 0; f < nfaces_owned; f++) flux[0][f] /= molar_rho_;
    }
  }

  // Development: miscalleneous
  algebraic_water_content_balance_ = fp_list_->get<bool>("algebraic water content balance", false);
  if (algebraic_water_content_balance_) {
    CompositeVectorSpace cvs; 
    cvs.SetMesh(mesh_)->SetGhosted(false)
        ->AddComponent("cell", AmanziMesh::CELL, 1)
        ->AddComponent("dpre", AmanziMesh::CELL, 1);
    cnls_limiter_ = Teuchos::rcp(new CompositeVector(cvs));
  }

  // Verbose output of initialization statistics.
  InitializeStatistics_();
}


/* ****************************************************************
* This completes initialization of common fields that were not 
* initialized by the state.
**************************************************************** */
void Richards_PK::InitializeFields_()
{
  Teuchos::OSTab tab = vo_->getOSTab();

  // set popular default values for missed fields.
  // -- viscosity: if not initialized, we constant value from state.
  if (S_->GetField("viscosity_liquid")->owner() == passwd_) {
    double mu = *S_->GetScalarData("fluid_viscosity");

    if (!S_->GetField("viscosity_liquid", passwd_)->initialized()) {
      S_->GetFieldData("viscosity_liquid", passwd_)->PutScalar(mu);
      S_->GetField("viscosity_liquid", passwd_)->set_initialized();

      if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM)
          *vo_->os() << "initialized viscosity_liquid to input value " << mu << std::endl;  
    }
  }

  if (S_->GetField(saturation_liquid_key_)->owner() == passwd_) {
    if (S_->HasField(saturation_liquid_key_)) {
      if (!S_->GetField(saturation_liquid_key_, passwd_)->initialized()) {
        S_->GetFieldData(saturation_liquid_key_, passwd_)->PutScalar(1.0);
        S_->GetField(saturation_liquid_key_, passwd_)->set_initialized();

        if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM)
            *vo_->os() << "initialized saturation_liquid to default value 1.0" << std::endl;  
      }
    }
  }

  InitializeFieldFromField_(prev_saturation_liquid_key_, saturation_liquid_key_, true);
  InitializeFieldFromField_("prev_water_content", "water_content", true);

  // set matrix fields assuming presure equilibrium
  // -- pressure
  if (S_->HasField("pressure_matrix")) {
    // if (!S_->GetField("pressure_matrix", passwd_)->initialized()) {
      const Epetra_MultiVector& p1 = *S_->GetFieldData(pressure_key_)->ViewComponent("cell");
      Epetra_MultiVector& p0 = *S_->GetFieldData("pressure_matrix", passwd_)->ViewComponent("cell");
      p0 = p1;

      S_->GetField("pressure_matrix", passwd_)->set_initialized();
      pressure_matrix_eval_->SetFieldAsChanged(S_.ptr());

      if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM)
          *vo_->os() << "initialized pressure_matrix to pressure" << std::endl;  
    // }
  }

  // -- water contents 
  if (S_->HasField("water_content_matrix")) {
    if (!S_->GetField("water_content_matrix", passwd_)->initialized()) {
      CalculateVWContentMatrix_();
      S_->GetField("water_content_matrix", passwd_)->set_initialized();

      if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM)
          *vo_->os() << "initialized water_content_matrix to VWContent(pressure_matrix)" << std::endl;  
    }
  }

  InitializeFieldFromField_("prev_water_content_matrix", "water_content_matrix", false);
}


/* ****************************************************************
* Auxiliary initialization technique.
**************************************************************** */
void Richards_PK::InitializeFieldFromField_(
    const std::string& field0, const std::string& field1, bool call_evaluator)
{
  if (S_->HasField(field0)) {
    if (!S_->GetField(field0, passwd_)->initialized()) {
      if (call_evaluator)
          S_->GetFieldEvaluator(field1)->HasFieldChanged(S_.ptr(), passwd_);

      const CompositeVector& f1 = *S_->GetFieldData(field1);
      CompositeVector& f0 = *S_->GetFieldData(field0, passwd_);
      f0 = f1;

      S_->GetField(field0, passwd_)->set_initialized();

      if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM)
          *vo_->os() << "initialized " << field0 << " to " << field1 << std::endl;
    }
  }
}


/* ******************************************************************
* Print the header for new time period.
****************************************************************** */
void Richards_PK::InitializeStatistics_()
{
  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    std::string ti_method_name = ti_list_->get<std::string>("time integration method");
    std::string pc_name = ti_list_->get<std::string>("preconditioner");

    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os()<< "\nEC:" << error_control_ 
              << " Upwind:" << op_matrix_diff_->little_k()
              << " PC:\"" << pc_name.c_str() << "\"" 
              << " TI:\"" << ti_method_name.c_str() << "\"" << std::endl
              << "matrix: " << op_matrix_->PrintDiagnostics() << std::endl
              << "precon: " << op_preconditioner_->PrintDiagnostics() << std::endl;

    int missed_tmp = missed_bc_faces_;
    int dirichlet_tmp = dirichlet_bc_faces_;
#ifdef HAVE_MPI
    mesh_->get_comm()->SumAll(&missed_tmp, &missed_bc_faces_, 1);
    mesh_->get_comm()->SumAll(&dirichlet_tmp, &dirichlet_bc_faces_, 1);
#endif

    *vo_->os() << "pressure BC assigned to " << dirichlet_bc_faces_ << " faces" << std::endl;
    *vo_->os() << "default (no-flow) BC assigned to " << missed_bc_faces_ << " faces" << std::endl << std::endl;

    VV_PrintHeadExtrema(*solution);
    VV_PrintSourceExtrema();

    if (vo_->getVerbLevel() >= Teuchos::VERB_EXTREME) {
      *vo_->os() << "\nPrinting WRM files..." << std::endl;
      relperm_->PlotWRMcurves();
    }

    *vo_->os() << vo_->color("green") << "Initialization of PK is complete, T=" 
               << units_.OutputTime(S_->time()) << vo_->reset() << std::endl << std::endl;
  }
}


/* ******************************************************************* 
* Performs one time step from time t_old to time t_new either for
* steady-state or transient simulation. If reinit=true, enforce 
* p-lambda constraints.
******************************************************************* */
bool Richards_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  dt_ = t_new - t_old;

  // initialize statistics
  ms_itrs_ = 0;
  ms_calls_ = 0;

  // save a copy of primary and conservative fields
  // -- pressure
  CompositeVector pressure_copy(*S_->GetFieldData(pressure_key_, passwd_));

  // -- saturations, swap prev <- current
  S_->GetFieldEvaluator(saturation_liquid_key_)->HasFieldChanged(S_.ptr(), "flow");
  const CompositeVector& sat = *S_->GetFieldData(saturation_liquid_key_);
  CompositeVector& sat_prev = *S_->GetFieldData(prev_saturation_liquid_key_, passwd_);

  CompositeVector sat_prev_copy(sat_prev);
  sat_prev = sat;

  // -- water_conten, swap prev <- current
  S_->GetFieldEvaluator("water_content")->HasFieldChanged(S_.ptr(), "flow");
  CompositeVector& wc = *S_->GetFieldData("water_content", "water_content");
  CompositeVector& wc_prev = *S_->GetFieldData("prev_water_content", passwd_);

  CompositeVector wc_prev_copy(wc_prev);
  wc_prev = wc;

  // -- field for multiscale models, save and swap
  Teuchos::RCP<CompositeVector> pressure_matrix_copy, wc_matrix_prev_copy;
  if (multiscale_porosity_) {
    pressure_matrix_copy = Teuchos::rcp(new CompositeVector(*S_->GetFieldData("pressure_matrix", passwd_)));

    CompositeVector& wc_matrix = *S_->GetFieldData("water_content_matrix", passwd_);
    CompositeVector& wc_matrix_prev = *S_->GetFieldData("prev_water_content_matrix", passwd_);

    wc_matrix_prev_copy = Teuchos::rcp(new CompositeVector(wc_matrix_prev));
    wc_matrix_prev = wc_matrix;
  }

  // enter subspace
  if (reinit && solution->HasComponent("face")) {
    EnforceConstraints(t_new, solution);

    if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
      VV_PrintHeadExtrema(*solution);
    }
  }

  // initialization
  if (num_itrs_ == 0) {
    Teuchos::RCP<TreeVector> udot = Teuchos::rcp(new TreeVector(*soln_));
    udot->PutScalar(0.0);
    bdf1_dae->SetInitialState(t_old, soln_, udot);

    UpdatePreconditioner(t_old, soln_, dt_);
    num_itrs_++;
  }

  // trying to make a step
  bool failed(false);
  failed = bdf1_dae->TimeStep(dt_, dt_next_, soln_);
  if (failed) {
    dt_ = dt_next_;

    // revover the original primary solution, pressure
    *S_->GetFieldData(pressure_key_, passwd_) = pressure_copy;
    pressure_eval_->SetFieldAsChanged(S_.ptr());

    // revover the original fields
    *S_->GetFieldData(prev_saturation_liquid_key_, passwd_) = sat_prev_copy;
    *S_->GetFieldData("prev_water_content", passwd_) = wc_prev_copy;

    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "Reverted pressure, prev_saturation_liquid, prev_water_content" << std::endl;

    if (multiscale_porosity_) {
      *S_->GetFieldData("pressure_matrix", passwd_) = *pressure_matrix_copy;
      *S_->GetFieldData("prev_water_content_matrix", passwd_) = *wc_matrix_prev_copy;

      *vo_->os() << "Reverted pressure_matrix, prev_water_content_matrix" << std::endl;
    }

    return failed;
  }

  // commit solution (should we do it here ?)
  bdf1_dae->CommitSolution(dt_, soln_);
  pressure_eval_->SetFieldAsChanged(S_.ptr());

  dt_tuple times(t_old, dt_);
  dT_history_.push_back(times);
  num_itrs_++;

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    VV_ReportWaterBalance(S_.ptr());
    VV_ReportMultiscale();
  }
  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    VV_ReportSeepageOutflow(S_.ptr(), dt_);
  }

  dt_ = dt_next_;
  
  return failed;
}


/* ******************************************************************
* Save internal data needed by time integration. Calculate temporarily
* the Darcy flux.
****************************************************************** */
void Richards_PK::CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S)
{
  // calculate Darcy flux.
  Teuchos::RCP<CompositeVector> darcy_flux = S->GetFieldData(darcy_flux_key_, passwd_);
  op_matrix_diff_->UpdateFlux(solution.ptr(), darcy_flux.ptr());

  Epetra_MultiVector& flux = *darcy_flux->ViewComponent("face", true);
  for (int f = 0; f < nfaces_owned; f++) flux[0][f] /= molar_rho_;
  *darcy_flux_copy->ViewComponent("face", true) = flux;

  // update time derivative
  *pdot_cells_prev = *pdot_cells;

  dt_ = dt_next_;
}


/* ******************************************************************
* Returns either known pressure face value or calculates it using
* the two-point flux approximation (FV) scheme.
****************************************************************** */
double Richards_PK::BoundaryFaceValue(int f, const CompositeVector& u)
{
  double face_value;

  if (u.HasComponent("face")) {
    const Epetra_MultiVector& u_face = *u.ViewComponent("face");
    face_value = u_face[0][f];
  } else {
    int c = BoundaryFaceGetCell(f);
    face_value = DeriveBoundaryFaceValue(f, u, wrm_->second[(*wrm_->first)[c]]);
  }
  return face_value;
}


/* ******************************************************************
* Calculates pressure value on the boundary using the two-point flux 
* approximation (FV) scheme.
****************************************************************** */
double Richards_PK::DeriveBoundaryFaceValue(
    int f, const CompositeVector& u, Teuchos::RCP<const WRM> wrm_model) 
{
  const std::vector<int>& bc_model = op_bc_->bc_model();
  const std::vector<double>& bc_value = op_bc_->bc_value();

  if (bc_model[f] == Operators::OPERATOR_BC_DIRICHLET) {
    return bc_value[f];
  } else {
    const Epetra_MultiVector& mu_cell = *S_->GetFieldData("viscosity_liquid")->ViewComponent("cell");
    const Epetra_MultiVector& u_cell = *u.ViewComponent("cell");
    AmanziMesh::Entity_ID_List cells;
    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    int c = cells[0];

    double pc_shift(atm_pressure_);   
    double trans_f = op_matrix_diff_->ComputeTransmissibility(f);
    double g_f = op_matrix_diff_->ComputeGravityFlux(f);
    double lmd = u_cell[0][c];
    int dir;
    mesh_->face_normal(f, false, c, &dir);
    double bnd_flux = dir * bc_value[f] / (molar_rho_ / mu_cell[0][c]);

    double max_val(atm_pressure_), min_val;
    if (bnd_flux <= 0.0) {
      min_val = u_cell[0][c];
    } else {
      min_val= u_cell[0][c] + (g_f - bnd_flux) / (dir * trans_f);
    }
    double eps = std::max(1.0e-4 * std::abs(bnd_flux), 1.0e-8);

    // std::cout<<"min_val "<<min_val<<" max_val "<<max_val<<" "<<" trans_f "<<trans_f<<"\n";
    // std::cout<<"g_f "<<g_f<<" bnd "<< bnd_flux <<" dir "<<dir<<"\n";
    // std::cout<<c <<"norm "<<n<<"\n";

    const KRelFn func = &WRM::k_relative;       
    Amanzi::BoundaryFaceSolver<WRM> bnd_solver(trans_f, g_f, u_cell[0][c], lmd, bnd_flux, dir, pc_shift, 
                                               min_val, max_val, eps, wrm_model, func);
    lmd = bnd_solver.FaceValue();

    return lmd;      
  }
}


/* ******************************************************************
* This is strange.
****************************************************************** */
void Richards_PK::VV_ReportMultiscale()
{
  if (multiscale_porosity_ && ms_calls_ && 
      vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "multiscale: NS:" << double(ms_itrs_) / ms_calls_ << std::endl;
  }
}


/* ******************************************************************
* This is strange.
****************************************************************** */
void Richards_PK::CalculateDiagnostics(const Teuchos::RCP<State>& S) {
  UpdateLocalFields_(S.ptr());
}

}  // namespace Flow
}  // namespace Amanzi

