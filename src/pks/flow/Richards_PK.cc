/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Neil Carlson (nnc@lanl.gov), 
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <vector>

// TPLs
#include "boost/math/tools/roots.hpp"
#include "boost/algorithm/string.hpp"
#include "Epetra_IntVector.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

// Amanzi
#include "BoundaryFlux.hh"
#include "dbc.hh"
#include "exceptions.hh"
#include "independent_variable_field_evaluator_fromfunction.hh"
#include "Mesh.hh"
#include "mfd3d_diffusion.hh"
#include "OperatorDefs.hh"
#include "OperatorDiffusionFactory.hh"
#include "PK_Utils.hh"
#include "Point.hh"
#include "primary_variable_field_evaluator.hh"
#include "UpwindFactory.hh"
#include "XMLParameterListWriter.hh"

// Flow
#include "DarcyVelocityEvaluator.hh"
#include "Flow_BC_Factory.hh"
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
    Flow_PK(),
    glist_(glist),
    soln_(soln)
{
  S_ = S;

  std::string pk_name = pk_tree.name();

  boost::iterator_range<std::string::iterator> res = boost::algorithm::find_last(pk_name,"->"); 
  if (res.end() - pk_name.end() != 0) boost::algorithm::erase_head(pk_name,  res.end() - pk_name.begin());

  
  // We need the flow list
  Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(glist, "PKs", true);
  Teuchos::RCP<Teuchos::ParameterList> flow_list = Teuchos::sublist(pk_list, pk_name, true);
  rp_list_ = Teuchos::sublist(flow_list, "Richards problem", true);
  
  // We also need miscaleneous sublists
  preconditioner_list_ = Teuchos::sublist(glist, "Preconditioners", true);
  linear_operator_list_ = Teuchos::sublist(glist, "Solvers", true);
  ti_list_ = Teuchos::sublist(rp_list_, "time integrator");

  vo_ = NULL;
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
  Teuchos::RCP<Teuchos::ParameterList> flow_list = Teuchos::sublist(pk_list, pk_list_name, true);
  rp_list_ = Teuchos::sublist(flow_list, "Richards problem", true);
 
  // We also need miscaleneous sublists
  preconditioner_list_ = Teuchos::sublist(glist, "Preconditioners", true);
  linear_operator_list_ = Teuchos::sublist(glist, "Solvers", true);
  ti_list_ = Teuchos::sublist(rp_list_, "time integrator");

  ms_itrs_ = 0;
  ms_calls_ = 0;

  vo_ = NULL;
}


/* ******************************************************************
* Clean memory.
****************************************************************** */
Richards_PK::~Richards_PK()
{
  if (bc_pressure != NULL) delete bc_pressure;
  if (bc_flux != NULL) delete bc_flux;
  if (bc_head != NULL) delete bc_head;
  if (bc_seepage != NULL) delete bc_seepage;

  for (int i = 0; i < srcs.size(); i++) {
    if (srcs[i] != NULL) delete srcs[i]; 
  }
  if (vo_ != NULL) delete vo_;
}


/* ******************************************************************
* Define structure of this PK. We request physical fields and their
* evaluators. Selection of a few models is available and driven by
* model factories, evaluator factories, and parameters of the list
* "physical models and assumptions".
****************************************************************** */
void Richards_PK::Setup()
{
  dt_ = 0.0;
  mesh_ = S_->GetMesh();
  dim = mesh_->space_dimension();

  Flow_PK::Setup();

  // Our decision can be affected by the list of models
  Teuchos::RCP<Teuchos::ParameterList> physical_models =
      Teuchos::sublist(rp_list_, "physical models and assumptions");
  std::string vwc_model = physical_models->get<std::string>("water content model", "constant density");
  std::string multiscale_model = physical_models->get<std::string>("multiscale model", "single porosity");

  // Require primary field for this PK, which is pressure
  std::vector<std::string> names;
  std::vector<AmanziMesh::Entity_kind> locations;
  std::vector<int> ndofs;

  Teuchos::RCP<Teuchos::ParameterList> list1 = Teuchos::sublist(rp_list_, "operators", true);
  Teuchos::RCP<Teuchos::ParameterList> list2 = Teuchos::sublist(list1, "diffusion operator", true);
  Teuchos::RCP<Teuchos::ParameterList> list3 = Teuchos::sublist(list2, "matrix", true);
  std::string name = list3->get<std::string>("discretization primary");

  names.push_back("cell");
  locations.push_back(AmanziMesh::CELL);
  ndofs.push_back(1);
  if (name != "fv: default") {
    names.push_back("face");
    locations.push_back(AmanziMesh::FACE);
    ndofs.push_back(1);
  }

  if (!S_->HasField("pressure")) {
    S_->RequireField("pressure", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponents(names, locations, ndofs);

    Teuchos::ParameterList elist;
    elist.set<std::string>("evaluator name", "pressure");
    pressure_eval_ = Teuchos::rcp(new PrimaryVariableFieldEvaluator(elist));
    S_->SetFieldEvaluator("pressure", pressure_eval_);
  }

  // Require conserved quantity.
  // -- water content
  if (!S_->HasField("water_content")) {
    S_->RequireField("water_content", "water_content")->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);

    Teuchos::ParameterList vwc_list;
    VWContentEvaluatorFactory fac;
    Teuchos::RCP<VWContentEvaluator> eval = fac.Create(vwc_model, vwc_list);
    S_->SetFieldEvaluator("water_content", eval);
  }

  // -- water content from the previous time step
  if (!S_->HasField("prev_water_content")) {
    S_->RequireField("prev_water_content", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    S_->GetField("prev_water_content", passwd_)->set_io_vis(false);
  }

  // -- multiscale extension: secondary (immobile water content)
  if (multiscale_model == "dual porosity") {
    if (!S_->HasField("pressure_matrix")) {
      S_->RequireField("pressure_matrix", passwd_)->SetMesh(mesh_)->SetGhosted(false)
        ->SetComponent("cell", AmanziMesh::CELL, 1);

      Teuchos::ParameterList elist;
      elist.set<std::string>("evaluator name", "pressure_matrix");
      pressure_matrix_eval_ = Teuchos::rcp(new PrimaryVariableFieldEvaluator(elist));
      S_->SetFieldEvaluator("pressure_matrix", pressure_matrix_eval_);
    }

    Teuchos::RCP<Teuchos::ParameterList> msp_list = Teuchos::sublist(rp_list_, "multiscale models", true);
    msp_ = CreateMultiscaleFlowPorosityPartition(mesh_, msp_list);

    if (!S_->HasField("water_content_matrix")) {
      S_->RequireField("water_content_matrix", passwd_)->SetMesh(mesh_)->SetGhosted(true)
        ->SetComponent("cell", AmanziMesh::CELL, 1);
    }
    if (!S_->HasField("prev_water_content_matrix")) {
      S_->RequireField("prev_water_content_matrix", passwd_)->SetMesh(mesh_)->SetGhosted(true)
        ->SetComponent("cell", AmanziMesh::CELL, 1);
      S_->GetField("prev_water_content_matrix", passwd_)->set_io_vis(false);
    }

    S_->RequireField("porosity_matrix", "porosity_matrix")->SetMesh(mesh_)->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, 1);

    Teuchos::ParameterList elist;
    elist.set<std::string>("porosity key", "porosity_matrix");
    elist.set<std::string>("pressure key", "pressure_matrix");
    Teuchos::RCP<PorosityModelPartition> pom = CreatePorosityModelPartition(mesh_, msp_list);
    Teuchos::RCP<PorosityModelEvaluator> eval = Teuchos::rcp(new PorosityModelEvaluator(elist, pom));
    S_->SetFieldEvaluator("porosity_matrix", eval);
  }

  // Require additional fields and evaluators for this PK.
  // -- darcy flux
  if (!S_->HasField("darcy_flux")) {
    S_->RequireField("darcy_flux", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("face", AmanziMesh::FACE, 1);

    Teuchos::ParameterList elist;
    elist.set<std::string>("evaluator name", "darcy_flux");
    darcy_flux_eval_ = Teuchos::rcp(new PrimaryVariableFieldEvaluator(elist));
    S_->SetFieldEvaluator("darcy_flux", darcy_flux_eval_);
  }

  // -- porosity
  if (!S_->HasField("porosity")) {
    S_->RequireField("porosity", "porosity")->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);

    Teuchos::RCP<Teuchos::ParameterList> physical_models = 
        Teuchos::sublist(rp_list_, "physical models and assumptions");
    std::string pom_name = physical_models->get<std::string>("porosity model", "constant porosity");

    if (pom_name == "compressible: pressure function") {
      Teuchos::RCP<Teuchos::ParameterList>
          pom_list = Teuchos::sublist(rp_list_, "porosity models", true);
      Teuchos::RCP<PorosityModelPartition> pom = CreatePorosityModelPartition(mesh_, pom_list);

      Teuchos::ParameterList elist;
      // elist.sublist("VerboseObject").set<std::string>("Verbosity Level", "extreme");
      Teuchos::RCP<PorosityModelEvaluator> eval = Teuchos::rcp(new PorosityModelEvaluator(elist, pom));
      S_->SetFieldEvaluator("porosity", eval);
    } else {
      S_->RequireFieldEvaluator("porosity");
    }
  }

  // -- viscosity: if not requested by any PK, we request its constant value.
  if (!S_->HasField("viscosity_liquid")) {
    if (!S_->HasField("fluid_viscosity")) {
      S_->RequireScalar("fluid_viscosity", passwd_);
    }
    S_->RequireField("viscosity_liquid", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    S_->GetField("viscosity_liquid", passwd_)->set_io_vis(false);
  }

  // -- model for liquid density is constant density unless specified otherwise
  //    in high-level PKs.
  if (!S_->HasField("molar_density_liquid")) {
    S_->RequireField("molar_density_liquid", "molar_density_liquid")->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);

    double rho = glist_->sublist("State").sublist("initial conditions")
                        .sublist("fluid_density").get<double>("value", 1000.0);
    double n_l = rho / CommonDefs::MOLAR_MASS_H2O;

    Teuchos::ParameterList& wc_eval = S_->FEList().sublist("molar_density_liquid");
    wc_eval.sublist("function").sublist("DOMAIN")
           .set<std::string>("region", "All")
           .set<std::string>("component","cell")
           .sublist("function").sublist("function-constant")
           .set<double>("value", n_l);
    wc_eval.set<std::string>("field evaluator type", "independent variable");

    S_->RequireFieldEvaluator("molar_density_liquid");
  }
  
  // -- saturation
  double patm = rp_list_->get<double>("atmospheric pressure", FLOW_PRESSURE_ATMOSPHERIC);

  Teuchos::RCP<Teuchos::ParameterList>
      wrm_list = Teuchos::sublist(rp_list_, "water retention models", true);
  wrm_ = CreateWRMPartition(mesh_, wrm_list);

  if (!S_->HasField("saturation_liquid")) {
    S_->RequireField("saturation_liquid", "saturation_liquid")->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);

    Teuchos::ParameterList elist;
    // elist.sublist("VerboseObject").set<std::string>("Verbosity Level", "extreme");
    Teuchos::RCP<WRMEvaluator> eval = Teuchos::rcp(new WRMEvaluator(elist, patm, wrm_));
    S_->SetFieldEvaluator("saturation_liquid", eval);
  }

  if (!S_->HasField("prev_saturation_liquid")) {
    S_->RequireField("prev_saturation_liquid", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    S_->GetField("prev_saturation_liquid", passwd_)->set_io_vis(false);
  }

  // Local fields and evaluators.
  // -- hydraulic head
  if (!S_->HasField("hydraulic_head")) {
    S_->RequireField("hydraulic_head", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  // -- Darcy velocity vector
  S_->RequireField("darcy_velocity", "darcy_velocity")->SetMesh(mesh_)->SetGhosted(true)
    ->SetComponent("cell", AmanziMesh::CELL, dim);

  Teuchos::ParameterList elist;
  Teuchos::RCP<DarcyVelocityEvaluator> eval = Teuchos::rcp(new DarcyVelocityEvaluator(elist));
  S_->SetFieldEvaluator("darcy_velocity", eval);
}


/* ******************************************************************
* This is long but simple subroutine. It goes through time integrator
* list and initializes various objects created during setup step.
* Some local objects needs to the most 
****************************************************************** */
void Richards_PK::Initialize()
{
  // times, initialization could be done on any non-zero interval.
  double t_old = S_->time(); 
  dt_ = ti_list_->get<double>("initial time step", 1.0);
  double t_new = t_old + dt_;

  dt_desirable_ = dt_;
  dt_next_ = dt_;

  // Initialize miscalleneous default parameters.
  error_control_ = FLOW_TI_ERROR_CONTROL_PRESSURE;

  // create verbosity object
  Teuchos::ParameterList vlist;
  vlist.sublist("VerboseObject") = rp_list_->sublist("VerboseObject");
  vo_ = new VerboseObject("FlowPK::Richards", vlist); 

  // Initilize various common data depending on mesh and state.
  Flow_PK::Initialize();

  // Create local evaluators. Initialize local fields.
  InitializeFields_();
  UpdateLocalFields_();

  // Initialize BCs and source terms.
  InitializeBCsSources_(*rp_list_);
  op_bc_ = Teuchos::rcp(new Operators::BCs(Operators::OPERATOR_BC_TYPE_FACE, bc_model, bc_value, bc_mixed));

  // Create relative permeability
  Teuchos::RCP<Teuchos::ParameterList> upw_list = Teuchos::sublist(rp_list_, "upwind", true);
  relperm_ = Teuchos::rcp(new RelPerm(*upw_list, mesh_, atm_pressure_, wrm_));

  CompositeVectorSpace cvs; 
  cvs.SetMesh(mesh_)->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::CELL, 1)
      ->AddComponent("face", AmanziMesh::FACE, 1);

  krel_ = Teuchos::rcp(new CompositeVector(cvs));
  dKdP_ = Teuchos::rcp(new CompositeVector(cvs));

  krel_upwind_method_ = FLOW_RELATIVE_PERM_NONE;
  krel_->PutScalarMasterAndGhosted(1.0);
  dKdP_->PutScalarMasterAndGhosted(0.0);

  // parameter which defines when update direction of update
  Operators::UpwindFactory<RelPerm> upwind_factory;
  upwind_ = upwind_factory.Create(mesh_, relperm_, *upw_list);

  std::string upw_upd = upw_list->get<std::string>("upwind update", "every timestep");
  if (upw_upd == "every nonlinear iteration") update_upwind = FLOW_UPWIND_UPDATE_ITERATION;
  else update_upwind = FLOW_UPWIND_UPDATE_TIMESTEP;  

  // models and assumptions
  // -- coupling with other physical PKs
  Teuchos::RCP<Teuchos::ParameterList> physical_models = 
      Teuchos::sublist(rp_list_, "physical models and assumptions");
  vapor_diffusion_ = physical_models->get<bool>("vapor diffusion", false);
  multiscale_porosity_ = (physical_models->get<std::string>(
      "multiscale model", "single porosity") != "single porosity");

  // Initialize actions on boundary condtions. 
  flux_units_ = molar_rho_ / rho_;

  ProcessShiftWaterTableList(*rp_list_);

  bc_pressure->Compute(t_new);
  bc_flux->Compute(t_new);
  bc_seepage->Compute(t_new);
  if (shift_water_table_.getRawPtr() == NULL)
    bc_head->Compute(t_new);
  else
    bc_head->ComputeShift(t_new, shift_water_table_->Values());

  // Process other fundamental structures.
  SetAbsolutePermeabilityTensor();

  // Select a proper matrix class. 
  const Teuchos::ParameterList& tmp_list = rp_list_->sublist("operators")
                                                    .sublist("diffusion operator");
  Teuchos::ParameterList oplist_matrix = tmp_list.sublist("matrix");
  Teuchos::ParameterList oplist_pc = tmp_list.sublist("preconditioner");

  std::string name = rp_list_->sublist("upwind").get<std::string>("relative permeability");
  std::string upw_method("standard: cell");
  if (name == "upwind: darcy velocity") {
    upw_method = "upwind: face";
  } else if (name == "upwind: gravity") {
    upw_method = "upwind: face";
  } else if (name == "upwind: amanzi") {
    upw_method = "divk: cell-face";
    // upw_method = "divk: face";
  } else if (name == "other: arithmetic average") {
    upw_method = "upwind: face";
  }
  oplist_matrix.set<std::string>("nonlinear coefficient", upw_method);
  oplist_pc.set<std::string>("nonlinear coefficient", upw_method);

  Operators::OperatorDiffusionFactory opfactory;
  op_matrix_diff_ = opfactory.Create(oplist_matrix, mesh_, op_bc_, rho_, gravity_);
  op_matrix_ = op_matrix_diff_->global_operator();
  op_preconditioner_diff_ = opfactory.Create(oplist_pc, mesh_, op_bc_, rho_, gravity_);
  op_preconditioner_ = op_preconditioner_diff_->global_operator();
  op_acc_ = Teuchos::rcp(new Operators::OperatorAccumulation(AmanziMesh::CELL, op_preconditioner_));

  if (vapor_diffusion_) {
    Teuchos::ParameterList oplist_vapor = tmp_list.sublist("vapor matrix");
    op_vapor_diff_ = opfactory.Create(oplist_vapor, mesh_, op_bc_);
    op_vapor_ = op_vapor_diff_->global_operator();
    op_preconditioner_->OpPushBack(op_vapor_diff_->local_matrices(),
                                   Operators::OPERATOR_PROPERTY_DATA_READ_ONLY);
  }

  // Create pointers to the primary flow field pressure.
  solution = S_->GetFieldData("pressure", passwd_);
  soln_->SetData(solution); 
  
  // Create auxiliary vectors for time history and error estimates.
  const Epetra_BlockMap& cmap_owned = mesh_->cell_map(false);
  pdot_cells_prev = Teuchos::rcp(new Epetra_Vector(cmap_owned));
  pdot_cells = Teuchos::rcp(new Epetra_Vector(cmap_owned));

  // Initialize two fields for upwind operators.
  darcy_flux_copy = Teuchos::rcp(new CompositeVector(*S_->GetFieldData("darcy_flux", passwd_)));
  InitializeUpwind_();

  // Other quantatities: injected water mass
  mass_bc = 0.0;
  seepage_mass_ = 0.0;
  mass_initial = 0.0;
  initialize_with_darcy_ = true;
  num_itrs_ = 0;
  
  // Conditional initialization of lambdas from pressures.
  CompositeVector& pressure = *S_->GetFieldData("pressure", passwd_);

  if (pressure.HasComponent("face")) {
    DeriveFaceValuesFromCellValues(*pressure.ViewComponent("cell"),
                                   *pressure.ViewComponent("face"));
  }

  // error control options
  ASSERT(ti_list_->isParameter("error control options"));

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
  ASSERT(ti_method_name == "BDF1");
  Teuchos::ParameterList& bdf1_list = ti_list_->sublist("BDF1");

  if (! bdf1_list.isSublist("VerboseObject"))
      bdf1_list.sublist("VerboseObject") = rp_list_->sublist("VerboseObject");

  bdf1_dae = Teuchos::rcp(new BDF1_TI<TreeVector, TreeVectorSpace>(*this, bdf1_list, soln_));

  // complete other steps
  // repeat upwind initialization, mainly for old MPC
  InitializeUpwind_();

  // initialize well modeling
  for (int i = 0; i < srcs.size(); ++i) {
    int type = srcs[i]->CollectActionsList();
    if (type & CommonDefs::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_PERMEABILITY) {
      PKUtils_CalculatePermeabilityFactorInWell(S_, Kxy);
    }
    srcs[i]->Compute(t_old, t_new, Kxy); 
  }

  // initialize matrix and preconditioner operators.
  // -- setup phase
  // -- molar density requires to rescale gravity later.
  op_matrix_->Init();
  Teuchos::RCP<std::vector<WhetStone::Tensor> > Kptr = Teuchos::rcpFromRef(K);
  op_matrix_diff_->SetBCs(op_bc_, op_bc_);
  op_matrix_diff_->Setup(Kptr, krel_, dKdP_);

  op_preconditioner_->Init();
  op_preconditioner_->SetBCs(op_bc_, op_bc_);
  op_preconditioner_diff_->Setup(Kptr, krel_, dKdP_);

  // -- assemble phase
  UpdateSourceBoundaryData(t_old, t_new, pressure);

  op_matrix_diff_->UpdateMatrices(Teuchos::null, solution.ptr());
  op_matrix_diff_->ApplyBCs(true, true);

  op_preconditioner_diff_->UpdateMatrices(darcy_flux_copy.ptr(), solution.ptr());
  op_preconditioner_diff_->UpdateMatricesNewtonCorrection(darcy_flux_copy.ptr(), solution.ptr(), molar_rho_);
  op_preconditioner_diff_->ApplyBCs(true, true);
  op_preconditioner_->SymbolicAssembleMatrix();

  if (vapor_diffusion_) {
    // op_vapor_diff_->SetBCs(op_bc_);
    op_vapor_diff_->SetScalarCoefficient(Teuchos::null, Teuchos::null);
  }

  // generic linear solver for all cases except for a few
  ASSERT(ti_list_->isParameter("linear solver"));
  solver_name_ = ti_list_->get<std::string>("linear solver");

  // preconditioner or encapsulated preconditioner
  ASSERT(ti_list_->isParameter("preconditioner"));
  preconditioner_name_ = ti_list_->get<std::string>("preconditioner");
  ASSERT(preconditioner_list_->isSublist(preconditioner_name_));
  
  op_pc_solver_ = op_preconditioner_;

  if (ti_list_->isParameter("preconditioner enhancement")) {
    std::string tmp_solver = ti_list_->get<std::string>("preconditioner enhancement");
    if (tmp_solver != "none") {
      ASSERT(linear_operator_list_->isSublist(tmp_solver));

      AmanziSolvers::LinearOperatorFactory<Operators::Operator, CompositeVector, CompositeVectorSpace> sfactory;
      op_pc_solver_ = sfactory.Create(tmp_solver, *linear_operator_list_, op_preconditioner_);
    }
  }
  
  // Optional step: calculate hydrostatic solution consistent with BCs
  // and clip it as requested. We have to do it only once at the beginning
  // of time period.
  if (ti_list_->isSublist("initialization") && initialize_with_darcy_ 
      && S_->position() == Amanzi::TIME_PERIOD_START) {
    initialize_with_darcy_ = false;
    Teuchos::ParameterList& ini_list = ti_list_->sublist("initialization");
 
    std::string solver_name_ini = ini_list.get<std::string>("method", "none");
    if (solver_name_ini == "saturated solver") {
      std::string name = ini_list.get<std::string>("linear solver");
      SolveFullySaturatedProblem(t_old, *solution, name);

      bool clip(false);
      double clip_saturation = ini_list.get<double>("clipping saturation value", -1.0);
      if (clip_saturation > 0.0) {
        double pmin = FLOW_PRESSURE_ATMOSPHERIC;
        Epetra_MultiVector& p = *solution->ViewComponent("cell");
        ClipHydrostaticPressure(pmin, clip_saturation, p);
        clip = true;
      }

      double clip_pressure = ini_list.get<double>("clipping pressure value", -1e+10);
      if (clip_pressure > -5 * FLOW_PRESSURE_ATMOSPHERIC) {
        Epetra_MultiVector& p = *solution->ViewComponent("cell");
        ClipHydrostaticPressure(clip_pressure, p);
        clip = true;
      }

      if (clip && vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
        Teuchos::OSTab tab = vo_->getOSTab();
        *vo_->os() << "Clipped pressure field." << std::endl;
      }

      if (clip && solution->HasComponent("face")) {
        Epetra_MultiVector& p = *solution->ViewComponent("cell");
        Epetra_MultiVector& lambda = *solution->ViewComponent("face", true);
        DeriveFaceValuesFromCellValues(p, lambda);
      }
    }
    else if (solver_name_ini == "picard") {
      AdvanceToSteadyState_Picard(ti_list_->sublist("initialization").sublist("picard parameters"));
    }
    pressure_eval_->SetFieldAsChanged(S_.ptr());

    // initialization is usually done at time 0, so we need to update other
    // fields such as prev_saturation_liquid
    S_->GetFieldEvaluator("saturation_liquid")->HasFieldChanged(S_.ptr(), "flow");
    CompositeVector& s_l = *S_->GetFieldData("saturation_liquid", "saturation_liquid");
    CompositeVector& s_l_prev = *S_->GetFieldData("prev_saturation_liquid", passwd_);
    s_l_prev = s_l;

    S_->GetFieldEvaluator("water_content")->HasFieldChanged(S_.ptr(), "flow");
    CompositeVector& wc = *S_->GetFieldData("water_content", "water_content");
    CompositeVector& wc_prev = *S_->GetFieldData("prev_water_content", passwd_);
    wc_prev = wc;

    // We start with pressure equilibrium
    if (multiscale_porosity_) {
      *S_->GetFieldData("pressure_matrix", passwd_)->ViewComponent("cell") =
          *S_->GetFieldData("pressure")->ViewComponent("cell");
      pressure_matrix_eval_->SetFieldAsChanged(S_.ptr());
    }
  }

  // Trigger update of secondary fields depending on the primary pressure.
  pressure_eval_->SetFieldAsChanged(S_.ptr());

  // derive mass flux (state may not have it at time 0)
  double tmp;
  darcy_flux_copy->Norm2(&tmp);
  if (tmp == 0.0) {
    op_matrix_diff_->UpdateFlux(*solution, *darcy_flux_copy);

    Epetra_MultiVector& flux = *darcy_flux_copy->ViewComponent("face", true);
    for (int f = 0; f < nfaces_owned; f++) flux[0][f] /= molar_rho_;
  }

  // subspace entering: re-initialize lambdas.
  if (ti_list_->isSublist("pressure-lambda constraints") && solution->HasComponent("face")) {
    solver_name_constraint_ = ti_list_->sublist("pressure-lambda constraints").get<std::string>("linear solver");

    if (S_->position() == Amanzi::TIME_PERIOD_START) {
      EnforceConstraints(t_new, solution);
      pressure_eval_->SetFieldAsChanged(S_.ptr());

      // update mass flux
      op_matrix_->Init();
      op_matrix_diff_->UpdateMatrices(Teuchos::null, solution.ptr());
      op_matrix_diff_->UpdateFlux(*solution, *darcy_flux_copy);

      // normalize to Darcy flux, m/s
      Epetra_MultiVector& flux = *darcy_flux_copy->ViewComponent("face", true);
      for (int f = 0; f < nfaces_owned; f++) flux[0][f] /= molar_rho_;

      InitializeUpwind_();
    }
  }

  // miscalleneous
  algebraic_water_content_balance_ = rp_list_->get<bool>("algebraic water content balance", false);
  if (algebraic_water_content_balance_) {
    CompositeVectorSpace cvs; 
    cvs.SetMesh(mesh_)->SetGhosted(false)
        ->AddComponent("cell", AmanziMesh::CELL, 1)
        ->AddComponent("dpre", AmanziMesh::CELL, 1);
    cnls_limiter_ = Teuchos::rcp(new CompositeVector(cvs));
  }

  // verbose output
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
          *vo_->os() << "initilized viscosity_liquid to input value " << mu << std::endl;  
    }
  }

  if (S_->GetField("saturation_liquid")->owner() == passwd_) {
    if (S_->HasField("saturation_liquid")) {
      if (!S_->GetField("saturation_liquid", passwd_)->initialized()) {
        S_->GetFieldData("saturation_liquid", passwd_)->PutScalar(1.0);
        S_->GetField("saturation_liquid", passwd_)->set_initialized();

        if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM)
            *vo_->os() << "initiliazed saturation_liquid to default value 1.0" << std::endl;  
      }
    }
  }

  InitializeFieldFromField_("prev_saturation_liquid", "saturation_liquid", true);
  InitializeFieldFromField_("prev_water_content", "water_content", true);

  // set matrix fields assuming presure equilibrium
  // -- pressure
  if (S_->HasField("pressure_matrix")) {
    // if (!S_->GetField("pressure_matrix", passwd_)->initialized()) {
      const Epetra_MultiVector& p1 = *S_->GetFieldData("pressure")->ViewComponent("cell");
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
          *vo_->os() << "initiliazed " << field0 << " to " << field1 << std::endl;
    }
  }
}


/* ******************************************************************
* Set defaults parameters. It could be called only once.
****************************************************************** */
void Richards_PK::InitializeUpwind_()
{
  // Create RCP pointer to upwind flux.
  if (relperm_->method() == FLOW_RELATIVE_PERM_UPWIND_DARCY_FLUX ||
      relperm_->method() == FLOW_RELATIVE_PERM_AMANZI_MFD) {
    darcy_flux_upwind = darcy_flux_copy;
  } else if (relperm_->method() == FLOW_RELATIVE_PERM_UPWIND_GRAVITY) {
    darcy_flux_upwind = Teuchos::rcp(new CompositeVector(*darcy_flux_copy));
    relperm_->ComputeGravityFlux(K, gravity_, darcy_flux_upwind);
  } else {
    darcy_flux_upwind = Teuchos::rcp(new CompositeVector(*darcy_flux_copy));
    darcy_flux_upwind->PutScalar(0.0);
  }
}


/* ******************************************************************
* Print the header for new time period.
****************************************************************** */
void Richards_PK::InitializeStatistics_()
{
  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    std::string ti_method_name = ti_list_->get<std::string>("time integration method");

    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << std::endl 
        << vo_->color("green") << "Initalization of PK is complete, T=" << S_->time()
        << " dT=" << dt_ << vo_->reset() << std::endl;
    *vo_->os()<< "EC:" << error_control_ 
              << " Upwind:" << relperm_->method() << op_matrix_diff_->little_k()
              << " PC:\"" << preconditioner_name_.c_str() << "\"" 
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
  }
}


/* ******************************************************************* 
* Performs one time step from time t_old to time t_new either for
* steady-state or transient sumulation. If reinit=true, enforce 
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
  CompositeVector pressure_copy(*S_->GetFieldData("pressure", passwd_));

  // -- saturations, swap prev <- current
  S_->GetFieldEvaluator("saturation_liquid")->HasFieldChanged(S_.ptr(), "flow");
  const CompositeVector& sat = *S_->GetFieldData("saturation_liquid");
  CompositeVector& sat_prev = *S_->GetFieldData("prev_saturation_liquid", passwd_);

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
    *S_->GetFieldData("pressure", passwd_) = pressure_copy;
    pressure_eval_->SetFieldAsChanged(S_.ptr());

    // revover the original fields
    *S_->GetFieldData("prev_saturation_liquid", passwd_) = sat_prev_copy;
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
    VV_ReportSeepageOutflow(S_.ptr());
  }

  dt_ = dt_next_;
  
  return failed;
}


/* ******************************************************************
* Save internal data needed by time integration. Calculate temporarily
* the Darcy flux.
****************************************************************** */
void Richards_PK::CommitStep(double t_old, double t_new)
{
  // calculate Darcy flux.
  CompositeVector& darcy_flux = *S_->GetFieldData("darcy_flux", passwd_);
  op_matrix_diff_->UpdateFlux(*solution, darcy_flux);

  Epetra_MultiVector& flux = *darcy_flux.ViewComponent("face", true);
  for (int f = 0; f < nfaces_owned; f++) flux[0][f] /= molar_rho_;
  *darcy_flux_copy->ViewComponent("face", true) = flux;

  // update time derivative
  *pdot_cells_prev = *pdot_cells;

  dt_ = dt_next_;
}


/* ******************************************************************
* A wrapper for updating boundary conditions.
****************************************************************** */
void Richards_PK::UpdateSourceBoundaryData(double t_old, double t_new, const CompositeVector& u)
{
  for (int i = 0; i < srcs.size(); ++i) {
    srcs[i]->Compute(t_old, t_new, Kxy); 
  }

  bc_pressure->Compute(t_new);
  bc_flux->Compute(t_new);
  bc_seepage->Compute(t_new);
  if (shift_water_table_.getRawPtr() == NULL)
    bc_head->Compute(t_new);
  else
    bc_head->ComputeShift(t_new, shift_water_table_->Values());

  ComputeBCs(u);
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
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
    int c = cells[0];

    double pc_shift = atm_pressure_;   
    double trans_f = op_matrix_diff_->ComputeTransmissibility(f);
    double g_f = op_matrix_diff_->ComputeGravityFlux(f);
    double lmd = u_cell[0][c];
    int dir;
    const AmanziGeometry::Point n = mesh_->face_normal(f, false, c, &dir);
    double bnd_flux = dir*bc_value[f] / (molar_rho_ / mu_cell[0][c]);

    double max_val = atm_pressure_;
    double min_val;
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
void Richards_PK::CalculateDiagnostics() {
  UpdateLocalFields_();
}

}  // namespace Flow
}  // namespace Amanzi

