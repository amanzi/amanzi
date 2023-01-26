/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
           Daniil Svyatskiy
*/

/*
  This is the mpc_pk component of the Amanzi code.

  Process kernel that couples flow and energy in matrix and fractures.
*/

#include "FlowEnergy_PK.hh"
#include "InverseFactory.hh"
#include "PDE_CouplingFlux.hh"
#include "PDE_DiffusionFracturedMatrix.hh"
#include "TreeOperator.hh"
#include "SuperMap.hh"
#include "UniqueLocalIndex.hh"

#include "FlowEnergyMatrixFracture_PK.hh"
#include "FractureInsertion.hh"
#include "PK_MPCStrong.hh"
#include "PK_Physical.hh"

namespace Amanzi {

using CV_t = CompositeVector;
using CVS_t = CompositeVectorSpace;

/* *******************************************************************
* Constructor
******************************************************************* */
FlowEnergyMatrixFracture_PK::FlowEnergyMatrixFracture_PK(
  Teuchos::ParameterList& pk_tree,
  const Teuchos::RCP<Teuchos::ParameterList>& glist,
  const Teuchos::RCP<State>& S,
  const Teuchos::RCP<TreeVector>& soln)
  : Amanzi::PK_MPC<PK_BDF>(pk_tree, glist, S, soln),
    Amanzi::PK_MPCStrong<PK_BDF>(pk_tree, glist, S, soln),
    glist_(glist)
{
  Teuchos::ParameterList vlist;
  vo_ = Teuchos::rcp(new VerboseObject("CoupledThermalFlow_PK", vlist));
  Teuchos::RCP<Teuchos::ParameterList> pks_list = Teuchos::sublist(glist, "PKs");
  if (pks_list->isSublist(name_)) {
    plist_ = Teuchos::sublist(pks_list, name_);
  } else {
    std::stringstream messagestream;
    messagestream << "There is no sublist for PK " << name_ << "in PKs list\n";
    Errors::Message message(messagestream.str());
    Exceptions::amanzi_throw(message);
  }

  preconditioner_list_ = Teuchos::sublist(glist, "preconditioners", true);
  linear_operator_list_ = Teuchos::sublist(glist, "solvers", true);
  ti_list_ = Teuchos::sublist(plist_, "time integrator", true);
}


/* *******************************************************************
* Physics-based setup of PK.
******************************************************************* */
void
FlowEnergyMatrixFracture_PK::Setup()
{
  mesh_domain_ = S_->GetMesh();
  mesh_fracture_ = S_->GetMesh("fracture");

  // keys
  normal_permeability_key_ = "fracture-normal_permeability";
  normal_conductivity_key_ = "fracture-normal_conductivity";

  matrix_vol_flowrate_key_ = "volumetric_flow_rate";
  fracture_vol_flowrate_key_ = "fracture-volumetric_flow_rate";

  // primary and secondary fields for matrix affected by non-uniform
  // distribution of DOFs, so we need to define them here
  // -- pressure
  auto cvs = Operators::CreateFracturedMatrixCVS(mesh_domain_, mesh_fracture_);
  if (!S_->HasRecord("pressure")) {
    *S_->Require<CV_t, CVS_t>("pressure", Tags::DEFAULT).SetMesh(mesh_domain_)->SetGhosted(true) =
      *cvs;
    AddDefaultPrimaryEvaluator_("pressure", Tags::DEFAULT);
  }

  if (!S_->HasRecord("temperature")) {
    *S_->Require<CV_t, CVS_t>("temperature", Tags::DEFAULT)
       .SetMesh(mesh_domain_)
       ->SetGhosted(true) = *cvs;
    AddDefaultPrimaryEvaluator_("temperature", Tags::DEFAULT);
  }

  // -- darcy flux
  if (!S_->HasRecord(matrix_vol_flowrate_key_)) {
    std::string name("face");
    auto mmap = cvs->Map("face", false);
    auto gmap = cvs->Map("face", true);
    S_->Require<CV_t, CVS_t>(matrix_vol_flowrate_key_, Tags::DEFAULT)
      .SetMesh(mesh_domain_)
      ->SetGhosted(true)
      ->SetComponent(name, AmanziMesh::Entity_kind::FACE, mmap, gmap, 1);
    AddDefaultPrimaryEvaluator_(matrix_vol_flowrate_key_, Tags::DEFAULT);
  }

  // -- darcy flux for fracture
  if (!S_->HasRecord(fracture_vol_flowrate_key_)) {
    auto cvs2 = Operators::CreateManifoldCVS(mesh_fracture_);
    *S_->Require<CV_t, CVS_t>(fracture_vol_flowrate_key_, Tags::DEFAULT)
       .SetMesh(mesh_fracture_)
       ->SetGhosted(true) = *cvs2;
  }

  // Require additional fields and evaluators
  if (!S_->HasRecord(normal_permeability_key_)) {
    S_->Require<CV_t, CVS_t>(normal_permeability_key_, Tags::DEFAULT)
      .SetMesh(mesh_fracture_)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  }

  if (!S_->HasRecord(normal_conductivity_key_)) {
    S_->Require<CV_t, CVS_t>(normal_conductivity_key_, Tags::DEFAULT, "state")
      .SetMesh(mesh_fracture_)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  }

  // inform dependent PKs about coupling
  // -- flow
  auto pks0 = plist_->get<Teuchos::Array<std::string>>("PKs order").toVector();
  auto pks1 = glist_->sublist("PKs")
                .sublist(pks0[0])
                .get<Teuchos::Array<std::string>>("PKs order")
                .toVector();
  auto pks2 = glist_->sublist("PKs")
                .sublist(pks0[1])
                .get<Teuchos::Array<std::string>>("PKs order")
                .toVector();

  auto& mflow = glist_->sublist("PKs").sublist(pks1[0]).sublist("physical models and assumptions");
  mflow.set<std::string>("coupled matrix fracture flow", "matrix");

  auto& fflow = glist_->sublist("PKs").sublist(pks2[0]).sublist("physical models and assumptions");
  fflow.set<std::string>("coupled matrix fracture flow", "fracture");

  // -- energy
  auto& menergy =
    glist_->sublist("PKs").sublist(pks1[1]).sublist("physical models and assumptions");
  menergy.set<std::string>("coupled matrix fracture energy", "matrix");

  auto& fenergy =
    glist_->sublist("PKs").sublist(pks2[1]).sublist("physical models and assumptions");
  fenergy.set<std::string>("coupled matrix fracture energy", "fracture");

  // process other PKs
  PK_MPCStrong<PK_BDF>::Setup();
}


/* *******************************************************************
* Initialization create a tree operator to assemble global matrix
******************************************************************* */
void
FlowEnergyMatrixFracture_PK::Initialize()
{
  PK_MPCStrong<PK_BDF>::Initialize();

  // diagonal blocks (0,0) and (2,2) in tree operator must be Darcy PKs
  // one reason is that Darcy_PK combines matrix and preconditioner
  for (int i = 0; i < 2; ++i) {
    std::string pk_name = (*Teuchos::rcp_dynamic_cast<FlowEnergy_PK>(sub_pks_[i])->begin())->name();
    AMANZI_ASSERT(pk_name == "darcy" || pk_name == "richards");
  }

  // diagonal blocks in the 4x4 tree operator are the FlowEnergy PKs.
  // These blocks go into the preconditioner (op_tree_pc). The coupling
  // terms are linear and go into both the preconditioner and matrix ops.
  auto pk_matrix = Teuchos::rcp_dynamic_cast<FlowEnergy_PK>(sub_pks_[0]);
  auto pk_fracture = Teuchos::rcp_dynamic_cast<FlowEnergy_PK>(sub_pks_[1]);

  //
  // NOTE: In this PK, the solution is 2-deep, with 2 sub-PKs (matrix,
  // fracture) each with 2 sub-pks of their own (flow,energy).  Currently
  // TreeOperator cannot handle this, so instead we must flatten the map.
  auto tvs = Teuchos::rcp(new TreeVectorSpace(solution_->Map()));
  op_tree_matrix_ = Teuchos::rcp(new Operators::TreeOperator(tvs));
  for (int i = 0; i < 2; ++i) {
    const auto& row = solution_->SubVector(i)->get_map();
    for (int j = 0; j < 2; ++j) {
      const auto& col = solution_->SubVector(j)->get_map();
      op_tree_matrix_->set_block(i, j, Teuchos::rcp(new Operators::TreeOperator(row, col)));
    }
  }

  op_tree_pc_ = Teuchos::rcp(new Operators::TreeOperator(tvs));
  op_tree_pc_->set_block(0, 0, pk_matrix->op_tree_pc()->Clone());
  op_tree_pc_->set_block(1, 1, pk_fracture->op_tree_pc()->Clone());
  for (int i = 0; i < 2; ++i) {
    const auto& row = solution_->SubVector(i)->get_map();
    for (int j = 0; j < 2; ++j) {
      const auto& col = solution_->SubVector(j)->get_map();
      if (i != j) op_tree_pc_->set_block(i, j, Teuchos::rcp(new Operators::TreeOperator(row, col)));
    }
  }

  // off-diagonal blocks are coupled PDEs
  // -- minimum composite vector spaces containing the coupling term
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_matrix = S_->GetMesh("domain");
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_fracture = S_->GetMesh("fracture");

  auto& mmap = solution_->SubVector(0)->SubVector(0)->Data()->ViewComponent("face", false)->Map();
  auto& gmap = solution_->SubVector(0)->SubVector(0)->Data()->ViewComponent("face", true)->Map();

  // -- indices transmissibimility coefficients for matrix-fracture flux
  const auto& kn = *S_->Get<CV_t>(normal_permeability_key_).ViewComponent("cell");
  const auto& tn = *S_->Get<CV_t>(normal_conductivity_key_).ViewComponent("cell");
  double gravity = norm(S_->Get<AmanziGeometry::Point>("gravity"));

  FractureInsertion fi(mesh_matrix, mesh_fracture);
  fi.InitMatrixFaceToFractureCell(Teuchos::rcpFromRef(mmap), Teuchos::rcpFromRef(gmap));

  // -- operators (flow)
  fi.SetValues(kn, 1.0 / gravity);
  AddCouplingFluxes_(fi.get_cvs_matrix(),
                     fi.get_cvs_fracture(),
                     fi.get_inds_matrix(),
                     fi.get_inds_fracture(),
                     fi.get_values(),
                     0,
                     op_tree_matrix_);

  // diagonal coupling terms were added in the previous call
  AddCouplingFluxes_(fi.get_cvs_matrix(),
                     fi.get_cvs_fracture(),
                     fi.get_inds_matrix(),
                     fi.get_inds_fracture(),
                     fi.get_values(),
                     0,
                     op_tree_pc_);

  // -- operators (energy)
  fi.SetValues(tn, 1.0);
  AddCouplingFluxes_(fi.get_cvs_matrix(),
                     fi.get_cvs_fracture(),
                     fi.get_inds_matrix(),
                     fi.get_inds_fracture(),
                     fi.get_values(),
                     1,
                     op_tree_matrix_);

  AddCouplingFluxes_(fi.get_cvs_matrix(),
                     fi.get_cvs_fracture(),
                     fi.get_inds_matrix(),
                     fi.get_inds_fracture(),
                     fi.get_values(),
                     1,
                     op_tree_pc_);

  // -- indices transmissibimility coefficients for matrix-fracture advective flux
  FractureInsertion fia(mesh_matrix, mesh_fracture);
  fia.InitMatrixCellToFractureCell();

  // -- operators (energy)
  auto cvs_matrix = Teuchos::rcp(new CompositeVectorSpace());
  cvs_matrix->SetMesh(mesh_matrix)->SetGhosted(true)->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  auto cvs_fracture = Teuchos::rcp(new CompositeVectorSpace());
  cvs_fracture->SetMesh(mesh_fracture)->SetGhosted(true)->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  adv_coupling_matrix_ = AddCouplingFluxes_(fia.get_cvs_matrix(),
                                            fia.get_cvs_fracture(),
                                            fia.get_inds_matrix(),
                                            fia.get_inds_fracture(),
                                            fia.get_values(),
                                            1,
                                            op_tree_matrix_);

  adv_coupling_pc_ = AddCouplingFluxes_(fia.get_cvs_matrix(),
                                        fia.get_cvs_fracture(),
                                        fia.get_inds_matrix(),
                                        fia.get_inds_fracture(),
                                        fia.get_values(),
                                        1,
                                        op_tree_pc_);

  // create global matrix
  op_tree_pc_->SymbolicAssembleMatrix();

  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "solution vector:\n";
    solution_->Print(*vo_->os(), false);
    *vo_->os() << "\nmatrix:" << std::endl
               << op_tree_matrix_->PrintDiagnostics() << "preconditioner:" << std::endl
               << op_tree_pc_->PrintDiagnostics() << std::endl
               << vo_->color("green") << "Initialization of PK is complete: my dT=" << get_dt()
               << vo_->reset() << std::endl
               << std::endl;
  }
}


/* *******************************************************************
* Performs one time step.
******************************************************************* */
bool
FlowEnergyMatrixFracture_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  auto counter = Teuchos::TimeMonitor::getNewCounter("Flow-Energy MPC PK");
  Teuchos::TimeMonitor tm(*counter);

  // make copy of evaluators
  std::vector<Key> names = { "saturation_liquid", "water_storage", "energy" };
  Teuchos::RCP<CompositeVector> copies[6];

  int k(0), nnames(names.size());
  for (int i = 0; i < nnames; ++i) {
    SwapEvaluatorField_(names[i], copies[k], copies[k + 1]);
    k += 2;
  }

  // make copy of primary unknowns
  // save a copy of solution, i.e. primary variables
  TreeVector solution_copy(*solution_);

  bool fail;
  try {
    fail = PK_MPCStrong<PK_BDF>::AdvanceStep(t_old, t_new, reinit);
  } catch (Errors::CutTimeStep& e) {
    if (vo_->os_OK(Teuchos::VERB_HIGH)) {
      Teuchos::OSTab tab = vo_->getOSTab();
      *vo_->os() << e.what() << std::endl;
    }
    fail = false;
  }

  if (fail) {
    std::string passwd("");
    k = 0;
    for (int i = 0; i < nnames; ++i) {
      if (S_->HasRecord(names[i])) {
        S_->GetW<CV_t>("prev_" + names[i], passwd) = *(copies[k]);
        S_->GetW<CV_t>("fracture-prev_" + names[i], passwd) = *(copies[k + 1]);
      }
      k += 2;
    }

    // recover the original solution
    *solution_ = solution_copy;
    ChangedSolution();

    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "Step failed. Restored [fracture-]{ " << names[0] << ", " << names[1] << ", "
               << names[2] << ", pressure, temperature }" << std::endl;
  }

  return fail;
}


/* *******************************************************************
* Residual evaluation
******************************************************************* */
void
FlowEnergyMatrixFracture_PK::FunctionalResidual(double t_old,
                                                double t_new,
                                                Teuchos::RCP<TreeVector> u_old,
                                                Teuchos::RCP<TreeVector> u_new,
                                                Teuchos::RCP<TreeVector> f)
{
  // generate local matrices and apply sources and boundary conditions
  PK_MPCStrong<PK_BDF>::FunctionalResidual(t_old, t_new, u_old, u_new, f);

  // add contribution of coupling terms to the residual
  TreeVector g(*f);

  UpdateCouplingFluxes_(adv_coupling_matrix_);

  op_tree_matrix_->Apply(*u_new, g);
  f->Update(1.0, g, 1.0);

  // convergence control
  f->NormInf(&residual_norm_);
}


/* *******************************************************************
* Preconditioner update
******************************************************************* */
void
FlowEnergyMatrixFracture_PK::UpdatePreconditioner(double t,
                                                  Teuchos::RCP<const TreeVector> up,
                                                  double dt)
{
  // generate local matrices and apply boundary conditions
  PK_MPCStrong<PK_BDF>::UpdatePreconditioner(t, up, dt);

  UpdateCouplingFluxes_(adv_coupling_pc_);

  std::string pc_name = ti_list_->get<std::string>("preconditioner");
  Teuchos::ParameterList pc_list = preconditioner_list_->sublist(pc_name);

  op_tree_pc_->AssembleMatrix();

  // block indices for preconditioner are (0, 1, 0, 1)
  auto smap = op_tree_pc_->get_row_supermap();
  // NOTE: this is freed by Hypre
  auto block_indices = Teuchos::rcp(new std::vector<int>(smap->Map()->NumMyElements()));

  std::vector<std::string> comps = { "face", "cell" };
  for (int n = 0; n < 4; ++n) {
    int id = (n == 0 || n == 2) ? 0 : 1;

    for (int k = 0; k < 2; ++k) {
      const auto& inds = smap->Indices(n, comps[k], 0);
      for (int i = 0; i != inds.size(); ++i) (*block_indices)[inds[i]] = id;
    }
  }
  auto block_ids = std::make_pair(2, block_indices);

  op_tree_pc_->set_coloring(2, block_indices);
  op_tree_pc_->set_inverse_parameters(pc_list);

  // create a stronger preconditioner by wrapping one inside an iterative solver
  if (ti_list_->isParameter("preconditioner enhancement")) {
    std::string solver_name = ti_list_->get<std::string>("preconditioner enhancement");
    AMANZI_ASSERT(linear_operator_list_->isSublist(solver_name));
    auto tmp_plist = linear_operator_list_->sublist(solver_name);
    op_pc_solver_ = AmanziSolvers::createIterativeMethod(tmp_plist, op_tree_pc_);
  } else {
    op_pc_solver_ = op_tree_pc_;
  }

  op_pc_solver_->InitializeInverse();
  op_pc_solver_->ComputeInverse();
}


/* *******************************************************************
* Application of preconditioner
******************************************************************* */
int
FlowEnergyMatrixFracture_PK::ApplyPreconditioner(Teuchos::RCP<const TreeVector> X,
                                                 Teuchos::RCP<TreeVector> Y)
{
  Y->PutScalar(0.0);
  return op_pc_solver_->ApplyInverse(*X, *Y);
}


/* *******************************************************************
* Coupling term: diffusion fluxes between matrix and fracture
******************************************************************* */
std::vector<Teuchos::RCP<Operators::PDE_CouplingFlux>>
FlowEnergyMatrixFracture_PK::AddCouplingFluxes_(
  const Teuchos::RCP<CompositeVectorSpace>& cvs_matrix,
  const Teuchos::RCP<CompositeVectorSpace>& cvs_fracture,
  std::shared_ptr<const std::vector<std::vector<int>>> inds_matrix,
  std::shared_ptr<const std::vector<std::vector<int>>> inds_fracture,
  std::shared_ptr<const std::vector<double>> values,
  int i,
  Teuchos::RCP<Operators::TreeOperator>& op_tree)
{
  Teuchos::ParameterList oplist;

  auto op00 = op_tree->get_block(0, 0)->get_operator_block(i, i);
  auto op01 = op_tree->get_block(0, 1)->get_operator_block(i, i);
  auto op10 = op_tree->get_block(1, 0)->get_operator_block(i, i);
  auto op11 = op_tree->get_block(1, 1)->get_operator_block(i, i);

  Teuchos::RCP<Operators::PDE_CouplingFlux> op_coupling00, op_coupling11;
  // add diagonal
  op_coupling00 = Teuchos::rcp(new Operators::PDE_CouplingFlux(
    oplist, cvs_matrix, cvs_matrix, inds_matrix, inds_matrix, op00));
  op_coupling00->Setup(values, 1.0);
  op_coupling00->UpdateMatrices(Teuchos::null, Teuchos::null);

  op_coupling11 = Teuchos::rcp(new Operators::PDE_CouplingFlux(
    oplist, cvs_fracture, cvs_fracture, inds_fracture, inds_fracture, op11));
  op_coupling11->Setup(values, 1.0);
  op_coupling11->UpdateMatrices(Teuchos::null, Teuchos::null);

  if (op00 == Teuchos::null)
    op_tree->get_block(0, 0)->set_operator_block(i, i, op_coupling00->global_operator());

  if (op11 == Teuchos::null)
    op_tree->get_block(1, 1)->set_operator_block(i, i, op_coupling11->global_operator());

  auto op_coupling01 = Teuchos::rcp(new Operators::PDE_CouplingFlux(
    oplist, cvs_matrix, cvs_fracture, inds_matrix, inds_fracture, op01));
  op_coupling01->Setup(values, -1.0);
  op_coupling01->UpdateMatrices(Teuchos::null, Teuchos::null);

  auto op_coupling10 = Teuchos::rcp(new Operators::PDE_CouplingFlux(
    oplist, cvs_fracture, cvs_matrix, inds_fracture, inds_matrix, op10));
  op_coupling10->Setup(values, -1.0);
  op_coupling10->UpdateMatrices(Teuchos::null, Teuchos::null);

  if (op10 == Teuchos::null)
    op_tree->get_block(1, 0)->set_operator_block(i, i, op_coupling10->global_operator());

  if (op01 == Teuchos::null)
    op_tree->get_block(0, 1)->set_operator_block(i, i, op_coupling01->global_operator());

  // return operators
  std::vector<Teuchos::RCP<Operators::PDE_CouplingFlux>> ops;
  ops.push_back(op_coupling00);
  ops.push_back(op_coupling01);
  ops.push_back(op_coupling10);
  ops.push_back(op_coupling11);

  return ops;
}


/* *******************************************************************
* Compute coupling fluxes
******************************************************************* */
void
FlowEnergyMatrixFracture_PK::UpdateCouplingFluxes_(
  const std::vector<Teuchos::RCP<Operators::PDE_CouplingFlux>>& adv_coupling)
{
  S_->Get<CV_t>("enthalpy").ScatterMasterToGhosted("cell");
  S_->Get<CV_t>("molar_density_liquid").ScatterMasterToGhosted("cell");
  S_->Get<CV_t>("temperature").ScatterMasterToGhosted("cell");
  S_->Get<CV_t>(matrix_vol_flowrate_key_).ScatterMasterToGhosted("face");

  // extract enthalpy fields
  S_->GetEvaluator("enthalpy").Update(*S_, "enthalpy");
  const auto& H_m = *S_->Get<CV_t>("enthalpy").ViewComponent("cell", true);
  const auto& T_m = *S_->Get<CV_t>("temperature").ViewComponent("cell", true);
  const auto& n_l_m = *S_->Get<CV_t>("molar_density_liquid").ViewComponent("cell", true);

  S_->GetEvaluator("fracture-enthalpy").Update(*S_, "fracture-enthalpy");
  const auto& H_f = *S_->Get<CV_t>("fracture-enthalpy").ViewComponent("cell", true);
  const auto& T_f = *S_->Get<CV_t>("fracture-temperature").ViewComponent("cell", true);
  const auto& n_l_f = *S_->Get<CV_t>("fracture-molar_density_liquid").ViewComponent("cell", true);

  // update coupling terms for advection
  int ncells_owned_f =
    mesh_fracture_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);
  auto values1 = std::make_shared<std::vector<double>>(2 * ncells_owned_f, 0.0);
  auto values2 = std::make_shared<std::vector<double>>(2 * ncells_owned_f, 0.0);

  int np(0), dir, shift;
  AmanziMesh::Entity_ID_List cells;
  const auto& flux = *S_->Get<CV_t>(matrix_vol_flowrate_key_).ViewComponent("face", true);
  const auto& mmap = flux.Map();
  for (int c = 0; c < ncells_owned_f; ++c) {
    int f = mesh_fracture_->getEntityParent(AmanziMesh::Entity_kind::CELL, c);
    int first = mmap.FirstPointInElement(f);

    cells = mesh_domain_->getFaceCells(f, AmanziMesh::Parallel_type::ALL);
    int ncells = cells.size();
    mesh_domain_->getFaceNormal(f, cells[0], &dir);
    shift = Operators::UniqueIndexFaceToCells(*mesh_domain_, f, cells[0]);

    for (int k = 0; k < ncells; ++k) {
      double tmp = flux[0][first + shift] * dir;

      // since we multiply by temperature, the model for the flux is
      // q (\eta H / T) * T for both matrix and preconditioner
      if (tmp > 0) {
        int c1 = cells[k];
        double factor = H_m[0][c1] * n_l_m[0][c1] / T_m[0][c1];
        (*values1)[np] = tmp * factor;
      } else {
        double factor = H_f[0][c] * n_l_f[0][c] / T_f[0][c];
        (*values2)[np] = -tmp * factor;
      }

      dir = -dir;
      shift = 1 - shift;
      np++;
    }
  }

  // setup coupling operators with new data
  adv_coupling[0]->Setup(values1, 1.0);
  adv_coupling[0]->UpdateMatrices(Teuchos::null, Teuchos::null);

  adv_coupling[1]->Setup(values2, -1.0);
  adv_coupling[1]->UpdateMatrices(Teuchos::null, Teuchos::null);

  adv_coupling[2]->Setup(values1, -1.0);
  adv_coupling[2]->UpdateMatrices(Teuchos::null, Teuchos::null);

  adv_coupling[3]->Setup(values2, 1.0);
  adv_coupling[3]->UpdateMatrices(Teuchos::null, Teuchos::null);
}


/* *******************************************************************
* Copy: Evaluator (BASE) -> Field (prev_BASE)
******************************************************************* */
void
FlowEnergyMatrixFracture_PK::SwapEvaluatorField_(const Key& key,
                                                 Teuchos::RCP<CompositeVector>& fdm_copy,
                                                 Teuchos::RCP<CompositeVector>& fdf_copy)
{
  if (!S_->HasRecord(key)) return;

  // matrix
  Key ev_key, fd_key;
  ev_key = key;
  fd_key = "prev_" + key;

  std::string passwd("");
  S_->GetEvaluator(ev_key).Update(*S_, passwd);
  {
    const auto& ev = S_->Get<CV_t>(ev_key);
    auto& fd = S_->GetW<CV_t>(fd_key, passwd);
    fdm_copy = Teuchos::rcp(new CompositeVector(fd));
    fd = ev;
  }

  // fracture
  ev_key = "fracture-" + key;
  fd_key = "fracture-prev_" + key;
  if (!S_->HasRecord(ev_key)) return;

  S_->GetEvaluator(ev_key).Update(*S_, passwd);
  {
    const auto& ev = S_->Get<CV_t>(ev_key);
    auto& fd = S_->GetW<CV_t>(fd_key, passwd);

    fdf_copy = Teuchos::rcp(new CompositeVector(fd));
    fd = ev;
  }
}


/* ******************************************************************
* Check solution and fields for convergence
****************************************************************** */
double
FlowEnergyMatrixFracture_PK::ErrorNorm(Teuchos::RCP<const TreeVector> u,
                                       Teuchos::RCP<const TreeVector> du)
{
  double error = PK_MPCStrong<PK_BDF>::ErrorNorm(u, du);

  // residual control for energy (note that we cannot do it at
  // the lower level due to coupling terms.
  const auto mesh_f = S_->GetMesh("fracture");
  const auto& mol_fc = *S_->Get<CV_t>("fracture-molar_density_liquid").ViewComponent("cell");

  int ncells = mol_fc.MyLength();
  double mean_energy, error_r(0.0), mass(0.0);

  for (int c = 0; c < ncells; ++c) {
    mass += mol_fc[0][c] * mesh_f->getCellVolume(c); // reference cell energy
  }
  if (ncells > 0) {
    mean_energy = 76.0 * mass / ncells;
    error_r = (residual_norm_ * dt_) / mean_energy;
  }

  error = std::max(error, error_r);

#ifdef HAVE_MPI
  double tmp = error;
  u->Comm()->MaxAll(&tmp, &error, 1);
#endif

  return error;
}


/* *******************************************************************
* This should be refactored, see for simialr function in PK_Physical
******************************************************************* */
void
FlowEnergyMatrixFracture_PK::AddDefaultPrimaryEvaluator_(const Key& key, const Tag& tag)
{
  Teuchos::ParameterList elist(key);
  elist.set<std::string>("tag", tag.get()).set<std::string>("evaluator name", key);
  auto eval = Teuchos::rcp(new EvaluatorPrimary<CompositeVector, CompositeVectorSpace>(elist));
  S_->SetEvaluator(key, tag, eval);
}

} // namespace Amanzi
