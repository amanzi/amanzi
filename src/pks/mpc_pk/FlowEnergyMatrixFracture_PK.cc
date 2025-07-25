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

#include "FlowMatrixFracture_PK.hh"
#include "EnergyMatrixFracture_PK.hh"
#include "InverseFactory.hh"
#include "PDE_CouplingFlux.hh"
#include "PDE_DiffusionFracturedMatrix.hh"
#include "TreeOperator.hh"
#include "StateArchive.hh"
#include "SuperMap.hh"
#include "UniqueLocalIndex.hh"

#include "FlowEnergyMatrixFracture_PK.hh"
#include "FractureInsertion.hh"
#include "FractureInsertion_Helper.hh"
#include "HeatDiffusionMatrixFracture.hh"
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

  Teuchos::ParameterList vlist;
  vlist.sublist("verbose object") = plist_->sublist("verbose object");
  vo_ = Teuchos::rcp(new VerboseObject("CoupledThermalFlow_PK", vlist));
}


/* *******************************************************************
* Physics-based setup of PK.
******************************************************************* */
void
FlowEnergyMatrixFracture_PK::Setup()
{
  mesh_domain_ = S_->GetMesh();
  mesh_fracture_ = S_->GetMesh("fracture");

  // primary and secondary fields for matrix affected by non-uniform
  // distribution of DOFs, so we need to define them here
  // -- flow: pressure and flux
  auto cvs = Operators::CreateFracturedMatrixCVS(mesh_domain_, mesh_fracture_);
  if (!S_->HasRecord("pressure")) {
    *S_->Require<CV_t, CVS_t>("pressure", Tags::DEFAULT).SetMesh(mesh_domain_)->SetGhosted(true) =
      *cvs;
    AddDefaultPrimaryEvaluator(S_, "pressure", Tags::DEFAULT);
  }

  // -- energy: temperature
  if (!S_->HasRecord("temperature")) {
    *S_->Require<CV_t, CVS_t>("temperature", Tags::DEFAULT)
       .SetMesh(mesh_domain_)
       ->SetGhosted(true) = *cvs;
    AddDefaultPrimaryEvaluator(S_, "temperature", Tags::DEFAULT);
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

  auto& fflow = glist_->sublist("PKs").sublist(pks1[1]).sublist("physical models and assumptions");
  fflow.set<std::string>("coupled matrix fracture flow", "fracture");

  // -- energy
  auto& menergy =
    glist_->sublist("PKs").sublist(pks2[0]).sublist("physical models and assumptions");
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

  auto mpc_flow = Teuchos::rcp_dynamic_cast<FlowMatrixFracture_PK>(sub_pks_[0]);
  auto mpc_energy = Teuchos::rcp_dynamic_cast<EnergyMatrixFracture_PK>(sub_pks_[1]);

  // NOTE: In this PK, the solution is 2-deep, with 2 sub-PKs (flow, energy)
  // each with 2 sub-pks of their own (matrix, fracture).  Currently
  // TreeOperator cannot handle this, so instead we must flatten the map.
  auto tvs = Teuchos::rcp(new TreeVectorSpace(solution_->Map()));
  op_tree_pc_ = Teuchos::rcp(new Operators::TreeOperator(tvs));

  op_tree_pc_->set_block(0, 0, mpc_flow->op_tree_pc()->Clone());
  op_tree_pc_->set_block(1, 1, mpc_energy->op_tree_pc()->Clone());
  /*
  for (int i = 0; i < 2; ++i) {
    const auto& row = solution_->SubVector(i)->get_map();
    for (int j = 0; j < 2; ++j) {
      const auto& col = solution_->SubVector(j)->get_map();
      if (i != j) op_tree_pc_->set_block(i, j, Teuchos::rcp(new Operators::TreeOperator(row, col)));
    }
  }
  */

  // create global matrix
  op_tree_pc_->SymbolicAssembleMatrix();

  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "solution vector:\n";
    solution_->Print(*vo_->os(), false);
    *vo_->os() << "\npreconditioner:" << std::endl
               << op_tree_pc_->PrintDiagnostics() << std::endl
               << vo_->color("green") << "Initialization of PK is complete: my dT=" << get_dt()
               << vo_->reset() << std::endl
               << std::endl;
  }
}


/* *******************************************************************
* Performs one timestep.
******************************************************************* */
bool
FlowEnergyMatrixFracture_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  auto counter = Teuchos::TimeMonitor::getNewCounter("Flow-Energy MPC PK");
  Teuchos::TimeMonitor tm(*counter);

  // make copy of conservative fields
  std::vector<Key> names = {
    "pressure", "temperature", "fracture-pressure", "fracture-temperature"
  };
  std::vector<Key> fields = {
    "saturation_liquid",          "water_storage",          "energy",
    "fracture-saturation_liquid", "fracture-water_storage", "fracture-energy"
  };
  StateArchive archive(S_, vo_);
  archive.Add(names, Tags::DEFAULT);
  archive.CopyFieldsToPrevFields(fields, "", true);

  bool fail;
  try {
    fail = PK_MPCStrong<PK_BDF>::AdvanceStep(t_old, t_new, reinit);
  } catch (Errors::CutTimestep& e) {
    if (vo_->os_OK(Teuchos::VERB_HIGH)) {
      Teuchos::OSTab tab = vo_->getOSTab();
      *vo_->os() << e.what() << std::endl;
    }
    fail = false;
  }
  if (fail) archive.Restore("");

  return fail;
}


/* *******************************************************************
* Residual evaluation
******************************************************************* */
void
FlowEnergyMatrixFracture_PK::FunctionalResidual(double t_old,
                                                double t_new,
                                                Teuchos::RCP<const TreeVector> u_old,
                                                Teuchos::RCP<TreeVector> u_new,
                                                Teuchos::RCP<TreeVector> f)
{
  // flow
  auto u_old0 = u_old->SubVector(0);
  auto u_new0 = u_new->SubVector(0);
  auto f0 = f->SubVector(0);
  sub_pks_[0]->FunctionalResidual(t_old, t_new, u_old0, u_new0, f0);

  // update molar flux
  std::vector<Key> keys({ "molar_flow_rate", "fracture-molar_flow_rate" });
  for (int i = 0; i < 2; ++i) {
    auto mol_flowrate = S_->GetPtrW<CV_t>(keys[i], Tags::DEFAULT, "");
    auto mpc = Teuchos::rcp_dynamic_cast<PK_MPC<PK_BDF>>(sub_pks_[0]);
    auto op0 = mpc->getSubPK(i)->my_pde(Operators::PDEType::PDE_DIFFUSION);
    op0->UpdateFlux(u_new0->SubVector(i)->Data().ptr(), mol_flowrate.ptr());

    if (Keys::getVarName(mpc->getSubPK(i)->name()) == "darcy") {
      double molar_mass = S_->Get<double>("const_fluid_molar_mass");
      mol_flowrate->Scale(1.0 / molar_mass);
    }
  }

  // energy
  auto u_old1 = u_old->SubVector(1);
  auto u_new1 = u_new->SubVector(1);
  auto f1 = f->SubVector(1);
  sub_pks_[1]->FunctionalResidual(t_old, t_new, u_old1, u_new1, f1);

  // convergence control
  f->SubVector(1)->SubVector(1)->NormInf(&residual_norm_);
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
      for (int i = 0; i != inds.size() ; ++i) (*block_indices)[inds[i]] = id;
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
  const auto& energy_fc = *S_->Get<CV_t>("fracture-energy").ViewComponent("cell");

  int ncells = energy_fc.MyLength();
  double mean_energy(0.0), error_r(0.0), energy(0.0);

  for (int c = 0; c < ncells; ++c) {
    energy += energy_fc[0][c] * mesh_f->getCellVolume(c); // reference cell energy
  }
  if (ncells > 0) {
    mean_energy = energy / ncells;
    error_r = (residual_norm_ * dt_) / mean_energy;
  }

  // removed since coupling fluxes are ignored
  // error = std::max(error, error_r);

#ifdef HAVE_MPI
  double tmp = error;
  u->Comm()->MaxAll(&tmp, &error, 1);
#endif

  return error;
}

} // namespace Amanzi
