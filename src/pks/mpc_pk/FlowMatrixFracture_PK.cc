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

  Process kernel that couples flow in matrix and fracture network.
*/

#include "CommonDefs.hh"
#include "EvaluatorPrimary.hh"
#include "InverseFactory.hh"
#include "PDE_CouplingFlux.hh"
#include "PDE_DiffusionFracturedMatrix.hh"
#include "TreeOperator.hh"

#include "FractureInsertion.hh"
#include "FlowMatrixFracture_PK.hh"
#include "PK_MPCStrong.hh"
#include "PK_Utils.hh"

namespace Amanzi {

using CV_t = CompositeVector;
using CVS_t = CompositeVectorSpace;

/* *******************************************************************
* Constructor
******************************************************************* */
FlowMatrixFracture_PK::FlowMatrixFracture_PK(Teuchos::ParameterList& pk_tree,
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
  vo_ = Teuchos::rcp(new VerboseObject("CoupledFlow_PK", vlist));
}


/* *******************************************************************
* Physics-based setup of PK.
******************************************************************* */
void
FlowMatrixFracture_PK::Setup()
{
  mesh_matrix_ = S_->GetMesh();
  mesh_fracture_ = S_->GetMesh("fracture");

  // primary and secondary fields for matrix affected by non-uniform
  // distribution of DOFs
  // -- pressure
  auto cvs = Operators::CreateFracturedMatrixCVS(mesh_matrix_, mesh_fracture_);
  if (!S_->HasRecord("pressure")) {
    *S_->Require<CV_t, CVS_t>("pressure", Tags::DEFAULT).SetMesh(mesh_matrix_)->SetGhosted(true) =
      *cvs;
    AddDefaultPrimaryEvaluator(S_, "pressure", Tags::DEFAULT);
  }

  // -- darcy flux
  if (!S_->HasRecord("volumetric_flow_rate")) {
    std::string name("face");
    auto mmap = cvs->Map("face", false);
    auto gmap = cvs->Map("face", true);
    S_->Require<CV_t, CVS_t>("volumetric_flow_rate", Tags::DEFAULT)
      .SetMesh(mesh_matrix_)
      ->SetGhosted(true)
      ->SetComponent(name, AmanziMesh::FACE, mmap, gmap, 1);
    AddDefaultPrimaryEvaluator(S_, "volumetric_flow_rate", Tags::DEFAULT);
  }

  // -- darcy flux for fracture
  if (!S_->HasRecord("fracture-volumetric_flow_rate")) {
    auto cvs2 = Operators::CreateManifoldCVS(mesh_fracture_);
    *S_->Require<CV_t, CVS_t>("fracture-volumetric_flow_rate", Tags::DEFAULT)
       .SetMesh(mesh_fracture_)
       ->SetGhosted(true) = *cvs2;
    AddDefaultPrimaryEvaluator(S_, "fracture-volumetric_flow_rate", Tags::DEFAULT);
  }

  // additional fields and evaluators related to coupling
  Key normal_permeability_key_("fracture-normal_permeability");
  if (!S_->HasRecord(normal_permeability_key_)) {
    S_->Require<CV_t, CVS_t>(normal_permeability_key_, Tags::DEFAULT)
      .SetMesh(mesh_fracture_)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  // inform dependent PKs about coupling
  // -- flow (matrix)
  std::vector<std::string> pks = plist_->get<Teuchos::Array<std::string>>("PKs order").toVector();
  Teuchos::ParameterList& mflow =
    glist_->sublist("PKs").sublist(pks[0]).sublist("physical models and assumptions");
  mflow.set<std::string>("coupled matrix fracture flow", "matrix")
    .set<bool>("use bulk modulus", true);

  // -- flow (fracture)
  Teuchos::ParameterList& fflow =
    glist_->sublist("PKs").sublist(pks[1]).sublist("physical models and assumptions");
  fflow.set<std::string>("coupled matrix fracture flow", "fracture");

  // modify time integrator (only for Darcy)
  // ti_list_->sublist("BDF1").set<bool>("freeze preconditioner", true);

  // process other PKs.
  PK_MPCStrong<PK_BDF>::Setup();
}


/* *******************************************************************
* Initialization create a tree operator to assemble global matrix
******************************************************************* */
void
FlowMatrixFracture_PK::Initialize()
{
  PK_MPCStrong<PK_BDF>::Initialize();

  // since solution's map could be anything, to create a global operator,
  // we have to rely on pk's operator structure.
  auto tvs = Teuchos::rcp(new TreeVectorSpace());

  for (const auto& pk : sub_pks_) {
    auto& cvs = pk->my_operator(Operators::OPERATOR_MATRIX)->get_domain_map();
    auto tmp = Teuchos::rcp(new TreeVectorSpace(cvs));
    tvs->PushBack(tmp);
  }

  // we assume that 0 and 1 correspond to matrix and fracture, respectively
  // to avoid modifying original operators, we clone them.
  auto op0 = sub_pks_[0]->my_operator(Operators::OPERATOR_MATRIX)->Clone();
  auto op1 = sub_pks_[1]->my_operator(Operators::OPERATOR_MATRIX)->Clone();

  // off-diagonal blocks are coupled PDEs
  // -- minimum composite vector spaces containing the coupling term
  auto& mmap = solution_->SubVector(0)->Data()->ViewComponent("face", false)->Map();
  auto& gmap = solution_->SubVector(0)->Data()->ViewComponent("face", true)->Map();

  // -- indices transmissibimility coefficients for matrix-fracture flux
  const auto& kn = *S_->Get<CV_t>("fracture-normal_permeability").ViewComponent("cell");
  double gravity = norm(S_->Get<AmanziGeometry::Point>("gravity"));

  FractureInsertion fi(mesh_matrix_, mesh_fracture_);
  fi.InitMatrixFaceToFractureCell(Teuchos::rcpFromRef(mmap), Teuchos::rcpFromRef(gmap));
  double scale = (sub_pks_[0]->name() == "darcy") ? 1.0 : 1.0 / CommonDefs::MOLAR_MASS_H2O;
  fi.SetValues(kn, scale / gravity);

  // -- operators
  Teuchos::ParameterList oplist;

  auto op_coupling00 = Teuchos::rcp(new Operators::PDE_CouplingFlux(oplist,
                                                                    fi.get_cvs_matrix(),
                                                                    fi.get_cvs_matrix(),
                                                                    fi.get_inds_matrix(),
                                                                    fi.get_inds_matrix(),
                                                                    Teuchos::null));
  op0->OpPushBack(op_coupling00->local_op());
  op_coupling00->Setup(fi.get_values(), 1.0);
  op_coupling00->UpdateMatrices(Teuchos::null, Teuchos::null);

  auto op_coupling01 = Teuchos::rcp(new Operators::PDE_CouplingFlux(oplist,
                                                                    fi.get_cvs_matrix(),
                                                                    fi.get_cvs_fracture(),
                                                                    fi.get_inds_matrix(),
                                                                    fi.get_inds_fracture()));
  op_coupling01->Setup(fi.get_values(), -1.0);
  op_coupling01->UpdateMatrices(Teuchos::null, Teuchos::null);

  auto op_coupling10 = Teuchos::rcp(new Operators::PDE_CouplingFlux(oplist,
                                                                    fi.get_cvs_fracture(),
                                                                    fi.get_cvs_matrix(),
                                                                    fi.get_inds_fracture(),
                                                                    fi.get_inds_matrix()));
  op_coupling10->Setup(fi.get_values(), -1.0);
  op_coupling10->UpdateMatrices(Teuchos::null, Teuchos::null);

  auto op_coupling11 = Teuchos::rcp(new Operators::PDE_CouplingFlux(oplist,
                                                                    fi.get_cvs_fracture(),
                                                                    fi.get_cvs_fracture(),
                                                                    fi.get_inds_fracture(),
                                                                    fi.get_inds_fracture(),
                                                                    Teuchos::null));
  op1->OpPushBack(op_coupling11->local_op());
  op_coupling11->Setup(fi.get_values(), 1.0);
  op_coupling11->UpdateMatrices(Teuchos::null, Teuchos::null);

  // create global matrix
  // -- tree matrix (for other MPCs)
  op_tree_matrix_ = Teuchos::rcp(new Operators::TreeOperator(tvs));
  op_tree_matrix_->set_operator_block(0, 0, op0);
  op_tree_matrix_->set_operator_block(1, 1, op1);
  op_tree_matrix_->set_operator_block(0, 1, op_coupling01->global_operator());
  op_tree_matrix_->set_operator_block(1, 0, op_coupling10->global_operator());

  // -- tree preconditioner
  auto pc0 = sub_pks_[0]->my_operator(Operators::OPERATOR_PRECONDITIONER_RAW)->Clone();
  auto pc1 = sub_pks_[1]->my_operator(Operators::OPERATOR_PRECONDITIONER_RAW)->Clone();

  pc0->OpPushBack(op_coupling00->local_op());
  pc1->OpPushBack(op_coupling11->local_op());

  op_tree_pc_ = Teuchos::rcp(new Operators::TreeOperator(tvs));
  op_tree_pc_->set_operator_block(0, 0, pc0);
  op_tree_pc_->set_operator_block(1, 1, pc1);
  op_tree_pc_->set_operator_block(0, 1, op_coupling01->global_operator());
  op_tree_pc_->set_operator_block(1, 0, op_coupling10->global_operator());

  // -- configure preconditioner
  sub_pks_[0]->my_pde(Operators::PDE_DIFFUSION)->ApplyBCs(true, true, true);

  std::string name = ti_list_->get<std::string>("preconditioner");
  std::string ls_name = ti_list_->get<std::string>("preconditioner enhancement", "none");
  auto inv_list = AmanziSolvers::mergePreconditionerSolverLists(
    name, *preconditioner_list_, ls_name, *linear_operator_list_, true);
  op_tree_pc_->set_inverse_parameters(inv_list);
  op_tree_pc_->InitializeInverse();

  // -- tree coupling matrix
  op_matrix_ = Teuchos::rcp(new Operators::TreeOperator(tvs));

  op_matrix_->set_operator_block(0, 0, op_coupling00->global_operator());
  op_matrix_->set_operator_block(1, 1, op_coupling11->global_operator());
  op_matrix_->set_operator_block(0, 1, op_coupling01->global_operator());
  op_matrix_->set_operator_block(1, 0, op_coupling10->global_operator());

  // stationary solve is modelled with large dt. To pick the correct
  // boundary conditions, dt is negative. This assumes that we are at
  // the beginning of simulation.
  if (ti_list_->isSublist("initialization")) {
    // bool wells_on = ti_list_->sublist("initialization").get<bool>("active wells", false);
    double dt(-1e+98), dt_solver;
    bool fail = time_stepper_->TimeStep(dt, dt_solver, solution_);
    if (fail) Exceptions::amanzi_throw("Solver for coupled flow did not converge.");
  }

  if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "coupling matrix:" << std::endl << op_matrix_->PrintDiagnostics() << std::endl;
    *vo_->os() << "preconditioner:" << std::endl
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
FlowMatrixFracture_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  // create copies of conservative fields
  std::vector<std::string> fields = { "saturation_liquid", "fracture-saturation_liquid" };
  if (sub_pks_[0]->name() == "richards") {
    fields.push_back("water_storage");
    fields.push_back("fracture-water_storage");
  }

  StateArchive archive(S_, vo_);
  archive.Add(fields, {}, {}, Tags::DEFAULT, name());
  archive.CopyFieldsToPrevFields("");

  bool fail = PK_MPCStrong<PK_BDF>::AdvanceStep(t_old, t_new, reinit);

  if (fail) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "Step failed." << std::endl;

    archive.Restore("");
  }

  // update some fields, we cannot move this to commit step due to "initialize"
  S_->GetW<CV_t>("fracture-prev_aperture", Tags::DEFAULT, "") =
    S_->Get<CV_t>("fracture-aperture", Tags::DEFAULT);

  return fail;
}


/* *******************************************************************
* Residual evaluation
******************************************************************* */
void
FlowMatrixFracture_PK::FunctionalResidual(double t_old,
                                          double t_new,
                                          Teuchos::RCP<TreeVector> u_old,
                                          Teuchos::RCP<TreeVector> u_new,
                                          Teuchos::RCP<TreeVector> f)
{
  PK_MPCStrong<PK_BDF>::FunctionalResidual(t_old, t_new, u_old, u_new, f);

  int ierr = op_matrix_->Apply(*u_new, *f, 1.0);
  AMANZI_ASSERT(!ierr);
}


/* *******************************************************************
* Preconditioner update
******************************************************************* */
void
FlowMatrixFracture_PK::UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> u, double dt)
{
  PK_MPCStrong<PK_BDF>::UpdatePreconditioner(t, u, dt);
  op_tree_pc_->ComputeInverse();
}


/* *******************************************************************
* Application of preconditioner
******************************************************************* */
int
FlowMatrixFracture_PK::ApplyPreconditioner(Teuchos::RCP<const TreeVector> X,
                                           Teuchos::RCP<TreeVector> Y)
{
  Y->PutScalar(0.0);
  int ok = op_tree_pc_->ApplyInverse(*X, *Y);
  return ok;
}

} // namespace Amanzi
