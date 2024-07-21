/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
           Daniil Svyatskiy
*/

/*
  MPC

  Process kernel that couples heat conduction in matrix and fracture network.
*/

#include "CommonDefs.hh"
#include "EvaluatorPrimary.hh"
#include "InverseFactory.hh"
#include "PDE_CouplingFlux.hh"
#include "PDE_DiffusionFracturedMatrix.hh"
#include "StateArchive.hh"
#include "TreeOperator.hh"

#include "EnergyMatrixFracture_PK.hh"
#include "FractureInsertion.hh"
#include "FractureInsertion_Helper.hh"
#include "HeatDiffusionMatrixFracture.hh"
#include "PK_MPCStrong.hh"

namespace Amanzi {

using CV_t = CompositeVector;
using CVS_t = CompositeVectorSpace;

/* *******************************************************************
* Constructor
******************************************************************* */
EnergyMatrixFracture_PK::EnergyMatrixFracture_PK(Teuchos::ParameterList& pk_tree,
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
    Errors::Message msg(messagestream.str());
    Exceptions::amanzi_throw(msg);
  }

  preconditioner_list_ = Teuchos::sublist(glist, "preconditioners", true);
  linear_operator_list_ = Teuchos::sublist(glist, "solvers", true);
  ti_list_ = Teuchos::sublist(plist_, "time integrator", true);

  Teuchos::ParameterList vlist;
  vlist.sublist("verbose object") = plist_->sublist("verbose object");
  vo_ = Teuchos::rcp(new VerboseObject("CoupledEnergy_PK", vlist));
}


/* *******************************************************************
* Physics-based setup of PK.
******************************************************************* */
void
EnergyMatrixFracture_PK::Setup()
{
  mesh_matrix_ = S_->GetMesh();
  mesh_fracture_ = S_->GetMesh("fracture");

  heat_diffusion_to_matrix_key_ = "fracture-heat_diffusion_to_matrix";

  // primary and secondary fields for matrix affected by non-uniform
  // distribution of DOFs
  // -- temperature
  auto cvs = Operators::CreateFracturedMatrixCVS(mesh_matrix_, mesh_fracture_);
  if (!S_->HasRecord("temperature")) {
    *S_->Require<CV_t, CVS_t>("temperature", Tags::DEFAULT)
       .SetMesh(mesh_matrix_)
       ->SetGhosted(true) = *cvs;
    AddDefaultPrimaryEvaluator(S_, "temperature", Tags::DEFAULT);
  }

  // -- molar flow rates
  if (!S_->HasRecord("molar_flow_rate")) {
    auto mmap = cvs->Map("face", false);
    auto gmap = cvs->Map("face", true);
    S_->Require<CV_t, CVS_t>("molar_flow_rate", Tags::DEFAULT)
      .SetMesh(mesh_matrix_)
      ->SetGhosted(true)
      ->SetComponent("face", AmanziMesh::FACE, mmap, gmap, 1);
    AddDefaultPrimaryEvaluator(S_, "molar_flow_rate", Tags::DEFAULT);
  }

  if (!S_->HasRecord("fracture-molar_flow_rate")) {
    auto cvs2 = Operators::CreateManifoldCVS(mesh_fracture_);
    *S_->Require<CV_t, CVS_t>("fracture-molar_flow_rate", Tags::DEFAULT)
       .SetMesh(mesh_fracture_)
       ->SetGhosted(true) = *cvs2;
    AddDefaultPrimaryEvaluator(S_, "fracture-molar_flow_rate", Tags::DEFAULT);
  }

  // additional fields and evaluators related to matrix-frcature coupling
  if (!S_->HasRecord(heat_diffusion_to_matrix_key_)) {
    S_->Require<CV_t, CVS_t>(
        heat_diffusion_to_matrix_key_, Tags::DEFAULT, heat_diffusion_to_matrix_key_)
      .SetMesh(mesh_fracture_)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 2);
    S_->RequireEvaluator(heat_diffusion_to_matrix_key_, Tags::DEFAULT);
  }

  // inform dependent PKs about coupling
  std::vector<std::string> pks = plist_->get<Teuchos::Array<std::string>>("PKs order").toVector();

  glist_->sublist("PKs")
    .sublist(pks[0])
    .sublist("physical models and assumptions")
    .set<std::string>("coupled matrix fracture energy", "matrix");

  glist_->sublist("PKs")
    .sublist(pks[1])
    .sublist("physical models and assumptions")
    .set<std::string>("coupled matrix fracture energy", "fracture");

  // process other PKs.
  PK_MPCStrong<PK_BDF>::Setup();
}


/* *******************************************************************
* Initialization create a tree operator to assemble global matrix
******************************************************************* */
void
EnergyMatrixFracture_PK::Initialize()
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
  auto eval = S_->GetEvaluatorPtr(heat_diffusion_to_matrix_key_, Tags::DEFAULT);
  auto eval_tmp = Teuchos::rcp_dynamic_cast<HeatDiffusionMatrixFracture>(eval);
  if (eval_tmp.get()) eval_tmp->Update(*S_, "coupled energy");

  const auto& kn = *S_->Get<CV_t>(heat_diffusion_to_matrix_key_).ViewComponent("cell");

  FractureInsertion fi(mesh_matrix_, mesh_fracture_);
  fi.InitMatrixFaceToFractureCell(Teuchos::rcpFromRef(mmap), Teuchos::rcpFromRef(gmap));
  fi.SetValues(kn, 1.0);

  // -- operators. They must have own pointers to global operators to be used as
  // building blocks later.
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

  // -- operators for enthalpy. They must have own pointers to global operators
  // to be used as building blocks later.
  FractureInsertion fia(mesh_matrix_, mesh_fracture_);
  fia.InitMatrixCellToFractureCell();

  auto op_enthalpy00 = Teuchos::rcp(new Operators::PDE_CouplingFlux(oplist,
                                                                    fia.get_cvs_matrix(),
                                                                    fia.get_cvs_matrix(),
                                                                    fia.get_inds_matrix(),
                                                                    fia.get_inds_matrix()));
  op0->OpPushBack(op_enthalpy00->local_op());
  op_coupling00->global_operator()->OpPushBack(op_enthalpy00->local_op());
  op_enthalpy00->Setup(fia.get_values(), 1.0);
  op_enthalpy00->UpdateMatrices(Teuchos::null, Teuchos::null);

  auto op_enthalpy01 = Teuchos::rcp(new Operators::PDE_CouplingFlux(oplist,
                                                                    fia.get_cvs_matrix(),
                                                                    fia.get_cvs_fracture(),
                                                                    fia.get_inds_matrix(),
                                                                    fia.get_inds_fracture()));
  op_coupling01->global_operator()->OpPushBack(op_enthalpy01->local_op());
  op_enthalpy01->Setup(fia.get_values(), -1.0);
  op_enthalpy01->UpdateMatrices(Teuchos::null, Teuchos::null);

  auto op_enthalpy10 = Teuchos::rcp(new Operators::PDE_CouplingFlux(oplist,
                                                                    fia.get_cvs_fracture(),
                                                                    fia.get_cvs_matrix(),
                                                                    fia.get_inds_fracture(),
                                                                    fia.get_inds_matrix()));
  op_coupling10->global_operator()->OpPushBack(op_enthalpy10->local_op());
  op_enthalpy10->Setup(fia.get_values(), -1.0);
  op_enthalpy10->UpdateMatrices(Teuchos::null, Teuchos::null);

  auto op_enthalpy11 = Teuchos::rcp(new Operators::PDE_CouplingFlux(oplist,
                                                                    fia.get_cvs_fracture(),
                                                                    fia.get_cvs_fracture(),
                                                                    fia.get_inds_fracture(),
                                                                    fia.get_inds_fracture()));
  op1->OpPushBack(op_enthalpy11->local_op());
  op_coupling11->global_operator()->OpPushBack(op_enthalpy11->local_op());
  op_enthalpy11->Setup(fia.get_values(), 1.0);
  op_enthalpy11->UpdateMatrices(Teuchos::null, Teuchos::null);

  adv_coupling_ops_.push_back(op_enthalpy00);
  adv_coupling_ops_.push_back(op_enthalpy01);
  adv_coupling_ops_.push_back(op_enthalpy10);
  adv_coupling_ops_.push_back(op_enthalpy11);

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
EnergyMatrixFracture_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  // create copies of conservative fields
  std::vector<std::string> fields = { "energy", "fracture-energy" };
  StateArchive archive(S_, vo_);
  archive.Add(fields, Tags::DEFAULT);

  bool fail = PK_MPCStrong<PK_BDF>::AdvanceStep(t_old, t_new, reinit);

  if (fail) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "Step failed." << std::endl;

    archive.Restore("");
  }

  return fail;
}


/* *******************************************************************
* Residual evaluation
******************************************************************* */
void
EnergyMatrixFracture_PK::FunctionalResidual(double t_old,
                                            double t_new,
                                            Teuchos::RCP<TreeVector> u_old,
                                            Teuchos::RCP<TreeVector> u_new,
                                            Teuchos::RCP<TreeVector> f)
{
  PK_MPCStrong<PK_BDF>::FunctionalResidual(t_old, t_new, u_old, u_new, f);

  UpdateEnthalpyCouplingFluxes(S_, mesh_matrix_, mesh_fracture_, adv_coupling_ops_, true);

  int ierr = op_matrix_->Apply(*u_new, *f, 1.0);
  AMANZI_ASSERT(!ierr);
}


/* *******************************************************************
* Preconditioner update
******************************************************************* */
void
EnergyMatrixFracture_PK::UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> u, double dt)
{
  PK_MPCStrong<PK_BDF>::UpdatePreconditioner(t, u, dt);
  op_tree_pc_->ComputeInverse();
}


/* *******************************************************************
* Application of preconditioner
******************************************************************* */
int
EnergyMatrixFracture_PK::ApplyPreconditioner(Teuchos::RCP<const TreeVector> X,
                                             Teuchos::RCP<TreeVector> Y)
{
  Y->PutScalar(0.0);
  return op_tree_pc_->ApplyInverse(*X, *Y);
}

} // namespace Amanzi
