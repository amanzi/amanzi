/*
  MPC PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
           Daniil Svyatskiy

  Process kernel that couples flow in matrix and fracture network.
*/

#include "CommonDefs.hh"
#include "EvaluatorPrimary.hh"
#include "InverseFactory.hh"
#include "Multiphase_PK.hh"
#include "PDE_CouplingFlux.hh"
#include "PDE_DiffusionFracturedMatrix.hh"
#include "TreeOperator.hh"

#include "FractureInsertion.hh"
#include "MultiphaseMatrixFracture_PK.hh"
#include "PK_MPCStrong.hh"
#include "PK_Utils.hh"

namespace Amanzi {

using CV_t = CompositeVector;
using CVS_t = CompositeVectorSpace;

/* *******************************************************************
* Constructor
******************************************************************* */
MultiphaseMatrixFracture_PK::MultiphaseMatrixFracture_PK(
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
  vo_ = Teuchos::rcp(new VerboseObject("CoupledMultiphase_PK", vlist));
}


/* *******************************************************************
* Physics-based setup of PK.
******************************************************************* */
void
MultiphaseMatrixFracture_PK::Setup()
{
  std::string passwd("");
  mesh_domain_ = S_->GetMesh();
  mesh_fracture_ = S_->GetMesh("fracture");

  // primary and secondary fields for matrix affected by non-uniform
  // distribution of DOFs
  // -- volumetric fluxes
  auto cvs = Operators::CreateFracturedMatrixCVS(mesh_domain_, mesh_fracture_);
  if (!S_->HasRecord("volumetric_flow_rate_liquid")) {
    std::string name("face");
    auto mmap = cvs->Map("face", false);
    auto gmap = cvs->Map("face", true);
    S_->Require<CV_t, CVS_t>("volumetric_flow_rate_liquid", Tags::DEFAULT, passwd)
      .SetMesh(mesh_domain_)
      ->SetGhosted(true)
      ->SetComponent(name, AmanziMesh::FACE, mmap, gmap, 1);
  }

  if (!S_->HasRecord("volumetric_flow_rate_gas")) {
    std::string name("face");
    auto mmap = cvs->Map("face", false);
    auto gmap = cvs->Map("face", true);
    S_->Require<CV_t, CVS_t>("volumetric_flow_rate_gas", Tags::DEFAULT)
      .SetMesh(mesh_domain_)
      ->SetGhosted(true)
      ->SetComponent(name, AmanziMesh::FACE, mmap, gmap, 1);
  }

  auto cvs2 = Operators::CreateManifoldCVS(mesh_fracture_);
  if (!S_->HasRecord("fracture-volumetric_flow_rate_liquid")) {
    *S_->Require<CV_t, CVS_t>("fracture-volumetric_flow_rate_liquid", Tags::DEFAULT)
       .SetMesh(mesh_fracture_)
       ->SetGhosted(true) = *cvs2;
  }

  if (!S_->HasRecord("fracture-volumetric_flow_rate_gas")) {
    *S_->Require<CV_t, CVS_t>("fracture-volumetric_flow_rate_gas", Tags::DEFAULT)
       .SetMesh(mesh_fracture_)
       ->SetGhosted(true) = *cvs2;
  }

  // Require additional fields and evaluators
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

  // process other PKs
  PK_MPCStrong<PK_BDF>::Setup();
}


/* *******************************************************************
* Initialization create a tree operator to assemble global matrix
******************************************************************* */
void
MultiphaseMatrixFracture_PK::Initialize()
{
  PK_MPCStrong<PK_BDF>::Initialize();

  // we assume that 0 and 1 correspond to matrix and fracture
  auto pk_matrix = Teuchos::rcp_dynamic_cast<Multiphase::Multiphase_PK>(sub_pks_[0]);
  auto pk_fracture = Teuchos::rcp_dynamic_cast<Multiphase::Multiphase_PK>(sub_pks_[1]);

  mesh_domain_ = S_->GetMesh("domain");
  mesh_fracture_ = S_->GetMesh("fracture");

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
  auto& mmap = solution_->SubVector(0)->Data()->ViewComponent("face", false)->Map();
  auto& gmap = solution_->SubVector(0)->Data()->ViewComponent("face", true)->Map();

  // -- indices transmissibimility coefficients for matrix-fracture flux
  const auto& kn = *S_->Get<CV_t>("fracture-normal_permeability").ViewComponent("cell");
  double gravity = norm(S_->Get<AmanziGeometry::Point>("gravity"));

  FractureInsertion fi(mesh_domain_, mesh_fracture_);
  fi.InitMatrixFaceToFractureCell(Teuchos::rcpFromRef(mmap), Teuchos::rcpFromRef(gmap));

  // -- operators
  fi.SetValues(kn, 1.0 / gravity);
  /*
  AddCouplingFluxes_(fi.get_cvs_matrix(), fi.get_cvs_fracture(),
                     fi.get_inds_matrix(), fi.get_inds_fracture(),
                     fi.get_values(), 0, op_tree_matrix_);
  */

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
MultiphaseMatrixFracture_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  // create copies of conservative fields
  std::vector<std::string> fields = { "prev_saturation_liquid", "fracture-prev_saturation_liquid" };
  if (sub_pks_[0]->name() == "richards") {
    fields.push_back("prev_water_storage");
    fields.push_back("fracture-prev_water_storage");
  }

  StateArchive archive(S_, vo_);
  archive.Add(fields, {}, {}, Tags::DEFAULT, name());
  archive.Swap("");

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
MultiphaseMatrixFracture_PK::FunctionalResidual(double t_old,
                                                double t_new,
                                                Teuchos::RCP<TreeVector> u_old,
                                                Teuchos::RCP<TreeVector> u_new,
                                                Teuchos::RCP<TreeVector> f)
{
  PK_MPCStrong<PK_BDF>::FunctionalResidual(t_old, t_new, u_old, u_new, f);

  int ierr = op_tree_matrix_->Apply(*u_new, *f, 1.0);
  AMANZI_ASSERT(!ierr);
}


/* *******************************************************************
* Preconditioner update
******************************************************************* */
void
MultiphaseMatrixFracture_PK::UpdatePreconditioner(double t,
                                                  Teuchos::RCP<const TreeVector> u,
                                                  double dt)
{
  PK_MPCStrong<PK_BDF>::UpdatePreconditioner(t, u, dt);
  op_tree_pc_->ComputeInverse();
}


/* *******************************************************************
* Application of preconditioner
******************************************************************* */
int
MultiphaseMatrixFracture_PK::ApplyPreconditioner(Teuchos::RCP<const TreeVector> X,
                                                 Teuchos::RCP<TreeVector> Y)
{
  Y->PutScalar(0.0);
  int ok = op_tree_pc_->ApplyInverse(*X, *Y);
  return ok;
}

} // namespace Amanzi
