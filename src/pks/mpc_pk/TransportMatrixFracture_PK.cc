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

  Process kernel that couples Transport in matrix and fracture.

  For molecular diffusion, we use the boundary layer or film theory (Fahien, 1983),
  according to which the mass transfer of solute per unit area at the fracture-matrix
  interface is given by kn (Cf - Cm). The mass transfer coefficient kn = Dm / h,
  where h is the effective layer thickness.
*/

// Amanzi
#include "FractureInsertion.hh"
#include "PDE_CouplingFlux.hh"
#include "PK_MPCStrong.hh"
#include "InverseFactory.hh"
#include "Transport_PK.hh"

#include "SoluteDiffusionMatrixFracture.hh"
#include "TransportMatrixFracture_PK.hh"

namespace Amanzi {

using CV_t = CompositeVector;
using CVS_t = CompositeVectorSpace;

/* *******************************************************************
* Constructor
******************************************************************* */
TransportMatrixFracture_PK::TransportMatrixFracture_PK(
  Teuchos::ParameterList& pk_tree,
  const Teuchos::RCP<Teuchos::ParameterList>& glist,
  const Teuchos::RCP<State>& S,
  const Teuchos::RCP<TreeVector>& soln)
  : Amanzi::PK(pk_tree, glist, S, soln), Amanzi::PK_MPCWeak(pk_tree, glist, S, soln), glist_(glist)
{
  Teuchos::RCP<Teuchos::ParameterList> pks_list = Teuchos::sublist(glist, "PKs");
  plist_ = Teuchos::sublist(pks_list, name_);

  Teuchos::ParameterList vlist;
  vlist.sublist("verbose object") = plist_->sublist("verbose object");
  vo_ = Teuchos::rcp(new VerboseObject("TransportMatrixFracture_PK", vlist));
}


/* *******************************************************************
* Physics-based setup of PK.
******************************************************************* */
void
TransportMatrixFracture_PK::Setup()
{
  mesh_domain_ = S_->GetMesh();
  mesh_fracture_ = S_->GetMesh("fracture");

  tcc_matrix_key_ = "total_component_concentration";
  tcc_fracture_key_ = "fracture-total_component_concentration";
  solute_diffusion_to_matrix_key_ = "fracture-solute_diffusion_to_matrix";

  // darcy fluxes use non-uniform distribution of DOFs
  // -- darcy flux for matrix
  if (!S_->HasRecord("volumetric_flow_rate")) {
    auto cvs = Operators::CreateFracturedMatrixCVS(mesh_domain_, mesh_fracture_);
    auto mmap = cvs->Map("face", false);
    auto gmap = cvs->Map("face", true);
    S_->Require<CV_t, CVS_t>("volumetric_flow_rate", Tags::DEFAULT, "transport")
      .SetMesh(mesh_domain_)
      ->SetGhosted(true)
      ->SetComponent("face", AmanziMesh::FACE, mmap, gmap, 1);
  }

  // -- darcy flux for fracture
  if (!S_->HasRecord("fracture-volumetric_flow_rate")) {
    auto cvs = Operators::CreateManifoldCVS(mesh_fracture_);
    *S_->Require<CV_t, CVS_t>("fracture-volumetric_flow_rate", Tags::DEFAULT, "transport")
       .SetMesh(mesh_fracture_)
       ->SetGhosted(true) = *cvs;
  }

  // Require additional fields and evaluators
  if (!S_->HasRecord(solute_diffusion_to_matrix_key_)) {
    S_->Require<CV_t, CVS_t>(
        solute_diffusion_to_matrix_key_, Tags::DEFAULT, solute_diffusion_to_matrix_key_)
      .SetMesh(mesh_fracture_)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 2);
    S_->RequireEvaluator(solute_diffusion_to_matrix_key_, Tags::DEFAULT);
  }

  // add boundary condition to transport in matrix list
  auto pks =
    glist_->sublist("PKs").sublist(name_).get<Teuchos::Array<std::string>>("PKs order").toVector();
  Teuchos::ParameterList& bclist = glist_->sublist("PKs")
                                     .sublist(pks[0])
                                     .sublist("boundary conditions")
                                     .sublist("concentration")
                                     .sublist("coupling")
                                     .sublist("BC coupling");

  Teuchos::Array<std::string> regs;
  regs.push_back("FRACTURE_NETWORK_INTERNAL");
  bclist.set<std::string>("spatial distribution method", "domain coupling")
    .set<std::string>("submodel", "field")
    .set<Teuchos::Array<std::string>>("regions", regs);

  bclist.sublist("boundary concentration")
    .set<std::string>("external field key", tcc_fracture_key_);

  // add source term to transport in fracture list
  Teuchos::ParameterList& srclist = glist_->sublist("PKs")
                                      .sublist(pks[1])
                                      .sublist("source terms")
                                      .sublist("concentration")
                                      .sublist("coupling")
                                      .sublist("fracture");

  regs.clear();
  regs.push_back("All");
  srclist.set<std::string>("spatial distribution method", "domain coupling")
    .set<std::string>("submodel", "rate")
    .set<Teuchos::Array<std::string>>("regions", regs);

  srclist.sublist("sink")
    .set<std::string>("external field key", tcc_matrix_key_)
    .set<std::string>("flux key", "volumetric_flow_rate");

  // setup the sub-PKs
  PK_MPCWeak::Setup();
}


/* *******************************************************************
* Reduce stable dt to avoid the 2-cycle behavior of the transport PK
******************************************************************* */
double
TransportMatrixFracture_PK::get_dt()
{
  return 0.999 * PK_MPCWeak::get_dt();
}


/* *******************************************************************
* Initialization create a tree operator to assemble global matrix
******************************************************************* */
void
TransportMatrixFracture_PK::Initialize()
{
  PK_MPCWeak::Initialize();

  flag_dispersion_ = false;
  for (const auto& pk : sub_pks_) {
    flag_dispersion_ |=
      Teuchos::rcp_dynamic_cast<Transport::Transport_PK>(pk)->get_flag_dispersion();
  }

  // decision to create this field could be possibly made during setup
  if (!flag_dispersion_)
    InitializeCVField(S_, *vo_, solute_diffusion_to_matrix_key_, Tags::DEFAULT, "state", 0.0);
}


/* *******************************************************************
* Performs one time step.
******************************************************************* */
bool
TransportMatrixFracture_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  bool fail = PK_MPCWeak::AdvanceStep(t_old, t_new, reinit);

  if (fail) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "Step failed." << std::endl;
  }

  if (!flag_dispersion_) return fail;

  std::string passwd("state");

  // we assume that 0 and 1 correspond to matrix and fracture
  auto& tcc_prev_m =
    *S_->GetW<CV_t>(tcc_matrix_key_, Tags::DEFAULT, passwd).ViewComponent("cell", true);
  auto& tcc_prev_f =
    *S_->GetW<CV_t>(tcc_fracture_key_, Tags::DEFAULT, passwd).ViewComponent("cell", true);

  auto& tcc_next_m =
    *S_->GetW<CV_t>(tcc_matrix_key_, Tags::COPY, passwd).ViewComponent("cell", true);
  auto& tcc_next_f =
    *S_->GetW<CV_t>(tcc_fracture_key_, Tags::COPY, passwd).ViewComponent("cell", true);

  auto pk0 = Teuchos::rcp_dynamic_cast<Transport::Transport_PK>(sub_pks_[0]);
  auto pk1 = Teuchos::rcp_dynamic_cast<Transport::Transport_PK>(sub_pks_[1]);

  int num_aqueous = tcc_next_m.NumVectors();

  for (int i = 0; i < num_aqueous; ++i) {
    auto op0 = pk0->DispersionSolver(tcc_prev_m, tcc_next_m, t_old, t_new, i)->Clone();
    auto op1 = pk1->DispersionSolver(tcc_prev_f, tcc_next_f, t_old, t_new, i)->Clone();

    auto eval = S_->GetEvaluatorPtr("fracture-solute_diffusion_to_matrix", Tags::DEFAULT);
    auto eval_tmp = Teuchos::rcp_dynamic_cast<SoluteDiffusionMatrixFracture>(eval);
    if (eval_tmp.get()) {
      eval_tmp->set_subvector(pk0->getDiffusion(i));
      eval_tmp->Update(*S_, "coupled transport");
    }

    // since solution's map could be anything, to create a global operator,
    // we have to rely on pk's operator structure.
    auto& cvs0 = op0->get_domain_map();
    auto& cvs1 = op1->get_domain_map();

    auto tvs = Teuchos::rcp(new TreeVectorSpace());
    tvs->PushBack(Teuchos::rcp(new TreeVectorSpace(cvs0)));
    tvs->PushBack(Teuchos::rcp(new TreeVectorSpace(cvs1)));

    // primary dispersion/diffusion operators
    op_dispersion_ = Teuchos::rcp(new Operators::TreeOperator(tvs));
    op_dispersion_->set_operator_block(0, 0, op0);
    op_dispersion_->set_operator_block(1, 1, op1);

    // coupling operators
    // -- indices transmissibimility coefficients for matrix-fracture flux
    auto mesh_matrix = S_->GetMesh("domain");
    auto mesh_fracture = S_->GetMesh("fracture");
    const auto& kn = *S_->Get<CV_t>("fracture-solute_diffusion_to_matrix").ViewComponent("cell");

    auto& mmap = *cvs0->Map("face", false);
    auto& gmap = *cvs0->Map("face", true);

    FractureInsertion fi(mesh_matrix, mesh_fracture);
    fi.InitMatrixFaceToFractureCell(Teuchos::rcpFromRef(mmap), Teuchos::rcpFromRef(gmap));
    fi.SetValues(kn, 1.0);

    // -- generate add interface operators to the tree operator
    Teuchos::ParameterList oplist;

    auto op_coupling00 = Teuchos::rcp(new Operators::PDE_CouplingFlux(oplist,
                                                                      fi.get_cvs_matrix(),
                                                                      fi.get_cvs_matrix(),
                                                                      fi.get_inds_matrix(),
                                                                      fi.get_inds_matrix(),
                                                                      op0));
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
                                                                      op1));
    op_coupling11->Setup(fi.get_values(), 1.0);
    op_coupling11->UpdateMatrices(Teuchos::null, Teuchos::null);

    op_dispersion_->set_operator_block(0, 1, op_coupling01->global_operator());
    op_dispersion_->set_operator_block(1, 0, op_coupling10->global_operator());

    // assemple and solve
    std::string pc_name = plist_->get<std::string>("preconditioner");
    std::string ls_name = plist_->get<std::string>("solver", "none");
    auto inv_list = AmanziSolvers::mergePreconditionerSolverLists(
      pc_name, glist_->sublist("preconditioners"), ls_name, glist_->sublist("solvers"), true);
    inv_list.setName(pc_name);

    op_dispersion_->set_inverse_parameters(inv_list);
    op_dispersion_->InitializeInverse();
    op_dispersion_->ComputeInverse();

    TreeVector rhs(tvs), sol(tvs);
    *rhs.SubVector(0)->Data() = *op0->rhs();
    *rhs.SubVector(1)->Data() = *op1->rhs();

    op_dispersion_->ApplyInverse(rhs, sol);

    // copy only cell component from potentially larger solution vector
    auto& tcc_m = *S_->GetW<CV_t>(tcc_matrix_key_, Tags::COPY, passwd).ViewComponent("cell");
    *tcc_m(i) = *(*sol.SubVector(0)->Data()->ViewComponent("cell"))(0);

    auto& tcc_f = *S_->GetW<CV_t>(tcc_fracture_key_, Tags::COPY, passwd).ViewComponent("cell");
    *tcc_f(i) = *(*sol.SubVector(1)->Data()->ViewComponent("cell"))(0);

    if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
      Teuchos::OSTab tab = vo_->getOSTab();
      pk0->VV_PrintSoluteExtrema(tcc_m, t_new - t_old, " (domain)");
      pk1->VV_PrintSoluteExtrema(tcc_f, t_new - t_old, " (fracture)");
    }
  }

  return fail;
}

} // namespace Amanzi
