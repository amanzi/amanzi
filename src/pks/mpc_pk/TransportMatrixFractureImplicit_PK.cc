/*
  This is the mpc_pk component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
           Daniil Svyatskiy

  Process kernel that couples Transport in matrix and fracture
  using implicit scheme.
*/

#include "LinearOperatorPCG.hh"
#include "Op_Diagonal.hh"
#include "PDE_CouplingFlux.hh"
#include "PDE_DiffusionFracturedMatrix.hh"
#include "Transport_PK.hh"
#include "TreeOperator.hh"

#include "TransportMatrixFractureImplicit_PK.hh"

namespace Amanzi {

/* ******************************************************************* 
* Constructor
******************************************************************* */
TransportMatrixFractureImplicit_PK::TransportMatrixFractureImplicit_PK(
    Teuchos::ParameterList& pk_tree,
    const Teuchos::RCP<Teuchos::ParameterList>& glist,
    const Teuchos::RCP<State>& S,
    const Teuchos::RCP<TreeVector>& soln)
  : glist_(glist),
    Amanzi::PK_MPC<PK_BDF>(pk_tree, glist, S, soln),
    Amanzi::PK_MPCStrong<PK_BDF>(pk_tree, glist, S, soln)
{
  Teuchos::ParameterList vlist;
  vo_ = Teuchos::rcp(new VerboseObject("TransportMatrixFractureImplicit_PK", vlist)); 
}


/* ******************************************************************* 
* Physics-based setup of PK.
******************************************************************* */
void TransportMatrixFractureImplicit_PK::Setup(const Teuchos::Ptr<State>& S)
{
  mesh_domain_ = S->GetMesh();
  mesh_fracture_ = S->GetMesh("fracture");

  // primary and secondary fields for matrix affected by non-uniform
  // distribution of DOFs
  auto cvs = Operators::CreateFracturedMatrixCVS(mesh_domain_, mesh_fracture_);

  // -- darcy flux
  if (!S->HasField("darcy_flux")) {
    auto mmap = cvs->Map("face", false);
    auto gmap = cvs->Map("face", true);
    S->RequireField("darcy_flux", "transport")->SetMesh(mesh_domain_)->SetGhosted(true) 
      ->SetComponent("face", AmanziMesh::FACE, mmap, gmap, 1);
  }

  PK_MPCStrong::Setup(S);

  // diagonal blocks in tree operator are the Transport Implicit PKs
  auto pk_matrix = Teuchos::rcp_dynamic_cast<Transport::TransportImplicit_PK>(sub_pks_[0]);
  auto pk_fracture = Teuchos::rcp_dynamic_cast<Transport::TransportImplicit_PK>(sub_pks_[1]);

  auto tvs = Teuchos::rcp(new TreeVectorSpace(solution_->Map()));
  op_tree_ = Teuchos::rcp(new Operators::TreeOperator(tvs));

  op_tree_->SetOperatorBlock(0, 0, pk_matrix->op());
  op_tree_->SetOperatorBlock(1, 1, pk_fracture->op());

  // off-diagonal blocks are coupled PDEs
  // -- minimum composite vector spaces containing the coupling term
  auto mesh_matrix = S_->GetMesh("domain");
  auto mesh_fracture = S_->GetMesh("fracture");

  auto& mmap = solution_->SubVector(0)->Data()->ViewComponent("face", false)->Map();
  auto& gmap = solution_->SubVector(0)->Data()->ViewComponent("face", true)->Map();
  int npoints_owned = mmap.NumMyPoints();

  auto cvs_matrix = Teuchos::rcp(new CompositeVectorSpace());
  auto cvs_fracture = Teuchos::rcp(new CompositeVectorSpace());

  cvs_matrix->SetMesh(mesh_matrix)->SetGhosted(true)
            ->AddComponent("face", AmanziMesh::FACE, Teuchos::rcpFromRef(mmap), Teuchos::rcpFromRef(gmap), 1);

  cvs_fracture->SetMesh(mesh_matrix)->SetGhosted(true)
              ->AddComponent("cell", AmanziMesh::CELL, 1);

  // -- indices transmissibimility coefficients for matrix-fracture flux
  const auto& flux = *S_->GetFieldData("darcy_flux")->ViewComponent("face");

  int ncells_owned_f = mesh_fracture->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  auto inds_matrix = std::make_shared<std::vector<std::vector<int> > >(npoints_owned);
  auto inds_fracture = std::make_shared<std::vector<std::vector<int> > >(npoints_owned);
  auto values = std::make_shared<std::vector<double> >(npoints_owned, 0.0);

  int np(0);
  for (int c = 0; c < ncells_owned_f; ++c) {
    int f = mesh_fracture->entity_get_parent(AmanziMesh::CELL, c);
    double area = mesh_fracture->cell_volume(c);
    int first = mmap.FirstPointInElement(f);
    int ndofs = mmap.ElementSize(f);

    for (int k = 0; k < ndofs; ++k) {
      (*inds_matrix)[np].resize(1);
      (*inds_fracture)[np].resize(1);
      (*inds_matrix)[np][0] = first + k;
      (*inds_fracture)[np][0] = c;

      (*values)[np] = flux[0][f] * area;
      np++;
    }
  }

  inds_matrix->resize(np);
  inds_fracture->resize(np);
  values->resize(np);

  // -- operators
  Teuchos::ParameterList oplist;

  op_coupling00_ = Teuchos::rcp(new Operators::PDE_CouplingFlux(
      oplist, cvs_matrix, cvs_matrix, inds_matrix, inds_matrix, pk_matrix->op()));
  op_coupling00_->Setup(values, 1.0);
  op_coupling00_->UpdateMatrices(Teuchos::null, Teuchos::null);

  op_coupling01_ = Teuchos::rcp(new Operators::PDE_CouplingFlux(
      oplist, cvs_matrix, cvs_fracture, inds_matrix, inds_fracture));
  op_coupling01_->Setup(values, -1.0);
  op_coupling01_->UpdateMatrices(Teuchos::null, Teuchos::null);

  op_coupling10_ = Teuchos::rcp(new Operators::PDE_CouplingFlux(
      oplist, cvs_fracture, cvs_matrix, inds_fracture, inds_matrix));
  op_coupling10_->Setup(values, -1.0);
  op_coupling10_->UpdateMatrices(Teuchos::null, Teuchos::null);

  op_coupling11_ = Teuchos::rcp(new Operators::PDE_CouplingFlux(
      oplist, cvs_fracture, cvs_fracture, inds_fracture, inds_fracture, pk_fracture->op()));
  op_coupling11_->Setup(values, 1.0);
  op_coupling11_->UpdateMatrices(Teuchos::null, Teuchos::null);  

  op_tree_->SetOperatorBlock(0, 1, op_coupling01_->global_operator());
  op_tree_->SetOperatorBlock(1, 0, op_coupling10_->global_operator());

  // create a global problem
  pk_matrix->op_adv()->ApplyBCs(true, true, true);

  op_tree_->SymbolicAssembleMatrix();

  // Test SPD properties of the matrix.
  // VerificationTV ver(op_tree_);
  // ver.CheckMatrixSPD();
}


/* ******************************************************************* 
* Initialization creates a tree operator to assemble global matrix
******************************************************************* */
void TransportMatrixFractureImplicit_PK::Initialize(const Teuchos::Ptr<State>& S)
{
  PK_MPCStrong<PK_BDF>::Initialize(S);

  auto tvs = Teuchos::rcp(new TreeVectorSpace(solution_->Map()));
  op_tree_ = Teuchos::rcp(new Operators::TreeOperator(tvs));
}


/* ******************************************************************* 
* Performs one time step.
******************************************************************* */
bool TransportMatrixFractureImplicit_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  double dt = t_new - t_old;

  auto pk_matrix = Teuchos::rcp_dynamic_cast<Transport::TransportImplicit_PK>(sub_pks_[0]);
  auto pk_fracture = Teuchos::rcp_dynamic_cast<Transport::TransportImplicit_PK>(sub_pks_[1]);

  // update coupling terms
  int np = Teuchos::rcp_dynamic_cast<Operators::Op_Diagonal>(op_coupling00_->local_op())->row_inds().size();
  auto values1 = std::make_shared<std::vector<double> >(np, 0.0);
  auto values2 = std::make_shared<std::vector<double> >(np, 0.0);

  int ncells_owned_f = mesh_fracture_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  const auto& flux = *S_->GetFieldData("darcy_flux")->ViewComponent("face");
  auto mmap = flux.Map();

  np = 0;
  for (int c = 0; c < ncells_owned_f; ++c) {
    int f = mesh_fracture_->entity_get_parent(AmanziMesh::CELL, c);
    double area = mesh_fracture_->cell_volume(c);
    int first = mmap.FirstPointInElement(f);
    int ndofs = mmap.ElementSize(f);

    for (int k = 0; k < ndofs; ++k) {
      double tmp = flux[0][first + k] * area * dt;
      if (tmp > 0) 
        (*values1)[np] = tmp;
      else 
        (*values2)[np] = tmp;
      np++;
    }
  }

  // assemble the operators
  pk_matrix->op_adv()->ApplyBCs(true, true, true);
  pk_fracture->op_adv()->ApplyBCs(true, true, true);

  op_coupling00_->Setup(values1, 1.0);
  op_coupling00_->UpdateMatrices(Teuchos::null, Teuchos::null);

  op_coupling01_->Setup(values1, -1.0);
  op_coupling01_->UpdateMatrices(Teuchos::null, Teuchos::null);

  op_coupling10_->Setup(values2, -1.0);
  op_coupling10_->UpdateMatrices(Teuchos::null, Teuchos::null);

  op_coupling11_->Setup(values2, 1.0);
  op_coupling11_->UpdateMatrices(Teuchos::null, Teuchos::null);  

  op_tree_->AssembleMatrix();

  // create solver
  Teuchos::ParameterList slist = glist_->sublist("solvers")
                                        .sublist("AztecOO CG").sublist("pcg parameters");
  AmanziSolvers::LinearOperatorPCG<Operators::TreeOperator, TreeVector, TreeVectorSpace>
      solver(op_tree_, op_tree_);
  solver.Init(slist);

  TreeVector rhs(*solution_);
  *rhs.SubVector(0)->Data() = *pk_matrix->op()->rhs();
  *rhs.SubVector(1)->Data() = *pk_fracture->op()->rhs();

  int ierr = solver.ApplyInverse(rhs, *solution_);

  // process error code
  bool fail = (ierr == 0);
  if (fail) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "Step failed." << std::endl;
  }

  return fail;
}

}  // namespace Amanzi

