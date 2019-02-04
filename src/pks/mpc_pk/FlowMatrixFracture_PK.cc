/*
  This is the mpc_pk component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
           Daniil Svyatskiy

  Process kernel that couples Flow flow in matrix and fractures.
*/

#include "Darcy_PK.hh"
#include "PDE_CouplingFlux.hh"
#include "PDE_DiffusionFracturedMatrix.hh"
#include "primary_variable_field_evaluator.hh"
#include "TreeOperator.hh"

#include "FlowMatrixFracture_PK.hh"
#include "PK_MPCStrong.hh"

namespace Amanzi {

/* ******************************************************************* 
* Constructor
******************************************************************* */
FlowMatrixFracture_PK::FlowMatrixFracture_PK(Teuchos::ParameterList& pk_tree,
                                             const Teuchos::RCP<Teuchos::ParameterList>& glist,
                                             const Teuchos::RCP<State>& S,
                                             const Teuchos::RCP<TreeVector>& soln) :
  glist_(glist), matrix_assembled_(false),
  Amanzi::PK_MPC<PK_BDF>(pk_tree, glist, S, soln),
  Amanzi::PK_MPCStrong<PK_BDF>(pk_tree, glist, S, soln)
{
  Teuchos::ParameterList vlist;
  vo_ =  Teuchos::rcp(new VerboseObject("MatrixFracture_PK", vlist));
  Teuchos::RCP<Teuchos::ParameterList> pks_list = Teuchos::sublist(glist, "PKs");
  if (pks_list->isSublist(name_)) {
    plist_ = Teuchos::sublist(pks_list, name_); 
  } else {
    std::stringstream messagestream;
    messagestream << "There is no sublist for PK "<<name_<<"in PKs list\n";
    Errors::Message message(messagestream.str());
    Exceptions::amanzi_throw(message);
  }
  
  dump_ = plist_->get<bool>("dump preconditioner", false);

  preconditioner_list_ = Teuchos::sublist(glist, "preconditioners", true);
  linear_operator_list_ = Teuchos::sublist(glist, "solvers", true);

}


/* ******************************************************************* 
* Physics-based setup of PK.
******************************************************************* */
void FlowMatrixFracture_PK::Setup(const Teuchos::Ptr<State>& S)
{
  mesh_domain_ = S->GetMesh();
  mesh_fracture_ = S->GetMesh("fracture");
  int dim = mesh_domain_->space_dimension();

  Teuchos::ParameterList& elist = S->FEList();

  // primary and secondary fields for matrix affected by non-uniform
  // distribution of DOFs
  // -- pressure
  auto cvs = Operators::CreateFracturedMatrixCVS(mesh_domain_, mesh_fracture_);
  if (!S->HasField("pressure")) {
    *S->RequireField("pressure", "flow")->SetMesh(mesh_domain_)->SetGhosted(true) = *cvs;
  }

  // -- darcy flux
  if (!S->HasField("darcy_flux")) {
    std::string name("face");
    auto mmap = cvs->Map("face", false);
    auto gmap = cvs->Map("face", true);
    S->RequireField("darcy_flux", "flow")->SetMesh(mesh_domain_)->SetGhosted(true) 
      ->SetComponent(name, mmap, gmap, 1);
  }

  // Require additional fields and evaluators

  // inform dependent PKs about coupling
  // -- flow (matrix)
  Teuchos::ParameterList& mflow = glist_->sublist("PKs").sublist("flow matrix")
                                         .sublist("Darcy problem")
                                         .sublist("physical models and assumptions");
  mflow.set<std::string>("coupled matrix fracture flow", "matrix");

  // -- flow (fracture)
  Teuchos::ParameterList& fflow = glist_->sublist("PKs").sublist("flow fracture")
                                         .sublist("Darcy problem")
                                         .sublist("physical models and assumptions");
  fflow.set<std::string>("coupled matrix fracture flow", "fracture");

  // process other PKs.
  PK_MPCStrong<PK_BDF>::Setup(S);
}


/* ******************************************************************* 
* Initialization create a tree operator to assemble global matrix
******************************************************************* */
void FlowMatrixFracture_PK::Initialize(const Teuchos::Ptr<State>& S)
{
  PK_MPCStrong<PK_BDF>::Initialize(S);

  // diagonal blocks in tree operator and the Darcy PKs
  auto pk_matrix = Teuchos::rcp_dynamic_cast<Flow::Darcy_PK>(sub_pks_[0]);
  auto pk_fracture = Teuchos::rcp_dynamic_cast<Flow::Darcy_PK>(sub_pks_[1]);

  auto tvs = Teuchos::rcp(new TreeVectorSpace(solution_->Map()));
  op_tree_ = Teuchos::rcp(new Operators::TreeOperator(tvs));
  //  op_tree_rhs_ =  Teuchos::rcp(new TreeVector(tvs));

  op_tree_->SetOperatorBlock(0, 0, pk_matrix->op());
  op_tree_->SetOperatorBlock(1, 1, pk_fracture->op());

  //op_tree_rhs_->PushBack(pk_matrix->op()->rhs());
  //op_tree_rhs_->PushBack(pk_fracture->op()->rhs());

  // off-diagonal blocks are coupled PDEs
  // -- minimum composite vector spaces containing the coupling term
  auto mesh_matrix = S_->GetMesh("domain");
  auto mesh_fracture = S_->GetMesh("fracture");

  auto& mmap = solution_->SubVector(0)->Data()->ViewComponent("face", false)->Map();
  auto& gmap = solution_->SubVector(0)->Data()->ViewComponent("face", true)->Map();
  int npoints_owned = mmap.NumMyPoints();

  auto cvs_matrix = Teuchos::rcp(new CompositeVectorSpace());
  auto cvs_fracture = Teuchos::rcp(new CompositeVectorSpace());

  std::string compname("face");
  cvs_matrix->SetMesh(mesh_matrix)->SetGhosted(true)
            ->AddComponent(compname, Teuchos::rcpFromRef(mmap), Teuchos::rcpFromRef(gmap), 1);

  cvs_fracture->SetMesh(mesh_matrix)->SetGhosted(true)
              ->AddComponent("cell", AmanziMesh::CELL, 1);

  // -- indices transmissibimility coefficients for matrix-fracture flux
  auto kn = *S_->GetFieldData("fracture-normal_permeability")->ViewComponent("cell");

  int ncells_owned_f = mesh_fracture->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  auto inds_matrix = std::make_shared<std::vector<std::vector<int> > >(npoints_owned);
  auto inds_fracture = std::make_shared<std::vector<std::vector<int> > >(npoints_owned);
  auto values = std::make_shared<std::vector<double> >(npoints_owned);

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

      (*values)[np] = kn[0][c] * area;
      np++;
    }
  }

  // -- operators
  Teuchos::ParameterList oplist;

  auto op_coupling00 = Teuchos::rcp(new Operators::PDE_CouplingFlux(
      oplist, cvs_matrix, cvs_matrix, inds_matrix, inds_matrix, pk_matrix->op()));
  op_coupling00->Setup(values, 1.0);

  auto op_coupling01 = Teuchos::rcp(new Operators::PDE_CouplingFlux(
      oplist, cvs_matrix, cvs_fracture, inds_matrix, inds_fracture));
  op_coupling01->Setup(values, -1.0);

  auto op_coupling10 = Teuchos::rcp(new Operators::PDE_CouplingFlux(
      oplist, cvs_fracture, cvs_matrix, inds_fracture, inds_matrix));
  op_coupling10->Setup(values, -1.0);

  auto op_coupling11 = Teuchos::rcp(new Operators::PDE_CouplingFlux(
      oplist, cvs_fracture, cvs_fracture, inds_fracture, inds_fracture, pk_fracture->op()));
  op_coupling11->Setup(values, 1.0);

  op_tree_->SetOperatorBlock(0, 1, op_coupling01->global_operator());
  op_tree_->SetOperatorBlock(1, 0, op_coupling10->global_operator());

  op_tree_->set_multi_domain(true);

  op_tree_->SymbolicAssembleMatrix();

  op_tree_->AssembleMatrix();

  matrix_assembled_ = true;

}


/* ******************************************************************* 
* Performs one time step.
******************************************************************* */
bool FlowMatrixFracture_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  bool fail = PK_MPCStrong<PK_BDF>::AdvanceStep(t_old, t_new, reinit);

  if (fail) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "Step failed." << std::endl;
  }

  return fail;
}

void FlowMatrixFracture_PK::FunctionalResidual(double t_old, double t_new,
                                               Teuchos::RCP<TreeVector> u_old,
                                               Teuchos::RCP<TreeVector> u_new,
                                               Teuchos::RCP<TreeVector> f){
  double norm, sol;
  
  int ierr = op_tree_->ApplyAssembled(*u_new, *f);
  AMANZI_ASSERT(!ierr);
  
  // diagonal blocks in tree operator and the Darcy PKs
  auto pk_matrix = Teuchos::rcp_dynamic_cast<Flow::Darcy_PK>(sub_pks_[0]);
  auto pk_fracture = Teuchos::rcp_dynamic_cast<Flow::Darcy_PK>(sub_pks_[1]);

  f->SubVector(0)->Data()->Update(-1, *pk_matrix->op()->rhs(), 1);
  f->SubVector(1)->Data()->Update(-1, *pk_fracture->op()->rhs(), 1);

  

}

void FlowMatrixFracture_PK::UpdatePreconditioner(double t,
                                                 Teuchos::RCP<const TreeVector> up,
                                                 double h, bool assemble) {

  if (assemble) {
    op_tree_->AssembleMatrix();
    if (dump_) {
      std::stringstream filename;
      filename << "FlowMatrixFracture_PC_.txt";
      EpetraExt::RowMatrixToMatlabFile(filename.str().c_str(), *op_tree_->A());
    }
    
    std::string name = plist_->get<std::string>("preconditioner");
    Teuchos::ParameterList pc_list = preconditioner_list_->sublist(name);

    op_tree_ -> InitPreconditioner(pc_list);
  }
  
}


int FlowMatrixFracture_PK::ApplyPreconditioner(Teuchos::RCP<const TreeVector> X, 
                             Teuchos::RCP<TreeVector> Y){

  Y->PutScalar(0.0);
  return op_tree_->ApplyInverse(*X, *Y);

}


}  // namespace Amanzi

