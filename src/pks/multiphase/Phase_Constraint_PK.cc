#include "Phase_Constraint_PK.hh"


namespace Amanzi {
namespace Multiphase {

Phase_Constraint_PK::Phase_Constraint_PK(Teuchos::ParameterList plist,
                                          const Teuchos::RCP<State> S):
                  S_(S), plist_(plist), passwd_("state") 
{
  mesh_ = S->GetMesh();
  dim_ = mesh_->space_dimension();
  henry_coef_ = plist_.sublist("EOS").sublist("Component 2").get<double>("coefficient value");
  //std::cout << "Henry coef: " << henry_coef_ << "\n";
  P_vap_ = plist_.sublist("EOS").sublist("Component 1").get<double>("coefficient value");
  //std::cout << "water vapor pressure: " << P_vap_ << "\n";
}

void Phase_Constraint_PK::Initialize() {
  // Initilize various common data depending on mesh and state.
  ncells_owned_ = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  ncells_wghost_ = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);

  nfaces_owned_ = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  nfaces_wghost_ = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  // Allocate memory for boundary data.
  bc_model_.resize(nfaces_wghost_, 0);
  bc_submodel_.resize(nfaces_wghost_, 0);
  bc_value_.resize(nfaces_wghost_, 0.0);
  bc_coef_.resize(nfaces_wghost_, 0.0);
  bc_mixed_.resize(nfaces_wghost_, 0.0);

  op_bc_ = Teuchos::rcp(new Operators::BCs(Operators::OPERATOR_BC_TYPE_FACE, bc_model_, bc_value_, bc_mixed_));

  Teuchos::ParameterList olist_adv = plist_.sublist("Constraints").sublist("operators").sublist("advection operator");
  op1_prec_ = Teuchos::rcp(new Operators::PDE_AdvectionUpwind(olist_adv, mesh_));
  op2_prec_ = Teuchos::rcp(new Operators::PDE_AdvectionUpwind(olist_adv, mesh_));
  op3_prec_ = Teuchos::rcp(new Operators::PDE_AdvectionUpwind(olist_adv, mesh_));
  op4_prec_ = Teuchos::rcp(new Operators::PDE_AdvectionUpwind(olist_adv, mesh_));
  op5_prec_ = Teuchos::rcp(new Operators::PDE_AdvectionUpwind(olist_adv, mesh_));
  op6_prec_ = Teuchos::rcp(new Operators::PDE_AdvectionUpwind(olist_adv, mesh_));
  op_p1_sat_prec_ = Teuchos::rcp(new Operators::PDE_AdvectionUpwind(olist_adv, mesh_));
  //op_p2_sat_prec_ = Teuchos::rcp(new Operators::PDE_AdvectionUpwind(olist_adv, mesh_));

  Teuchos::ParameterList& wrm_list = plist_.sublist("Water retention models");
  //std::cout << wrm_list << "\n";
  capillary_pressure_ = Teuchos::rcp(new CapillaryPressure(mesh_));
  capillary_pressure_->Init(wrm_list); 
}

void Phase_Constraint_PK::Functional(double t_old, double t_new, 
                          Teuchos::RCP<TreeVector> u_old,
                          Teuchos::RCP<TreeVector> u_new,
                          Teuchos::RCP<TreeVector> f) 
{
  Teuchos::RCP<const CompositeVector> pressure_w = u_new->SubVector(0)->Data();
  Teuchos::RCP<const CompositeVector> saturation_n = u_new->SubVector(1)->Data();
  Teuchos::RCP<const CompositeVector> fuga_h = u_new->SubVector(3)->Data();
  Teuchos::RCP<const CompositeVector> fuga_w = u_new->SubVector(2)->Data();
  //std::cout << "pressure_w: " << *pressure_w->ViewComponent("cell") << "\n";
  //std::cout << "saturation_n: " << *saturation_n->ViewComponent("cell") << "\n";
  //std::cout << "fuga_h: " << *fuga_h->ViewComponent("cell") << "\n"; 
  //std::cout << "fuga_w: " << *fuga_w->ViewComponent("cell") << "\n"; 
  //std::cout << "Henry coef: " << henry_coef_ << "\n";

  // compute saturation non-wetting phase and pressure wetting phase  
  Teuchos::RCP<CompositeVector> saturation_n_copy = Teuchos::rcp(new CompositeVector(*saturation_n));
  Teuchos::RCP<CompositeVector> saturation_w = Teuchos::rcp(new CompositeVector(*saturation_n));
  saturation_w->Scale(-1.0);
  saturation_w->Shift(1.0);

  // new nonwet pressure
  Teuchos::RCP<CompositeVector> pressure_n = Teuchos::rcp(new CompositeVector(*pressure_w));
  capillary_pressure_->Compute(*saturation_w);  
  pressure_n->Update(1.0, *capillary_pressure_->Pc(), 1.0);
  //std::cout << "pressure_n: " << *pressure_n->ViewComponent("cell") << "\n";

  // compute molar fractions of hydrogen component
  Teuchos::RCP<CompositeVector> xhl = Teuchos::rcp(new CompositeVector(*fuga_h));
  //std::cout << "xhl: " << *xhl->ViewComponent("cell") << "\n";
  Teuchos::RCP<CompositeVector> xhg = Teuchos::rcp(new CompositeVector(*fuga_h));
  xhg->ReciprocalMultiply(1.0, *pressure_n, *fuga_h, 0.0);
  xhl->Scale(1.0/henry_coef_);
  //std::cout << "xhl: " << *xhl->ViewComponent("cell") << "\n";
  //std::cout << "xhg: " << *xhg->ViewComponent("cell") << "\n";

  // compute molar fractions of water component
  Teuchos::RCP<CompositeVector> xwl = Teuchos::rcp(new CompositeVector(*fuga_w));
  Teuchos::RCP<CompositeVector> xwg = Teuchos::rcp(new CompositeVector(*fuga_w));
  xwg->ReciprocalMultiply(1.0, *pressure_n, *fuga_w, 0.0);
  xwl->Scale(1.0/P_vap_);
  //std::cout << "xwl: " << *xwl->ViewComponent("cell") << "\n";
  //std::cout << "xwg: " << *xwg->ViewComponent("cell") << "\n";

  // create variables for active sets
  Teuchos::RCP<CompositeVector> active_p1 = Teuchos::rcp(new CompositeVector(xwl->Map()));
  Teuchos::RCP<CompositeVector> inactive_p1 = Teuchos::rcp(new CompositeVector(xwl->Map()));
  active_p1->PutScalar(0.0);
  inactive_p1->PutScalar(1.0);

  Teuchos::RCP<CompositeVector> active_p2 = Teuchos::rcp(new CompositeVector(*xwg));
  Teuchos::RCP<CompositeVector> inactive_p2 = Teuchos::rcp(new CompositeVector(*xwg));
  active_p2->PutScalar(0.0);
  inactive_p2->PutScalar(1.0);

  // complementarity functions
  // identify active set for liquid phase
  Epetra_MultiVector& xwl_c = *xwl->ViewComponent("cell");
  Epetra_MultiVector& xhl_c = *xhl->ViewComponent("cell");
  Epetra_MultiVector& active_p1_c = *active_p1->ViewComponent("cell");
  Epetra_MultiVector& inactive_p1_c = *inactive_p1->ViewComponent("cell");
  Epetra_MultiVector& sat_w_c = *saturation_w->ViewComponent("cell");
  for (int c = 0; c < xwl_c.MyLength(); c++) {
    if (sat_w_c[0][c] - (1.0 - xwl_c[0][c] - xhl_c[0][c]) > 1e-10) {
      active_p1_c[0][c] = 1.0;
      inactive_p1_c[0][c] = 0.0;
    }
  }

  // identify active set for gas phase
  Epetra_MultiVector& xwg_c = *xwg->ViewComponent("cell");
  Epetra_MultiVector& xhg_c = *xhg->ViewComponent("cell");
  Epetra_MultiVector& active_p2_c = *active_p2->ViewComponent("cell");
  Epetra_MultiVector& inactive_p2_c = *inactive_p2->ViewComponent("cell");
  Epetra_MultiVector& sat_n_c = *saturation_n_copy->ViewComponent("cell");
  for (int c = 0; c < xwg_c.MyLength(); c++) {
    if (sat_n_c[0][c] - (1.0 - xwg_c[0][c] - xhg_c[0][c]) > 1e-10) {
      active_p2_c[0][c] = 1.0;
      inactive_p2_c[0][c] = 0.0;
    }
  }

  // residual constraint equation of liquid phase
  Epetra_MultiVector& res_l_c = *f->SubVector(3)->Data()->ViewComponent("cell");
  for (int c = 0; c < res_l_c.MyLength(); c++) {
    res_l_c[0][c] = std::min(sat_w_c[0][c], 1.0 - xwl_c[0][c] - xhl_c[0][c]);
  }

  // reisudal constraint of gas phase
  Epetra_MultiVector& res_g_c = *f->SubVector(2)->Data()->ViewComponent("cell");
  for (int c = 0; c < res_g_c.MyLength(); c++) {
    res_g_c[0][c] = std::min(sat_n_c[0][c], 1.0 - xwg_c[0][c] - xhg_c[0][c]);
  }
  //std::cout << "Residual liquid phase constraint: " << *f->SubVector(3)->Data()->ViewComponent("cell") << "\n";
  //std::cout << "Residual gas phase constraint: " << *f->SubVector(2)->Data()->ViewComponent("cell") << "\n";
}

void Phase_Constraint_PK::UpdatePreconditioner(double T0, Teuchos::RCP<const TreeVector> u, double dTp) {
  Teuchos::RCP<const CompositeVector> pressure_w = u->SubVector(0)->Data();
  Teuchos::RCP<const CompositeVector> saturation_n = u->SubVector(1)->Data();
  Teuchos::RCP<const CompositeVector> fuga_h = u->SubVector(3)->Data();
  Teuchos::RCP<const CompositeVector> fuga_w = u->SubVector(2)->Data();
  //std::cout << "pressure_w: " << *pressure_w->ViewComponent("cell") << "\n";
  //std::cout << "saturation_n: " << *saturation_n->ViewComponent("cell") << "\n";
  //std::cout << "fuga_h: " << *fuga_h->ViewComponent("cell") << "\n"; 
  //std::cout << "fuga_w: " << *fuga_w->ViewComponent("cell") << "\n"; 

  // compute saturation non-wetting phase and pressure wetting phase  
  Teuchos::RCP<CompositeVector> saturation_n_copy = Teuchos::rcp(new CompositeVector(*saturation_n));
  Teuchos::RCP<CompositeVector> saturation_w = Teuchos::rcp(new CompositeVector(*saturation_n));
  saturation_w->Scale(-1.0);
  saturation_w->Shift(1.0);

  // new nonwet pressure
  Teuchos::RCP<CompositeVector> pressure_n = Teuchos::rcp(new CompositeVector(*pressure_w));
  capillary_pressure_->Compute(*saturation_w);  
  capillary_pressure_->dPc_dS()->Scale(-1.0);
  pressure_n->Update(1.0, *capillary_pressure_->Pc(), 1.0);

  // compute molar fractions of hydrogen component
  Teuchos::RCP<CompositeVector> xhl = Teuchos::rcp(new CompositeVector(*fuga_h));
  //std::cout << "xhl: " << *xhl->ViewComponent("cell") << "\n";
  Teuchos::RCP<CompositeVector> xhg = Teuchos::rcp(new CompositeVector(*fuga_h));
  xhg->ReciprocalMultiply(1.0, *pressure_n, *fuga_h, 0.0);
  xhl->Scale(1.0/henry_coef_);
  //std::cout << "xhl: " << *xhl->ViewComponent("cell") << "\n";
  //std::cout << "xhg: " << *xhg->ViewComponent("cell") << "\n";

  // compute molar fractions of water component
  Teuchos::RCP<CompositeVector> xwl = Teuchos::rcp(new CompositeVector(*fuga_w));
  Teuchos::RCP<CompositeVector> xwg = Teuchos::rcp(new CompositeVector(*fuga_w));
  xwg->ReciprocalMultiply(1.0, *pressure_n, *fuga_w, 0.0);
  xwl->Scale(1.0/P_vap_);
  //std::cout << "xwl: " << *xwl->ViewComponent("cell") << "\n";
  //std::cout << "xwg: " << *xwg->ViewComponent("cell") << "\n";

  // create variables for active sets
  Teuchos::RCP<CompositeVector> active_p1 = Teuchos::rcp(new CompositeVector(*xwl));
  Teuchos::RCP<CompositeVector> inactive_p1 = Teuchos::rcp(new CompositeVector(*xwl));
  active_p1->PutScalar(0.0);
  inactive_p1->PutScalar(1.0);

  Teuchos::RCP<CompositeVector> active_p2 = Teuchos::rcp(new CompositeVector(*xwg));
  Teuchos::RCP<CompositeVector> inactive_p2 = Teuchos::rcp(new CompositeVector(*xwg));
  active_p2->PutScalar(0.0);
  inactive_p2->PutScalar(1.0);

  // complementarity functions
  // liquid phase
  Epetra_MultiVector& xwl_c = *xwl->ViewComponent("cell");
  Epetra_MultiVector& xhl_c = *xhl->ViewComponent("cell");
  Epetra_MultiVector& active_p1_c = *active_p1->ViewComponent("cell");
  Epetra_MultiVector& inactive_p1_c = *inactive_p1->ViewComponent("cell");
  Epetra_MultiVector& sat_w_c = *saturation_w->ViewComponent("cell");
  for (int c = 0; c < xwl_c.MyLength(); c++) {
    if (sat_w_c[0][c] - (1.0 - xwl_c[0][c] - xhl_c[0][c]) > 1e-10) {
      active_p1_c[0][c] = 1.0;
      inactive_p1_c[0][c] = 0.0;
    }
  }

  // gas phase
  Epetra_MultiVector& xwg_c = *xwg->ViewComponent("cell");
  Epetra_MultiVector& xhg_c = *xhg->ViewComponent("cell");
  Epetra_MultiVector& active_p2_c = *active_p2->ViewComponent("cell");
  Epetra_MultiVector& inactive_p2_c = *inactive_p2->ViewComponent("cell");
  Epetra_MultiVector& sat_n_c = *saturation_n_copy->ViewComponent("cell");
  for (int c2 = 0; c2 < xwg_c.MyLength(); c2++) {
    if (sat_n_c[0][c2] - (1.0 - xwg_c[0][c2] - xhg_c[0][c2]) > 1e-10) {
      active_p2_c[0][c2] = 1.0;
      inactive_p2_c[0][c2] = 0.0;
    }
  }

  /*
  std::cout << "active set water phase: " << active_p1_c << "\n";
  std::cout << "inactive set water phase: " << inactive_p1_c << "\n";
  std::cout << "active set gas phase: " << active_p2_c << "\n";
  std::cout << "inactive set gas phase: " << inactive_p2_c << "\n";
  */

  Teuchos::RCP<CompositeVector> coef_p1_f1 = Teuchos::rcp(new CompositeVector(*pressure_n));
  Teuchos::RCP<CompositeVector> coef_p1_f2 = Teuchos::rcp(new CompositeVector(*pressure_n));
  Teuchos::RCP<CompositeVector> coef_p2_f1 = Teuchos::rcp(new CompositeVector(*pressure_n));
  Teuchos::RCP<CompositeVector> coef_p2_f2 = Teuchos::rcp(new CompositeVector(*pressure_n));
  Teuchos::RCP<CompositeVector> coef_p2_p1 = Teuchos::rcp(new CompositeVector(*fuga_w));
  Teuchos::RCP<CompositeVector> coef_p2_s2 = Teuchos::rcp(new CompositeVector(*fuga_w));
  coef_p1_f1->PutScalar(-1.0);
  coef_p1_f2->PutScalar(-1.0);
  coef_p2_f1->PutScalar(-1.0);
  coef_p2_f2->PutScalar(-1.0);

  coef_p1_f2->Scale(1.0/henry_coef_);
  coef_p1_f1->Scale(1.0/P_vap_);
  coef_p2_f1->ReciprocalMultiply(1.0, *pressure_n, *coef_p2_f1, 0.0);
  coef_p2_f2->ReciprocalMultiply(1.0, *pressure_n, *coef_p2_f2, 0.0);
  coef_p2_p1->Update(1.0, *fuga_h, 1.0);
  coef_p2_p1->ReciprocalMultiply(1.0, *pressure_n, *coef_p2_p1, 0.0);
  coef_p2_p1->ReciprocalMultiply(1.0, *pressure_n, *coef_p2_p1, 0.0);
  coef_p2_s2->Update(1.0, *fuga_h, 1.0);
  coef_p2_s2->ReciprocalMultiply(1.0, *pressure_n, *coef_p2_s2, 0.0);
  coef_p2_s2->ReciprocalMultiply(1.0, *pressure_n, *coef_p2_s2, 0.0);
  // Multiply here might cause assert same CVS to fail if coef_p2_s2 does not have face components
  coef_p2_s2->Multiply(1.0, *capillary_pressure_->dPc_dS(), *coef_p2_s2, 0.0);

  Epetra_MultiVector& coef_p1_f1_c = *coef_p1_f1->ViewComponent("cell");
  Epetra_MultiVector& coef_p1_f2_c = *coef_p1_f2->ViewComponent("cell");
  Epetra_MultiVector& coef_p2_f1_c = *coef_p2_f1->ViewComponent("cell");
  Epetra_MultiVector& coef_p2_f2_c = *coef_p2_f2->ViewComponent("cell");
  Epetra_MultiVector& coef_p2_p1_c = *coef_p2_p1->ViewComponent("cell");
  Epetra_MultiVector& coef_p2_s2_c = *coef_p2_s2->ViewComponent("cell");
  for (int c = 0; c < ncells_wghost_; c++) {
    double cell_volume = mesh_->cell_volume(c);
    coef_p1_f1_c[0][c] /= cell_volume;
    coef_p2_f1_c[0][c] /= cell_volume;
    coef_p1_f2_c[0][c] /= cell_volume;
    coef_p2_f2_c[0][c] /= cell_volume;
    coef_p2_p1_c[0][c] /= cell_volume;
    coef_p2_s2_c[0][c] /= cell_volume;
    inactive_p1_c[0][c] /= cell_volume;
    inactive_p2_c[0][c] /= cell_volume;
  }
  coef_p1_f1->Multiply(1.0, *coef_p1_f1, *active_p1, 0.0);
  coef_p1_f2->Multiply(1.0, *coef_p1_f2, *active_p1, 0.0);
  coef_p2_f1->Multiply(1.0, *coef_p2_f1, *active_p2, 0.0);
  coef_p2_f2->Multiply(1.0, *coef_p2_f2, *active_p2, 0.0);
  coef_p2_p1->Multiply(1.0, *coef_p2_p1, *active_p2, 0.0);
  coef_p2_s2->Multiply(1.0, *coef_p2_s2, *active_p2, 0.0);

  //std::cout << "coef_p1_f1: " << *coef_p1_f1->ViewComponent("cell") << "\n";

  // Create data structures to initialize operators
  CompositeVectorSpace cvs; 
  cvs.SetMesh(mesh_);
  cvs.SetGhosted(true);
  cvs.SetComponent("cell", AmanziMesh::CELL, 1);
  cvs.SetOwned(false);
  cvs.AddComponent("face", AmanziMesh::FACE, 1);

  Teuchos::RCP<CompositeVector> velocity = Teuchos::rcp(new CompositeVector(cvs));
  velocity->PutScalarMasterAndGhosted(0.0);

  // block for dH1/dfg_1
  op1_prec_->global_operator()->Init();
  op1_prec_->Setup(*velocity);
  op1_prec_->UpdateMatrices(*velocity);
  op1_prec_->ApplyBCs(op_bc_, true);
  op1_acc_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, op1_prec_->global_operator())); 
  op1_acc_->AddAccumulationTerm(*saturation_n, *coef_p1_f1, 1.0, "cell");

  // block for dH1/dfg_2
  op2_prec_->global_operator()->Init();
  op2_prec_->Setup(*velocity);
  op2_prec_->UpdateMatrices(*velocity);
  op2_prec_->ApplyBCs(op_bc_, true);
  op2_acc_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, op2_prec_->global_operator())); 
  op2_acc_->AddAccumulationTerm(*saturation_n, *coef_p1_f2, 1.0, "cell");

  // block for dH2/dfg_1
  op3_prec_->global_operator()->Init();
  op3_prec_->Setup(*velocity);
  op3_prec_->UpdateMatrices(*velocity);
  op3_prec_->ApplyBCs(op_bc_, true);
  op3_acc_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, op3_prec_->global_operator())); 
  op3_acc_->AddAccumulationTerm(*saturation_n, *coef_p2_f1, 1.0, "cell");

  // block for dH2/dfg_2
  op4_prec_->global_operator()->Init();
  op4_prec_->Setup(*velocity);
  op4_prec_->UpdateMatrices(*velocity);
  op4_prec_->ApplyBCs(op_bc_, true);
  op4_acc_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, op4_prec_->global_operator())); 
  op4_acc_->AddAccumulationTerm(*saturation_n, *coef_p2_f2, 1.0, "cell");

  op5_prec_->global_operator()->Init();
  op5_prec_->Setup(*velocity);
  op5_prec_->UpdateMatrices(*velocity);
  op5_prec_->ApplyBCs(op_bc_, true);
  op5_acc_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, op5_prec_->global_operator()));
  op5_acc_->AddAccumulationTerm(*saturation_n, *coef_p2_p1, 1.0, "cell");

  op6_prec_->global_operator()->Init();
  op6_prec_->Setup(*velocity);
  op6_prec_->UpdateMatrices(*velocity);
  op6_prec_->ApplyBCs(op_bc_, true);
  op6_acc_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, op6_prec_->global_operator()));
  op6_acc_->AddAccumulationTerm(*saturation_n, *coef_p2_s2, 1.0, "cell");

  op_p1_sat_prec_->global_operator()->Init();
  op_p1_sat_prec_->Setup(*velocity);
  op_p1_sat_prec_->UpdateMatrices(*velocity);
  op_p1_sat_prec_->ApplyBCs(op_bc_, true);
  op_p1_sat_acc_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, op_p1_sat_prec_->global_operator()));
  inactive_p1->Scale(-1.0);
  op_p1_sat_acc_->AddAccumulationTerm(*saturation_n, *inactive_p1, 1.0, "cell"); 

  Teuchos::ParameterList olist_adv = plist_.sublist("Constraints").sublist("operators").sublist("advection operator");
  op_p2_sat_prec_ = Teuchos::rcp(new Operators::PDE_AdvectionUpwind(olist_adv, op6_prec_->global_operator()));
  op_p2_sat_prec_->Setup(*velocity);
  op_p2_sat_prec_->UpdateMatrices(*velocity);
  op_p2_sat_prec_->ApplyBCs(op_bc_, true);
  op_p2_sat_acc_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, op_p2_sat_prec_->global_operator()));
  op_p2_sat_acc_->AddAccumulationTerm(*saturation_n, *inactive_p2, 1.0, "cell");

  /*
  //std::cout << "Assemble Saturation 1\n";
  op_p1_sat_prec_->global_operator()->SymbolicAssembleMatrix();
  op_p1_sat_prec_->global_operator()->AssembleMatrix();
  //std::cout << *op_p1_sat_prec_->global_operator()->A() << "\n";

  //std::cout << "Assemble fugacity liquid block 1 \n";
  op1_prec_->global_operator()->SymbolicAssembleMatrix();
  op1_prec_->global_operator()->AssembleMatrix();
  //std::cout << *op1_prec_->global_operator()->A() << "\n";

  //std::cout << "Assemble fugacity liquid block 2 \n";
  op2_prec_->global_operator()->SymbolicAssembleMatrix();
  op2_prec_->global_operator()->AssembleMatrix();
  //std::cout << *op2_prec_->global_operator()->A() << "\n";

  //std::cout << "Assemble Saturation 2\n";
  op_p2_sat_prec_->global_operator()->SymbolicAssembleMatrix();
  op_p2_sat_prec_->global_operator()->AssembleMatrix();
  //std::cout << *op_p2_sat_prec_->global_operator()->A() << "\n";

  //std::cout << "Assemble Matrix 1\n";
  op1_prec_->global_operator()->SymbolicAssembleMatrix();
  op1_prec_->global_operator()->AssembleMatrix();
  //std::cout << *op1_prec_->global_operator()->A() << "\n";

  //std::cout << "Assemble Matrix 2\n";
  op2_prec_->global_operator()->SymbolicAssembleMatrix();
  op2_prec_->global_operator()->AssembleMatrix();
  //std::cout << *op2_prec_->global_operator()->A() << "\n";

  //std::cout << "Assemble Matrix 3\n";
  op3_prec_->global_operator()->SymbolicAssembleMatrix();
  op3_prec_->global_operator()->AssembleMatrix();
  //std::cout << *op3_prec_->global_operator()->A() << "\n";

  //std::cout << "Assemble Matrix 4\n";
  op4_prec_->global_operator()->SymbolicAssembleMatrix();
  op4_prec_->global_operator()->AssembleMatrix();
  //std::cout << *op4_prec_->global_operator()->A() << "\n";

  //std::cout << "Assemble Matrix 5\n";
  op5_prec_->global_operator()->SymbolicAssembleMatrix();
  op5_prec_->global_operator()->AssembleMatrix();
  //std::cout << *op5_prec_->global_operator()->A() << "\n";

  //std::cout << "Assemble Matrix 6\n";
  op6_prec_->global_operator()->SymbolicAssembleMatrix();
  op6_prec_->global_operator()->AssembleMatrix();
  //std::cout << *op6_prec_->global_operator()->A() << "\n";
  */
}

}
}
