#include "GasConstraint.hh"


namespace Amanzi {
namespace Multiphase {

GasConstraint::GasConstraint(Teuchos::ParameterList plist,
                                          const Teuchos::RCP<State> S):
                  S_(S), plist_(plist), passwd_("state") 
{
  mesh_ = S->GetMesh();
  //henry_coef_ = plist_.sublist("EOS").sublist("Component 2").get<double>("coefficient value");
  mu_ = 0.0;
  ncp_type_ = plist_.get<std::string>("NCP function", "min");
  H_ = 7.65e-6;
  M_h_ = 2e-3;
}

void GasConstraint::Initialize()
{
  ncells_owned_ = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  ncells_wghost_ = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);

  nfaces_owned_ = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  nfaces_wghost_ = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  // Allocate memory for boundary data.
  bc_model_.resize(nfaces_wghost_, 0);
  bc_value_.resize(nfaces_wghost_, 0.0);
  bc_mixed_.resize(nfaces_wghost_, 0.0);

  op_bc_ = Teuchos::rcp(new Operators::BCs(mesh_, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));

  // preconditioners
  op1_acc_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, mesh_)); 
  op2_acc_= Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, mesh_));
  op3_acc_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, mesh_));
  op1_acc_tmp_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, mesh_)); 
  op2_acc_tmp_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, mesh_));
  op3_acc_tmp_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, mesh_));

  Teuchos::ParameterList& wrm_list = plist_.sublist("Water retention models");
  capillary_pressure_ = Teuchos::rcp(new CapillaryPressureOld(mesh_));
  capillary_pressure_->Init(wrm_list); 
}

void GasConstraint::FunctionalResidual(double t_old, double t_new, 
                                       Teuchos::RCP<TreeVector> u_old,
                                       Teuchos::RCP<TreeVector> u_new,
                                       Teuchos::RCP<TreeVector> f) 
{
  Teuchos::RCP<const CompositeVector> pressure_w = u_new->SubVector(0)->Data();
  Teuchos::RCP<const CompositeVector> saturation_w = u_new->SubVector(1)->Data();
  Teuchos::RCP<const CompositeVector> rhl = u_new->SubVector(2)->Data();
  double Cg = H_ * M_h_;

  // new nonwet pressure
  Teuchos::RCP<CompositeVector> pressure_n = Teuchos::rcp(new CompositeVector(*pressure_w));
  capillary_pressure_->Compute(*saturation_w);  
  pressure_n->Update(1.0, *capillary_pressure_->Pc(), 1.0);

  // create variables for active sets
  Teuchos::RCP<CompositeVector> active_p_g = Teuchos::rcp(new CompositeVector(rhl->Map()));
  Teuchos::RCP<CompositeVector> inactive_p_g = Teuchos::rcp(new CompositeVector(rhl->Map()));
  active_p_g->PutScalar(0.0);
  inactive_p_g->PutScalar(1.0);

  // identify active set for gas phase
  const Epetra_MultiVector& rhl_c = *rhl->ViewComponent("cell");
  Epetra_MultiVector& active_p_g_c = *active_p_g->ViewComponent("cell");
  Epetra_MultiVector& inactive_p_g_c = *inactive_p_g->ViewComponent("cell");
  const Epetra_MultiVector& p2_c = *pressure_n->ViewComponent("cell");
  const Epetra_MultiVector& sat_w_c = *saturation_w->ViewComponent("cell");
  for (int c = 0; c < rhl_c.MyLength(); c++) {
    if (1.0 - sat_w_c[0][c] - (Cg * p2_c[0][c] - rhl_c[0][c]) > 1e-12) {
      active_p_g_c[0][c] = 1.0;
      inactive_p_g_c[0][c] = 0.0;
    }
  }

  // reisudal constraint of gas phase
  Epetra_MultiVector& res_g_c = *f->Data()->ViewComponent("cell");
  if (ncp_type_ == "min") {
    for (int c = 0; c < res_g_c.MyLength(); c++) {
      res_g_c[0][c] = std::min(1.0 - sat_w_c[0][c], Cg * p2_c[0][c] - rhl_c[0][c]);
    }
  } else if (ncp_type_ == "fischer-burmeister") {
    for (int c = 0; c < res_g_c.MyLength(); c++) {
      res_g_c[0][c] = pow(pow(1.0 - sat_w_c[0][c], 2.0) + pow(Cg * p2_c[0][c] - rhl_c[0][c], 2.0), 0.5) - 
        (1.0 - sat_w_c[0][c]) - (Cg * p2_c[0][c] - rhl_c[0][c]);
    }
  } else {
    Errors::Message msg;
    msg << "Reduced2p2c_PK: unknown NCP function";
    Exceptions::amanzi_throw(msg);
  }
}

void GasConstraint::UpdatePreconditioner(double T0, Teuchos::RCP<const TreeVector> u, double dTp)
{
  Teuchos::RCP<const CompositeVector> pressure_w = u->SubVector(0)->Data();
  Teuchos::RCP<const CompositeVector> saturation_w = u->SubVector(1)->Data();
  Teuchos::RCP<const CompositeVector> rhl = u->SubVector(2)->Data();
  double Cg = H_ * M_h_;

  // new nonwet pressure
  Teuchos::RCP<CompositeVector> pressure_n = Teuchos::rcp(new CompositeVector(*pressure_w));
  capillary_pressure_->Compute(*saturation_w);  
  pressure_n->Update(1.0, *capillary_pressure_->Pc(), 1.0);

  // create variables for active sets
  Teuchos::RCP<CompositeVector> active_p_g = Teuchos::rcp(new CompositeVector(rhl->Map()));
  Teuchos::RCP<CompositeVector> inactive_p_g = Teuchos::rcp(new CompositeVector(rhl->Map()));
  active_p_g->PutScalar(0.0);
  inactive_p_g->PutScalar(1.0);

  // identify active set for gas phase
  const Epetra_MultiVector& rhl_c = *rhl->ViewComponent("cell");
  Epetra_MultiVector& active_p_g_c = *active_p_g->ViewComponent("cell");
  Epetra_MultiVector& inactive_p_g_c = *inactive_p_g->ViewComponent("cell");
  const Epetra_MultiVector& p2_c = *pressure_n->ViewComponent("cell");
  const Epetra_MultiVector& sat_w_c = *saturation_w->ViewComponent("cell");
  const Epetra_MultiVector& dPc_dS_c = *capillary_pressure_->dPc_dS()->ViewComponent("cell");
  cnt_ = 0; // number of cells without gas
  int ncells = rhl_c.MyLength();
  for (int c = 0; c < ncells; c++) {
    if (std::abs(1.0 - sat_w_c[0][c]) - std::abs(Cg * p2_c[0][c] - rhl_c[0][c]) > 1e-12) {
      active_p_g_c[0][c] = 1.0;
      inactive_p_g_c[0][c] = 0.0;
    } else {
      cnt_++;
    }
  }

  Teuchos::RCP<CompositeVector> coef_df_dp = Teuchos::rcp(new CompositeVector(*pressure_n));
  Teuchos::RCP<CompositeVector> coef_df_ds = Teuchos::rcp(new CompositeVector(*pressure_n));
  Teuchos::RCP<CompositeVector> coef_df_dr = Teuchos::rcp(new CompositeVector(*pressure_n));
  Epetra_MultiVector& coef_df_dp_c = *coef_df_dp->ViewComponent("cell");
  Epetra_MultiVector& coef_df_ds_c = *coef_df_ds->ViewComponent("cell");
  Epetra_MultiVector& coef_df_dr_c = *coef_df_dr->ViewComponent("cell");

  /* vectors for the experimental smoothed preconditioner */
  Teuchos::RCP<CompositeVector> coef_df_dp_tmp = Teuchos::rcp(new CompositeVector(*pressure_n));
  Teuchos::RCP<CompositeVector> coef_df_ds_tmp = Teuchos::rcp(new CompositeVector(*pressure_n));
  Teuchos::RCP<CompositeVector> coef_df_dr_tmp = Teuchos::rcp(new CompositeVector(*pressure_n));
  Epetra_MultiVector& coef_df_dp_tmp_c = *coef_df_dp_tmp->ViewComponent("cell");
  Epetra_MultiVector& coef_df_ds_tmp_c = *coef_df_ds_tmp->ViewComponent("cell");
  Epetra_MultiVector& coef_df_dr_tmp_c = *coef_df_dr_tmp->ViewComponent("cell");

  if (ncp_type_ == "min") {
    for (int c = 0; c < ncells; c++) {
      double cell_volume = mesh_->cell_volume(c);
      int is_inactive = round(inactive_p_g_c[0][c]);
      if (is_inactive) {
        coef_df_dp_c[0][c] = 0.0;
        coef_df_ds_c[0][c] =  - 1.0 / cell_volume;
        coef_df_dr_c[0][c] = 0.0;
      } else {
        coef_df_dp_c[0][c] = Cg / cell_volume;
        coef_df_ds_c[0][c] =  Cg * dPc_dS_c[0][c] / cell_volume;
        coef_df_dr_c[0][c] = - 1.0 / cell_volume;
      }
    }
  } else if (ncp_type_ == "fischer-burmeister") {
    // the nonlinear ncp function phi(a,b) = sqrt(a^2 + b^2) - a - b
    // with a = S_n = 1.0 - S_w and b = Cg * P_n - rhl
    //if (std::abs(mu_) < 1e-15) {
      for (int c = 0; c < ncells; c++) {
        double cell_volume = mesh_->cell_volume(c);
        double func_a = 1.0 - sat_w_c[0][c];
        double func_b = Cg * p2_c[0][c] - rhl_c[0][c];
        int is_inactive = round(inactive_p_g_c[0][c]);
        if (is_inactive) {
          coef_df_dp_c[0][c] = Cg * (func_b * pow(pow(func_b, 2.0) + 2.0 * mu_, -0.5) - 1.0) / cell_volume;
          coef_df_ds_c[0][c] = (Cg * dPc_dS_c[0][c] * (func_b * pow(pow(func_b, 2.0) + 2.0 * mu_, -0.5) - 1.0) + 1.0) / cell_volume;
          coef_df_dr_c[0][c] = - (func_b * pow(pow(func_b, 2.0) + 2.0 * mu_, -0.5) - 1.0) / cell_volume;
        } else {
          coef_df_dp_c[0][c] = - Cg / cell_volume;
          coef_df_ds_c[0][c] = - ((func_a * pow(pow(func_a, 2.0) + 2.0 * mu_, -0.5) - 1.0) + Cg * dPc_dS_c[0][c]) / cell_volume;
          coef_df_dr_c[0][c] = 1.0 / cell_volume;
        }
      }
    //} else {
    //  for (int c = 0; c < ncells; c++) {
    //    double cell_volume = mesh_->cell_volume(c);
    //    double func_a = 1.0 - sat_w_c[0][c];
    //    double func_b = Cg * p2_c[0][c] - rhl_c[0][c];
    //    double factor = pow(pow(func_b, 2.0) + pow(func_a, 2.0) + 2.0 * mu_, -0.5);
    //    coef_df_dp_c[0][c] = Cg * (func_b * factor - 1.0) / cell_volume;
    //    coef_df_ds_c[0][c] = (Cg * dPc_dS_c[0][c] * (func_b * factor - 1.0) - func_a * factor + 1.0) / cell_volume;
    //    coef_df_dr_c[0][c] = - (func_b * factor - 1.0) / cell_volume;
    //  }
    //}
  } else {
    Errors::Message msg;
    msg << "Reduced2p2c_PK: unknown NCP function";
    Exceptions::amanzi_throw(msg);
  }

  // count the number of zero diagonal values for FB function
  if (ncp_type_ == "fischer-burmeister") {
    cnt_ = 0;
    for (int c = 0; c < ncells; c++) {
      if (std::abs(coef_df_dr_c[0][c]) < 1e-15) {
        cnt_++;
      }
    }
  }

  // identify the indices for zero diagonal values
  int cnt_idx = 0;
  inactive_gas_idx_ = new int[cnt_];
  if (ncp_type_ == "min") {
    for (int c = 0; c < ncells; c++) {
      int is_inactive = round(inactive_p_g_c[0][c]);
      if (is_inactive) {
        inactive_gas_idx_[cnt_idx] = 3*(c + 1);
        cnt_idx++;
      }
    }
  }
  else if (ncp_type_ == "fischer-burmeister") {
    for (int c = 0; c < ncells; c++) {
      if (std::abs(coef_df_dr_c[0][c]) < 1e-15) {
        inactive_gas_idx_[cnt_idx] = 3*(c + 1);
        cnt_idx++;
      }
    }
  }
  AMANZI_ASSERT(cnt_idx == cnt_);

  // block for dH1/dfg_1
  op1_acc_->global_operator()->Init();
  op1_acc_->AddAccumulationDelta(*saturation_w, *coef_df_dp, *coef_df_dp, 1.0, "cell");

  op2_acc_->global_operator()->Init();
  op2_acc_->AddAccumulationDelta(*saturation_w, *coef_df_ds, *coef_df_ds, 1.0, "cell");

  // block for dH2/dfg_1
  op3_acc_->global_operator()->Init();
  op3_acc_->AddAccumulationDelta(*saturation_w, *coef_df_dr, *coef_df_dr, 1.0, "cell");

  /* experiment with smoothed preconditioner with parameter mu */
  // block for dH1/dfg_1
  op1_acc_tmp_->global_operator()->Init();
  op1_acc_tmp_->AddAccumulationDelta(*saturation_w, *coef_df_dp_tmp, *coef_df_dp_tmp, 1.0, "cell");

  op2_acc_tmp_->global_operator()->Init();
  op2_acc_tmp_->AddAccumulationDelta(*saturation_w, *coef_df_ds_tmp, *coef_df_ds_tmp, 1.0, "cell");

  // block for dH2/dfg_1
  op3_acc_tmp_->global_operator()->Init();
  op3_acc_tmp_->AddAccumulationDelta(*saturation_w, *coef_df_dr_tmp, *coef_df_dr_tmp, 1.0, "cell");
}

}
}
