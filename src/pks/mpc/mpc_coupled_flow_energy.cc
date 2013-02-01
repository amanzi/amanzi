/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Derived MPC for flow and energy.  This couples using a block-diagonal coupler.
------------------------------------------------------------------------- */

#include "field_evaluator.hh"
#include "mpc_coupled_flow_energy.hh"
#include "strong_mpc.hh"
#include "Epetra_FECrsGraph.h"

#include "errors.hh"
#include "exceptions.hh"

#include "Trilinos_Util_CrsMatrixGallery.h"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"

namespace Amanzi {

#define DEBUG_FLAG 1

RegisteredPKFactory<MPCCoupledFlowEnergy> MPCCoupledFlowEnergy::reg_("energy-flow preconditioner coupled");


void MPCCoupledFlowEnergy::initialize(const Teuchos::Ptr<State>& S) {
  StrongMPC::initialize(S);

  is_matrix_constructed = false;
  decoupled = false;

  coupled_pc_ = plist_.sublist("Coupled PC");

  std::string prec_name = coupled_pc_.get<std::string>("preconditioner", "ILU");
  if (prec_name == "ML") {
    prec_method_ = TRILINOS_ML;
  } else if (prec_name == "ILU") {
    prec_method_ = TRILINOS_ILU;
  } else if (prec_name == "Block ILU") {
    prec_method_ = TRILINOS_BLOCK_ILU;
#ifdef HAVE_HYPRE
  } else if (prec_name == "HYPRE AMG") {
    prec_method_ = HYPRE_AMG;
  } else if (prec_name == "HYPRE Euclid") {
    prec_method_ = HYPRE_EUCLID;
  } else if (prec_name == "HYPRE ParaSails") {
    prec_method_ = HYPRE_EUCLID;
#endif
  } else {
#ifdef HAVE_HYPRE
    Errors::Message msg("Matrix_MFD: The specified preconditioner "+prec_name+" is not supported, we only support ML, ILU, HYPRE AMG, HYPRE Euclid, and HYPRE ParaSails");
#else
    Errors::Message msg("Matrix_MFD: The specified preconditioner "+prec_name+" is not supported, we only support ML, and ILU");
#endif
    Exceptions::amanzi_throw(msg);
  }


  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData("pressure");
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData("temperature");

  ASSERT(temp->size("cell") == pres->size("cell"));

  D_pT_ = Teuchos::rcp(new Epetra_MultiVector(*pres->ViewComponent("cell", false)));
  D_Tp_ = Teuchos::rcp(new Epetra_MultiVector(*temp->ViewComponent("cell", false)));

  flow_pk = Teuchos::rcp_dynamic_cast<Amanzi::Flow::Richards>(sub_pks_[0]);
  energy_pk = Teuchos::rcp_dynamic_cast<Amanzi::Energy::TwoPhase>(sub_pks_[1]);

  SymbolicAssembleGlobalMatrices_(S);
};

void MPCCoupledFlowEnergy::SymbolicAssembleGlobalMatrices_(const Teuchos::Ptr<State>& S) {
  int ierr(0);
  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData("pressure");

  Teuchos::RCP<const Epetra_BlockMap> cmap = pres->map("cell", false);
  Teuchos::RCP<const Epetra_BlockMap> fmap = pres->map("face", false);
  fmap_wghost = pres->map("face", true);

  int num_global = (*fmap).NumGlobalPoints(); //std::cout<<"num_global "<<num_global<<std::endl;
  int nummypoint = (*fmap).NumMyPoints(); //std::cout<<"num_mine "<<nummypoint<<std::endl;

  double_fmap = Teuchos::rcp(new Epetra_BlockMap((*fmap).NumGlobalPoints(),
          (*fmap).NumMyPoints(), (*fmap).MyGlobalElements(), 2,
          (*fmap).IndexBase(), (*fmap).Comm() ));

  double_cmap = Teuchos::rcp(new Epetra_BlockMap((*cmap).NumGlobalPoints(),
          (*cmap).NumMyPoints(), (*cmap).MyGlobalElements(), 2,
          (*cmap).IndexBase(), (*cmap).Comm() ));


  double_fmap_wghost = Teuchos::rcp(new Epetra_BlockMap((*fmap_wghost).NumGlobalPoints(),
          (*fmap_wghost).NumMyPoints(), (*fmap_wghost).MyGlobalElements(), 2,
          (*fmap_wghost).IndexBase(), (*fmap_wghost).Comm() ));

//   int avg_entries_row = (mesh_->space_dimension() == 2) ? MFD_QUAD_FACES : MFD_HEX_FACES;
  int avg_entries_row = 6;
  Epetra_CrsGraph cf_graph(Copy, *double_cmap, *double_fmap_wghost, avg_entries_row, false);
  Epetra_FECrsGraph ff_graph(Copy, *double_fmap, 2*avg_entries_row - 1, false);
  //  Epetra_FECrsGraph ff_graph(Copy, *double_fmap_wghost, 2*avg_entries_row - 1, false);

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;
  const int MFD_MAX_FACES = 14;
  int faces_LID[MFD_MAX_FACES];  // Contigious memory is required.
  int faces_GID[MFD_MAX_FACES];

  int ncells = flow_pk->mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c!=ncells; ++c) {
    flow_pk->mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    for (int n=0; n!=nfaces; ++n) {
      faces_LID[n] = faces[n];
      faces_GID[n] = double_fmap_wghost->GID(faces_LID[n]);
    }
    cf_graph.InsertMyIndices(c, nfaces, faces_LID);
    ierr = ff_graph.InsertGlobalIndices(nfaces, faces_GID, nfaces, faces_GID);
    ASSERT(!ierr);
  }
  cf_graph.FillComplete(*double_fmap, *double_cmap);
  ierr = ff_graph.GlobalAssemble();  // Symbolic graph is complete.
  ASSERT(!ierr);

  A2f2p_ = Teuchos::rcp(new Epetra_VbrMatrix(Copy, cf_graph));
  P2f2f_ = Teuchos::rcp(new Epetra_FEVbrMatrix(Copy, ff_graph, false));
  ierr = P2f2f_->GlobalAssemble();
  ASSERT(!ierr);

  InitPreconditioner_(coupled_pc_);
}

void MPCCoupledFlowEnergy::update_precon(double t, Teuchos::RCP<const TreeVector> up,
        double h) {
  StrongMPC::update_precon(t, up, h);
  int ierr(0);

  // d pressure_residual / d T
  // -- update the accumulation derivative
  S_next_->GetFieldEvaluator("water_content")
      ->HasFieldDerivativeChanged(S_next_.ptr(), "richards_pk", "temperature");

  // -- get the accumulation deriv
  const Epetra_MultiVector& dwc_dT =
      *S_next_->GetFieldData("dwater_content_dtemperature")->ViewComponent("cell",false);
  int ncells = dwc_dT.MyLength();
  for (int c=0; c!=ncells; ++c) {
    if (!decoupled) {
      (*D_pT_)[0][c] = dwc_dT[0][c] / h;
    } else {
      (*D_pT_)[0][c] = 0.;
    }
  }

  // d temperature_residual / d p
  // -- update the accumulation derivative
  S_next_->GetFieldEvaluator("energy")
      ->HasFieldDerivativeChanged(S_next_.ptr(), "energy_pk", "pressure");

  // -- get the accumulation deriv
  const Epetra_MultiVector& de_dp =
      *S_next_->GetFieldData("denergy_dpressure")->ViewComponent("cell",false);
  for (int c=0; c!=ncells; ++c) {
    if (!decoupled)
      (*D_Tp_)[0][c] = de_dp[0][c] / h;
    else
      (*D_Tp_)[0][c] = 0.;
  }

  ComputeShurComplementPK_();

  if (prec_method_ == TRILINOS_ML) {
    if (ml_prec_->IsPreconditionerComputed()) ml_prec_->DestroyPreconditioner();
    ml_prec_->SetParameterList(ml_plist_);
    ml_prec_->ComputePreconditioner();

  } else if (prec_method_ == TRILINOS_ILU) {
    ilu_prec_ = Teuchos::rcp(new Ifpack_ILU(&*P2f2f_));
    ilu_prec_->SetParameters(ilu_plist_);
    ilu_prec_->Initialize();
    ilu_prec_->Compute();

  } else if (prec_method_ == TRILINOS_BLOCK_ILU) {
    Ifpack factory;
    std::string prectype("ILU");
    int ovl = ifp_plist_.get<int>("overlap",0);
    ifp_plist_.set<std::string>("schwarz: combine mode","Add");
    ifp_prec_ = Teuchos::rcp(factory.Create(prectype, &*P2f2f_, ovl));
    ifp_prec_->SetParameters(ifp_plist_);
    ifp_prec_->Initialize();
    ifp_prec_->Compute();
#ifdef HAVE_HYPRE
  } else if (prec_method_ == HYPRE_AMG) {
    IfpHypre_Sff_ = Teuchos::rcp(new Ifpack_Hypre(&*P2f2f_));
    Teuchos::RCP<FunctionParameter> functs[8];
    functs[0] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetCoarsenType, 0));
    functs[1] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetPrintLevel, 0)); 
    functs[2] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetNumSweeps, hypre_nsmooth_));
    functs[3] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetMaxIter, hypre_ncycles_));
    functs[4] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetRelaxType, 6)); 
    functs[5] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetStrongThreshold, hypre_strong_threshold_)); 
    functs[6] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetTol, hypre_tol_)); 
    functs[7] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetCycleType, 1));  
    
    Teuchos::ParameterList hypre_list;
    hypre_list.set("Preconditioner", BoomerAMG);
    hypre_list.set("SolveOrPrecondition", Preconditioner);
    hypre_list.set("SetPreconditioner", true);
    hypre_list.set("NumFunctions", 8);
    hypre_list.set<Teuchos::RCP<FunctionParameter>*>("Functions", functs); 
    
    IfpHypre_Sff_->SetParameters(hypre_list);
    IfpHypre_Sff_->Initialize();
    IfpHypre_Sff_->Compute();
  } else if (prec_method_ == HYPRE_EUCLID) {
    IfpHypre_Sff_ = Teuchos::rcp(new Ifpack_Hypre(&*P2f2f_));

    Teuchos::ParameterList hypre_list;
    hypre_list.set("Preconditioner", Euclid);
    hypre_list.set("SolveOrPrecondition", Preconditioner);
    hypre_list.set("SetPreconditioner", true);
    hypre_list.set("NumFunctions", 0);
    
    IfpHypre_Sff_->SetParameters(hypre_list);
    IfpHypre_Sff_->Initialize();
    IfpHypre_Sff_->Compute();    
  } else if (prec_method_ == HYPRE_PARASAILS) {
    IfpHypre_Sff_ = Teuchos::rcp(new Ifpack_Hypre(&*P2f2f_));

    Teuchos::ParameterList hypre_list;
    hypre_list.set("Preconditioner", ParaSails);
    hypre_list.set("SolveOrPrecondition", Preconditioner);
    hypre_list.set("SetPreconditioner", true);
    hypre_list.set("NumFunctions", 0);
    
    IfpHypre_Sff_->SetParameters(hypre_list);
    IfpHypre_Sff_->Initialize();
    IfpHypre_Sff_->Compute();    
#endif
  }
};

void MPCCoupledFlowEnergy::ComputeShurComplementPK_(){
  int ierr(0);

  Teuchos::RCP<const CompositeVector> pres = S_->GetFieldData("pressure");
  Teuchos::RCP<const Epetra_BlockMap> cmap = pres->map("cell", false);
  Teuchos::RCP<const Epetra_BlockMap> fmap = pres->map("face", false);

  int ncells = pres->size("cell",false);
  int cell_GID ; 

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;
  const int MFD_MAX_FACES = 14;
  int faces_LID[MFD_MAX_FACES];  // Contigious memory is required.
  int faces_GID[MFD_MAX_FACES];

  std::vector<Teuchos::SerialDenseMatrix<int, double> >& flow_Aff =
      flow_pk->preconditioner_->Aff_cells();
  std::vector<double>& flow_Acc = flow_pk->preconditioner_->Acc_cells();
  std::vector<Epetra_SerialDenseVector>& flow_Afc = flow_pk->preconditioner_->Afc_cells();
  std::vector<Epetra_SerialDenseVector>& flow_Acf = flow_pk->preconditioner_->Acf_cells();

  std::vector<Teuchos::SerialDenseMatrix<int, double> >& energy_Aff =
      energy_pk->preconditioner_->Aff_cells();
  std::vector<double>& energy_Acc = energy_pk->preconditioner_->Acc_cells();
  std::vector<Epetra_SerialDenseVector>& energy_Afc = energy_pk->preconditioner_->Afc_cells();
  std::vector<Epetra_SerialDenseVector>& energy_Acf = energy_pk->preconditioner_->Acf_cells();

  Teuchos::SerialDenseMatrix<int, double> Couple_Inv(2, 2);
  Epetra_SerialDenseMatrix Values(2, 2);
  int nfaces = 6;

  Epetra_SerialDenseMatrix Schur(2*nfaces, 2*nfaces);
  Epetra_SerialDenseMatrix Apf(2, 2*nfaces);

  Cell_Couple_Inv_.clear();

  if (is_matrix_constructed) P2f2f_->PutScalar(0.0);

  for (int c=0; c < ncells; ++c){
    cell_GID = cmap->GID(c);
    flow_pk->mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);

    double det_cell_couple = flow_Acc[c] * energy_Acc[c] - (*D_pT_)[0][c] * (*D_Tp_)[0][c];

    if (det_cell_couple != 0.) {
      Couple_Inv(0, 0) = energy_Acc[c]/det_cell_couple;
      Couple_Inv(1, 1) = flow_Acc[c]/det_cell_couple;
      Couple_Inv(0, 1) = -(*D_pT_)[0][c]/det_cell_couple;
      Couple_Inv(1, 0) = -(*D_Tp_)[0][c]/det_cell_couple;
    } else {
      Errors::Message m("Division by zero: determinant of Cell_Couple is zero");
      Exceptions::amanzi_throw(m);
    }

    for (int i=0; i!=nfaces; ++i) {
      for (int j=0; j!=nfaces; ++j) {
        Schur(i, j) = flow_Aff[c](i, j) - flow_Afc[c](i)*Couple_Inv(0, 0)*flow_Acf[c](j);
      }
    }

    for (int i=0; i!=nfaces; ++i) {
      for (int j=0; j!=nfaces; ++j) {
        Schur(nfaces + i, nfaces + j) = energy_Aff[c](i, j) - energy_Afc[c](i)*Couple_Inv(1, 1)*energy_Acf[c](j);
      }
    }

    for (int i=0; i!=nfaces; ++i) {
      for (int j=0; j!=nfaces; ++j) {
        Schur(i, nfaces + j) = - flow_Afc[c](i)*Couple_Inv(0, 1)*energy_Acf[c](j);
      }
    }

    for (int i=0; i!=nfaces; ++i) {
      for (int j=0; j!=nfaces; ++j) {
        Schur(nfaces + i, j) = - energy_Afc[c](i)*Couple_Inv(1, 0)*flow_Acf[c](j);
      }
    }

    for (int i=0; i!=nfaces; ++i) {
      Apf(0,i) =           flow_Acf[c](i);
      Apf(0,i + nfaces)  = 0;
      Apf(1,i) =           0;
      Apf(1,i + nfaces)  = energy_Afc[c](i);
    }

    Cell_Couple_Inv_.push_back(Couple_Inv); 

    int NumEntries = nfaces;
    for (int i=0; i!=nfaces; ++i) {
      faces_LID[i] = faces[i];
      faces_GID[i] = fmap_wghost->GID(faces_LID[i]);
    }

    for (int i=0; i!=nfaces; ++i) {
      ierr = P2f2f_->BeginSumIntoGlobalValues(faces_GID[i], NumEntries, faces_GID);
      ASSERT(!ierr);

      for (int j=0; j!=nfaces; ++j){
        Values(0,0) = Schur(i,j);
        Values(0,1) = Schur(i,j + nfaces);
        Values(1,0) = Schur(i + nfaces,j);
        Values(1,1) = Schur(i+ nfaces,j+ nfaces);

        //ierr = P2f2f_->SubmitBlockEntry(Values);
        ierr = P2f2f_->SubmitBlockEntry(Values.A(), Values.LDA(), Values.M(), Values.N());
        ASSERT(!ierr);
      }

      ierr = P2f2f_->EndSubmitEntries();
      ASSERT(!ierr);
    }

    ierr = A2f2p_->BeginReplaceGlobalValues(cell_GID, NumEntries, faces_GID);
    ASSERT(!ierr);
    for (int i=0; i!=nfaces; ++i) {
      Values(0,0) = Apf(0,i);
      Values(0,1) = Apf(0,i + nfaces);
      Values(1,0) = Apf(1,i);
      Values(1,1) = Apf(1,i + nfaces);
      ierr = A2f2p_->SubmitBlockEntry(Values);
      ASSERT(!ierr);
    }
    ierr = A2f2p_->EndSubmitEntries();
    ASSERT(!ierr);

  }

  ierr = A2f2p_->FillComplete(*double_fmap, *double_cmap);
  ASSERT(!ierr);

  ierr = P2f2f_->GlobalAssemble();
  ASSERT(!ierr);


  is_matrix_constructed = true;
}


// -----------------------------------------------------------------------------
// Apply the preconditioner to u and return the result in Pu.
// -----------------------------------------------------------------------------
void MPCCoupledFlowEnergy::precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) {

  Teuchos::RCP<const TreeVector> pressure = u->SubVector("flow");
  Teuchos::RCP<const CompositeVector> pressure_data = pressure->data();

  Teuchos::RCP<const TreeVector> temp = u->SubVector("energy");
  Teuchos::RCP<const CompositeVector> temp_data = temp->data();

  Epetra_MultiVector pres_c = *(pressure_data->ViewComponent("cell", false));
  Epetra_MultiVector temp_c = *(temp_data->ViewComponent("cell", false));

  Epetra_MultiVector pres_f = *pressure_data->ViewComponent("face", false);
  Epetra_MultiVector temp_f = *temp_data->ViewComponent("face", false);

  // Temporary cell and face vectors.
  Epetra_MultiVector Tc(*double_cmap, 1);
  Epetra_MultiVector Tf(*double_fmap, 1);

  Epetra_MultiVector Pc(*double_cmap, 1);
  Epetra_MultiVector Pf(*double_fmap, 1);

  Epetra_MultiVector test_res(*double_fmap, 1);

  int ierr(0);

  int ncells = pres_c.MyLength();
  int nfaces = pres_f.MyLength();
  double val = 0.;
 for (int c=0; c!=ncells; ++c){
   double p_c = pres_c[0][c];
   double t_c = temp_c[0][c];
   val = -(Cell_Couple_Inv_[c](0,0)*p_c + Cell_Couple_Inv_[c](0,1)*t_c);
   ierr = Tc.ReplaceMyValue(c, 0, 0, val);
   ASSERT(!ierr);
   val = -(Cell_Couple_Inv_[c](1,0)*p_c + Cell_Couple_Inv_[c](1,1)*t_c);
   ierr = Tc.ReplaceMyValue(c, 1, 0, val);
   ASSERT(!ierr);
 }

 ierr = A2f2p_->Multiply(true, Tc, Tf); // It performs the required parallel communications.
 ASSERT(!ierr);

 for (int f=0; f!=nfaces; ++f){
   ierr = Tf.SumIntoMyValue(f, 0, 0, pres_f[0][f]);
   ASSERT(!ierr);
   ierr = Tf.SumIntoMyValue(f, 1, 0, temp_f[0][f]);
   ASSERT(!ierr);
 }

 if (prec_method_ == TRILINOS_ML) {
   ierr = ml_prec_->ApplyInverse(Tf, Pf);
 } else if (prec_method_ == TRILINOS_ILU) {
   ierr = ilu_prec_->ApplyInverse(Tf, Pf);
 } else if (prec_method_ == TRILINOS_BLOCK_ILU) {
   ierr = ifp_prec_->ApplyInverse(Tf, Pf);
#ifdef HAVE_HYPRE
 } else if (prec_method_ == HYPRE_AMG || prec_method_ == HYPRE_EUCLID || prec_method_ == HYPRE_PARASAILS) {
   ierr = IfpHypre_Sff_->ApplyInverse(Tf, Pf);
#endif
 } else {
   ASSERT(0);
 }
 ASSERT(!ierr);


 Epetra_MultiVector& Ppressure_c = *Pu->SubVector("flow")->data()->ViewComponent("cell",false);
 Epetra_MultiVector& Ppressure_f = *Pu->SubVector("flow")->data()->ViewComponent("face",false);
 Epetra_MultiVector& Ptemp_c = *Pu->SubVector("energy")->data()->ViewComponent("cell",false);
 Epetra_MultiVector& Ptemp_f = *Pu->SubVector("energy")->data()->ViewComponent("face",false);
 Teuchos::RCP<TreeVector> Ppressure = Pu->SubVector("flow");
 Teuchos::RCP<CompositeVector> Ppressure_data = Ppressure->data();
 Teuchos::RCP<TreeVector> Ptemp = Pu->SubVector("energy");
 Teuchos::RCP<CompositeVector> Ptemp_data = Ptemp->data();

 ierr = A2f2p_->Multiply(false, Pf, Pc);
 ASSERT(!ierr);

 for (int c=0; c!=ncells; ++c){
   double p_c = -pres_c[0][c];
   double t_c = -temp_c[0][c];
   Pc.SumIntoMyValue(c, 0, 0, p_c);
   Pc.SumIntoMyValue(c, 1, 0, t_c);
 }

 Pc.Scale(-1.0);

 for (int c=0; c!=ncells; ++c){
   Ppressure_c[0][c] = Cell_Couple_Inv_[c](0,0)*Pc[0][2*c] + Cell_Couple_Inv_[c](0,1)*Pc[0][2*c + 1];
   Ptemp_c[0][c] = Cell_Couple_Inv_[c](1,0)*Pc[0][2*c] + Cell_Couple_Inv_[c](1,1)*Pc[0][2*c + 1];
 }

 for (int f=0; f!=nfaces; ++f){
   Ppressure_f[0][f] = Pf[0][2*f];
   Ptemp_f[0][f] = Pf[0][2*f + 1];
 }

 if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
   Teuchos::OSTab tab = getOSTab();
   *out_ << "Preconditioned Updates:" << std::endl;
   *out_ << "  Pp0: " << Ppressure_c[0][0] << " " << Ppressure_f[0][3] << std::endl;
   *out_ << "  Pp1: " << Ppressure_c[0][99] << " " << Ppressure_f[0][500] << std::endl;
   *out_ << "  PT0: " << Ptemp_c[0][0] << " " << Ptemp_f[0][3] << std::endl;
   *out_ << "  PT1: " << Ptemp_c[0][99] << " " << Ptemp_f[0][500] << std::endl;
 }
};

void MPCCoupledFlowEnergy::InitPreconditioner_(Teuchos::ParameterList& prec_plist) {
  Epetra_RowMatrix *Pff = dynamic_cast<Epetra_RowMatrix*>(&*P2f2f_);  

  if (prec_method_ == TRILINOS_ML) {
    ml_plist_ =  prec_plist.sublist("ML Parameters");
    ml_prec_ = Teuchos::rcp(new ML_Epetra::MultiLevelPreconditioner(*P2f2f_, ml_plist_, false));
  } else if (prec_method_ == TRILINOS_ILU) {
    ilu_plist_ = prec_plist.sublist("ILU Parameters");
  } else if (prec_method_ == TRILINOS_BLOCK_ILU) {
    ifp_plist_ = prec_plist.sublist("Block ILU Parameters");
#ifdef HAVE_HYPRE
  } else if (prec_method_ == HYPRE_AMG) {
    // read some boomer amg parameters
    hypre_plist_ = prec_plist.sublist("HYPRE AMG Parameters");
    hypre_ncycles_ = hypre_plist_.get<int>("number of cycles",5);
    hypre_nsmooth_ = hypre_plist_.get<int>("number of smoothing iterations",3);
    hypre_tol_ = hypre_plist_.get<double>("tolerance",0.0);
    hypre_strong_threshold_ = hypre_plist_.get<double>("strong threshold",0.25);
  } else if (prec_method_ == HYPRE_EUCLID) {
    hypre_plist_ = prec_plist.sublist("HYPRE Euclid Parameters");
  } else if (prec_method_ == HYPRE_PARASAILS) {
    hypre_plist_ = prec_plist.sublist("HYPRE ParaSails Parameters");
#endif
  }
}


} //namespace
