/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Derived MPC for flow and energy.  This couples using a block-diagonal coupler.
------------------------------------------------------------------------- */

#include "field_evaluator.hh"
#include "mpc_prec_coupled_flow_energy.hh"
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
     
    

//   std::string methodstring = mpc_plist_.get<std::string>("preconditioning approach");
//   if (methodstring == "diagonal") {
//     method_ = PRECON_DIAGONAL;
//   } else if (methodstring == "upper triangular") {
//     method_ = PRECON_UPPER_TRIANGULAR;
//   } else if (methodstring == "lower triangular") {
//     method_ = PRECON_LOWER_TRIANGULAR;
//   } else if (methodstring == "alternating") {
//     method_ = PRECON_ALTERNATING;
//   } else if (methodstring == "accumulation") {
//     method_ = PRECON_ACCUMULATION;
//   } else {
//     ASSERT(0);
//   }

  damping_ = plist_.get<double>("preconditioner damping", 1.0);
  
  is_matrix_constructed = false;
  decoupled = false;
  
  coupled_pc_ = plist_.sublist("Coupled PC");
//   ml_plist_ = coupled_pc_.sublist("ML Parameters");

  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData("pressure");
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData("temperature");

  ASSERT(temp->size("cell") == pres->size("cell"));

  D_pT_ = Teuchos::rcp(new Epetra_MultiVector(*pres->ViewComponent("cell", false)));

  D_Tp_ = Teuchos::rcp(new Epetra_MultiVector(*temp->ViewComponent("cell", false)));
  
  SymbolicAssembleGlobalMatrices(S);
  
  flow_pk = Teuchos::rcp_dynamic_cast<Amanzi::Flow::Richards>(sub_pks_[0]);
  energy_pk = Teuchos::rcp_dynamic_cast<Amanzi::Energy::TwoPhase>(sub_pks_[1]);
  
  std::cout<<"Welcome to MPCCoupledFlowEnergy\n";

//   exit(0);

};

void MPCCoupledFlowEnergy::SymbolicAssembleGlobalMatrices (const Teuchos::Ptr<State>& S) 
{
    Teuchos::RCP<const CompositeVector> pres = S->GetFieldData("pressure");
//       Teuchos::RCP<const CompositeVector> temp = S->GetFieldData("temperature");
      
    Teuchos::RCP<const Epetra_BlockMap> cmap = pres->map("cell", false);
    Teuchos::RCP<const Epetra_BlockMap> fmap = pres->map("face", false);
    fmap_wghost = pres->map("face", true);
    
    std::cout<<"Cells "<<pres->size("cell")<<" Faces "<<(*fmap).NumGlobalPoints()<<std::endl;
    
    int num_global = (*fmap).NumGlobalPoints(); std::cout<<"num_global "<<num_global<<std::endl;
    int nummypoint = (*fmap).NumMyPoints(); std::cout<<"num_global "<<nummypoint<<std::endl;
//     int *myelem = (*fmap).MyGlobalElements();
    
    
    
    double_fmap = Teuchos::rcp(new Epetra_BlockMap((*fmap).NumGlobalPoints(),
                                                       (*fmap).NumMyPoints(),
                                                       (*fmap).MyGlobalElements(),
                                                        2,
                                                       (*fmap).IndexBase(),
                                                       (*fmap).Comm()
                                    ));
    

    double_cmap = Teuchos::rcp(new Epetra_BlockMap((*cmap).NumGlobalPoints(),
                                                       (*cmap).NumMyPoints(),
                                                       (*cmap).MyGlobalElements(),
                                                        2,
                                                       (*cmap).IndexBase(),
                                                       (*cmap).Comm()
                                  ));



    double_fmap_wghost = Teuchos::rcp(new Epetra_BlockMap((*fmap_wghost).NumGlobalPoints(),
                                                            (*fmap_wghost).NumMyPoints(),
                                                            (*fmap_wghost).MyGlobalElements(),
                                                             2,
                                                            (*fmap_wghost).IndexBase(),
                                                            (*fmap_wghost).Comm()
                                                           ));

//   int avg_entries_row = (mesh_->space_dimension() == 2) ? MFD_QUAD_FACES : MFD_HEX_FACES;
    int avg_entries_row = 6;

    Epetra_CrsGraph cf_graph(Copy, *double_cmap, *double_fmap_wghost, avg_entries_row, false); 
    Epetra_FECrsGraph ff_graph(Copy, *double_fmap, 2*avg_entries_row - 1, false);
//     Epetra_FECrsGraph ff_graph_test(Copy, *fmap, 2*avg_entries_row - 1, false);
// 
    AmanziMesh::Entity_ID_List faces;
    std::vector<int> dirs;
    const int MFD_MAX_FACES = 14;
    int faces_LID[MFD_MAX_FACES];  // Contigious memory is required.
    int faces_GID[MFD_MAX_FACES];

    Teuchos::RCP<const AmanziMesh::Mesh> mesh = S->GetMesh();

    int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

    for (int c=0; c!=ncells; ++c) {
        mesh->cell_get_faces_and_dirs(c, &faces, &dirs);
        int nfaces = faces.size();

        for (int n=0; n!=nfaces; ++n) {
                faces_LID[n] = faces[n];
                faces_GID[n] = double_fmap_wghost->GID(faces_LID[n]);
//                 std::cout<<faces_LID[n]<<" "<<faces_GID[n]<<std::endl;
        }
//         std::cout<<std::endl;
        cf_graph.InsertMyIndices(c, nfaces, faces_LID);
        ff_graph.InsertGlobalIndices(nfaces, faces_GID, nfaces, faces_GID);
//         ff_graph_test.InsertGlobalIndices(nfaces, faces_GID, nfaces, faces_GID);
    }
    cf_graph.FillComplete(*double_fmap, *double_cmap);
    ff_graph.GlobalAssemble();  // Symbolic graph is complete.
//     ff_graph_test.GlobalAssemble();
    
    
// 
//   // create global matrices

//      A2f2p_ = Teuchos::rcp(new Epetra_VbrMatrix(Copy, cf_graph));
//      P2f2f_ = Teuchos::rcp(new Epetra_FEVbrMatrix(Copy, ff_graph));
//      P2f2f_ = Teuchos::rcp(new Epetra_FEVbrMatrix(Copy, *double_fmap, 1));     

//      A2f2p_ = Teuchos::rcp(new Epetra_VbrMatrix(Copy, *double_cmap, 6));
        A2f2p_ = Teuchos::rcp(new Epetra_VbrMatrix(Copy, cf_graph));
//      P2f2f_ = Teuchos::rcp(new Epetra_FEVbrMatrix(Copy, *double_fmap, 11))
        P2f2f_ = Teuchos::rcp(new Epetra_FEVbrMatrix(Copy, ff_graph));
        
//         Epetra_FECrsMatrix test( Copy, ff_graph_test);  
//         test.GlobalAssemble();
//         test.PutScalar(0.0);
     
//      cout<<ff_graph<<endl;
        P2f2f_->GlobalAssemble();
        
        
        
//         P2f2f_->PutScalar(0.0);

     
     std::cout << "double_fmap.NumMyElements() "<<double_fmap->NumMyElements()<<std::endl;
     
//      exit(0);
     
     
     InitPreconditioner(coupled_pc_);


//      P2f2f_->GlobalAssemble();
     
//         Epetra_FEVbrMatrix test_crs(Copy, *double_fmap, *double_fmap, 1);
//         Epetra_FEVbrMatrix test_crs(Copy, *double_fmap, 1);
//         
//         Epetra_Vector h1(*double_fmap);
//         Epetra_Vector h2(*double_fmap);
//         
//         
// //         double test_values[4];
//         for (int i=0;i<250;i++){
//                 double test_values[4];
//                 test_values[0]=10;test_values[1]=0;test_values[2]=0;test_values[3]=10;
//                 faces_GID[0] = i;
//                 P2f2f_->BeginInsertGlobalValues(faces_GID[0],1,faces_GID);
//                 P2f2f_->SubmitBlockEntry(test_values, 2, 2, 2);
//                 P2f2f_->EndSubmitEntries();
// //                 test_crs.InsertGlobalValues(faces_GID[0], 1, &test_values, faces_GID);
//         }
// //         test_crs.FillComplete();
//         P2f2f_->GlobalAssemble();
//         h1.PutScalar(1.0);
//         
//         P2f2f_->Multiply(false, h1, h2);
//         
//         cout<<h2<<endl;
//         cout<<test_crs;
//         exit(0);

//      exit(0);
}

void MPCCoupledFlowEnergy::update_precon(double t, Teuchos::RCP<const TreeVector> up,
        double h) {
        
        StrongMPC::update_precon(t, up, h);
        
       // d pressure_residual / d T
        // -- update the accumulation derivative
        S_next_->GetFieldEvaluator("water_content")
                ->HasFieldDerivativeChanged(S_next_.ptr(), "richards_pk", "temperature");

        // -- get the accumulation deriv
        Teuchos::RCP<const CompositeVector> dwc_dT =
                S_next_->GetFieldData("dwater_content_dtemperature");
        for (int c=0; c!=dwc_dT->size("cell"); ++c) {
                if (!decoupled)
                        (*D_pT_)[0][c] = (*dwc_dT)("cell",c) / h;
                else
                        (*D_pT_)[0][c] = 0.;
        }

        // d temperature_residual / d p
        // -- update the accumulation derivative
        S_next_->GetFieldEvaluator("energy")
                ->HasFieldDerivativeChanged(S_next_.ptr(), "energy_pk", "pressure");

        // -- get the accumulation deriv
        Teuchos::RCP<const CompositeVector> de_dp =
                S_next_->GetFieldData("denergy_dpressure");
        for (int c=0; c!=de_dp->size("cell"); ++c) {
                if (!decoupled)
                        (*D_Tp_)[0][c] = (*de_dp)("cell",c) / h;
                else 
                        (*D_Tp_)[0][c] = 0.;
        }
        
        ASSERT(de_dp->size("cell") == dwc_dT->size("cell"));
        
        int num_cell = de_dp->size("cell");
        
        ComputeShurComplementPK();
        

        if (ml_prec_->IsPreconditionerComputed()) ml_prec_->DestroyPreconditioner();
        ml_prec_->SetParameterList(ml_plist_);
        ml_prec_->ComputePreconditioner();
 

//         exit(0);
};

void MPCCoupledFlowEnergy::ComputeShurComplementPK(){
        
        Teuchos::RCP<const CompositeVector> pres = S_->GetFieldData("pressure");
//       Teuchos::RCP<const CompositeVector> temp = S->GetFieldData("temperature");
      
        Teuchos::RCP<const Epetra_BlockMap> cmap = pres->map("cell", false);
        Teuchos::RCP<const Epetra_BlockMap> fmap = pres->map("face", false);
        
//         Teuchos::RCP<const Epetra_Map> test_map = pres->map("face", false);
        
        Teuchos::RCP<const AmanziMesh::Mesh> mesh = S_->GetMesh();
        
        int ncells = pres->size("cell");
        int cell_GID ; 
        
        AmanziMesh::Entity_ID_List faces;
        std::vector<int> dirs;
        const int MFD_MAX_FACES = 14;
        int faces_LID[MFD_MAX_FACES];  // Contigious memory is required.
        int faces_GID[MFD_MAX_FACES];
        
        
        
//         Epetra_FEVbrMatrix test_crs(Copy, *double_fmap, 1);
//         
//         Epetra_Vector h1(*double_fmap);
//         Epetra_Vector h2(*double_fmap);
//         
//         
// //         double test_values[4];
//         for (int i=0;i<251;i++){
//                 double test_values[4];
//                 Epetra_SerialDenseMatrix TTV(2,2);
//                 test_values[0]=10;test_values[1]=0;test_values[2]=0;test_values[3]=10;
//                 TTV(0,0)=i;TTV(1,0)=0;TTV(0,1)=0;TTV(1,1)=i;
//                 faces_GID[0] = i;
//                 P2f2f_->BeginInsertGlobalValues(faces_GID[0],1,faces_GID);
// //                 P2f2f_->SubmitBlockEntry(test_values, 2, 2, 2);
//                 P2f2f_->SubmitBlockEntry(TTV);
//                 P2f2f_->EndSubmitEntries();
// //                 test_crs.InsertGlobalValues(faces_GID[0], 1, &test_values, faces_GID);
//         }
// //         test_crs.FillComplete();
//         P2f2f_->GlobalAssemble();
//         h1.PutScalar(1.0);
//         
//         P2f2f_->Multiply(false, h1, h2);
//         
// //         cout<<*P2f2f_<<endl;
//         cout<<h2<<endl;
// //     
//         exit(0);

        std::vector<Teuchos::SerialDenseMatrix<int, double> >& flow_Aff = flow_pk->preconditioner_->Aff_cells();
        std::vector<double>& flow_Acc = flow_pk->preconditioner_->Acc_cells();
        std::vector<Epetra_SerialDenseVector>& flow_Afc = flow_pk->preconditioner_->Afc_cells();
        std::vector<Epetra_SerialDenseVector>& flow_Acf = flow_pk->preconditioner_->Acf_cells();
        std::vector<Teuchos::SerialDenseMatrix<int, double> >& energy_Aff = energy_pk->preconditioner_->Aff_cells();
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
//         cout<<(*P2f2f_);
//         exit(0);
        
        for (int c=0; c < ncells; ++c){
                
                cell_GID = cmap->GID(c);
                mesh->cell_get_faces_and_dirs(c, &faces, &dirs);
                
//                 int nfaces = faces_LID.size();
//                 std::cout<<nfaces<<std::endl;
                
                double det_cell_couple = flow_Acc[c] * energy_Acc[c] - (*D_pT_)[0][c] * (*D_Tp_)[0][c];
                
                                
                if (det_cell_couple != 0.) {
                        Couple_Inv(0, 0) = energy_Acc[c]/det_cell_couple;
                        Couple_Inv(1, 1) = flow_Acc[c]/det_cell_couple;
                        Couple_Inv(0, 1) = -(*D_pT_)[0][c]/det_cell_couple;
                        Couple_Inv(1, 0) = -(*D_Tp_)[0][c]/det_cell_couple;
                }
                else {
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
//                         Apf(0,i) =           Couple_Inv(0, 0)*flow_Acf[c](i);
//                         Apf(0,i + nfaces)  = Couple_Inv(0, 1)*energy_Afc[c](i);
//                         Apf(1,i) =           Couple_Inv(1, 0)*flow_Acf[c](i);
//                         Apf(1,i + nfaces)  = Couple_Inv(1, 1)*energy_Afc[c](i);
                        Apf(0,i) =           flow_Acf[c](i);
                        Apf(0,i + nfaces)  = 0;//energy_Afc[c](i);
                        Apf(1,i) =           0;
                        Apf(1,i + nfaces)  = energy_Afc[c](i);
                }
               
                Cell_Couple_Inv_.push_back(Couple_Inv); 
 
//                 std::cout << "Cell "<<c<<std::endl;
//                 std::cout << Schur<<std::endl<<std::endl;
                int NumEntries = nfaces;
                for (int i=0; i!=nfaces; ++i) {
                        faces_LID[i] = faces[i];
                        faces_GID[i] = fmap_wghost->GID(faces_LID[i]);
                }

                for (int i=0; i!=nfaces; ++i) {
                        P2f2f_->BeginSumIntoGlobalValues(faces_GID[i], NumEntries, faces_GID);
//                         cout<<NumEntries<<": "<<faces_GID[i]<<": ";
//                         for (int j=0; j!=nfaces; ++j) cout<<faces_GID[j]<<" ";
//                         cout<<endl;
                        
                        for (int j=0; j!=nfaces; ++j){
                                Values(0,0) = Schur(i,j);
                                Values(0,1) = Schur(i,j + nfaces);
                                Values(1,0) = Schur(i + nfaces,j);
                                Values(1,1) = Schur(i+ nfaces,j+ nfaces);
//                                 Values(0,0) = 10.;
//                                 Values(0,1) = 0.;
//                                 Values(1,0) = 0.;
//                                 Values(1,1) = 10.;
                                P2f2f_->SubmitBlockEntry(Values);
                        }
                        try {
                                P2f2f_->EndSubmitEntries();
                        }
                        catch (int err){
                          cout<<"An Exception occured. P2f2f_ Exception Nr. "<<err<<endl;
                          exit(0);
                        }
                }

//                 cout<<"cell_GID "<<cell_GID<<endl;

                A2f2p_->BeginReplaceGlobalValues(cell_GID, NumEntries, faces_GID);
                for (int i=0; i!=nfaces; ++i) {
                        Values(0,0) = Apf(0,i);
                        Values(0,1) = Apf(0,i + nfaces);
                        Values(1,0) = Apf(1,i);
                        Values(1,1) = Apf(1,i + nfaces);
//                         cout<<"Values "<<c<<"\n"<<Values<<endl;
                        A2f2p_->SubmitBlockEntry(Values);
                }
                try {
                        A2f2p_->EndSubmitEntries();
                 }
                        catch (int err){
                          cout<<"An Exception occured. A2f2p_ Exception Nr. "<<err<<endl;
                          exit(0);
                 }       

        }
 
        A2f2p_->FillComplete(*double_fmap, *double_cmap);
         
        P2f2f_->GlobalAssemble();
        
        is_matrix_constructed = true;

//         cout<<(*P2f2f_);
//         exit(0);

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
  
  
  int ierr;
  
  
//   Epetra_MultiVector cl((u->data)->ViewComponent("cell", false));
//   Epetra_MultiVector fc(*Y->ViewComponent("face", false));

//   Teuchos::RCP<const CompositeVector> pres = (u->SubVector("pressure"))->data();
//   Teuchos::RCP<const CompositeVector> temp = (u->SubVector("temperature"))->data();
  
//   Teuchos::RCP<const CompositeVector> pp = 
//   Teuchos::RCP<Epetra_MultiVector> lp = Teuchos::rcp(new Epetra_MultiVector(*pres->ViewComponent("face", false))); 
//   Teuchos::RCP<Epetra_MultiVector> tt = Teuchos::rcp(new Epetra_MultiVector(*temp->ViewComponent("cell", false))); 
//   Teuchos::RCP<Epetra_MultiVector> lt = Teuchos::rcp(new Epetra_MultiVector(*temp->ViewComponent("face", false))); 
  
  int ncells = pres_c.MyLength();
  int nfaces = pres_f.MyLength();
  double val = 0.;
  

          
//   cout<<temp_c<<endl;
  
  for (int c=0; c < ncells; ++c){
          double p_c = pres_c[0][c];
          double t_c = temp_c[0][c];
          val = -(Cell_Couple_Inv_[c](0,0)*p_c + Cell_Couple_Inv_[c](0,1)*t_c);
          Tc.ReplaceMyValue(c, 0, 0, val);
          val = -(Cell_Couple_Inv_[c](1,0)*p_c + Cell_Couple_Inv_[c](1,1)*t_c);
          Tc.ReplaceMyValue(c, 1, 0, val);
  }
  
//   cout<<Tc<<endl;
//   for (int c=0; c < ncells; ++c){
//           cout<<c<<" "<<Tc[0][2*c]<<endl;
//   }
// //   exit(0);
//   
//   cout<<(*A2f2p_)<<endl;


   A2f2p_->Multiply(true, Tc, Tf); // It performs the required parallel communications.
  
//    cout<<"Tf\n"<<endl;
//    for (int f=0; f < nfaces; ++f){
//            cout<<f<<" "<<Tf[0][2*f+1]<<" "<< temp_f[0][f]<<endl;
//    }
// //    cout<<Tf<<endl;
//    
//    exit(0);
 
  for (int f=0; f < nfaces; ++f){
  		Tf.SumIntoMyValue(f, 0, 0, pres_f[0][f]); 
  		Tf.SumIntoMyValue(f, 1, 0, temp_f[0][f]); 
  }  


//    cout<<"Tf\n"<<endl;
//    for (int f=0; f < nfaces; ++f){
//            cout<<f<<" "<<Tf[0][2*f+1]<<endl;
//    }
// //    cout<<Tf<<endl;
//    
//    exit(0);
  
  ierr = ml_prec_->ApplyInverse(Tf, Pf);
  
  
//   P2f2f_->Multiply(false, Pf, test_res);
//   test_res.Update(-1,Tf, 1);
//   
//   double res;
//   
//   test_res.Norm2(&res);
//   
//   cout<<"Residual : "<<res<<endl;
  
//     cout<<"P2f2f_\n"<<*P2f2f_<<endl;
    
//    cout<<"Pf\n"<<endl;
//    for (int f=0; f < nfaces; ++f){
//            cout<<f<<" "<<Pf[0][2*f+1]<<endl;
//    }
    
//     exit(0);

  Teuchos::RCP<TreeVector> Ppressure = Pu->SubVector("flow");
  Teuchos::RCP<CompositeVector> Ppressure_data = Ppressure->data();
  
  Teuchos::RCP<TreeVector> Ptemp = Pu->SubVector("energy");
  Teuchos::RCP<CompositeVector> Ptemp_data = Ptemp->data();
  
//   for (int c=0;c < ncells; ++c) (*Ppressure_data)("cell",c) = 100.;
//   
//   cout<<*Ppressure_data->ViewComponent("cell", false)<<endl;
  
//   exit(0);
  
//   Teuchos::RCP<Epetra_MultiVector> Ppres_c = Ppressure_data->ViewComponent("cell", false);
//   Teuchos::RCP<Epetra_MultiVector> Ptemp_c = temp_data->ViewComponent("cell", false);
//   
//   Teuchos::RCP<Epetra_MultiVector> Ppres_f = pressure_data->ViewComponent("face", false);
//   Teuchos::RCP<Epetra_MultiVector> Ptemp_f = temp_data->ViewComponent("face", false);
  
//   cout<<"Pf\n"<<endl;
//    for (int f=0; f < nfaces; ++f){
//            cout<<f<<" "<<Pf[0][2*f]<<endl;
//    }
  
//   cout<<*A2f2p_<<endl;
  
  ierr = A2f2p_->Multiply(false, Pf, Pc);
  
//   cout<<"Pc\n"<<endl;
//    for (int c=0; c < ncells; ++c){
//            cout<<c<<" "<<Pc[0][2*c]<<endl;
//    }
//    exit(0);
  
  
  for (int c=0; c < ncells; ++c){
          double p_c = -pres_c[0][c];
          double t_c = -temp_c[0][c];
          Pc.SumIntoMyValue(c, 0, 0, p_c);
          Pc.SumIntoMyValue(c, 1, 0, t_c);
  }

  Pc.Scale(-1.0);
  
//    cout<<"Pc\n"<<endl;
//    for (int c=0; c < ncells; ++c){
//            cout<<c<<" "<<Pc[0][2*c]<<endl;
//    }
//    exit(0);

  for (int c=0; c < ncells; ++c){
        (*Ppressure_data)("cell",c) = Cell_Couple_Inv_[c](0,0)*Pc[0][2*c] + Cell_Couple_Inv_[c](0,1)*Pc[0][2*c + 1];
        (*Ptemp_data)("cell",c) = Cell_Couple_Inv_[c](1,0)*Pc[0][2*c] + Cell_Couple_Inv_[c](1,1)*Pc[0][2*c + 1];
  }
  
  for (int f=0; f < nfaces; ++f){
       (*Ppressure_data)("face",f) = Pf[0][2*f];
       (*Ptemp_data)("face",f) = Pf[0][2*f + 1];   
  }
  
//   Ppressure_data->SetComponent("cell", Ppres_c);
//   Ppressure_data->SetComponent("face", Ppres_f);
  
//   for (int c=0; c < ncells; ++c){
//        cout<<pres_c[0][c]<<" "<<(*Ppres_c)[0][c]<<endl;
//   }
//     cout<<*(((Pu->SubVector("flow"))->data())->ViewComponent("cell", false));
    
// exit(0);

};

void MPCCoupledFlowEnergy::InitPreconditioner(Teuchos::ParameterList& prec_plist) {
        
   int NumPDEEqns = 5;

   int i = 900;
   
   cout<<"Start ML example\n";
   
      
   Epetra_RowMatrix *Pff = dynamic_cast<Epetra_RowMatrix*>(&*P2f2f_);  
   
//   if (prec_method_ == TRILINOS_ML) {
    ml_plist_ =  prec_plist.sublist("ML Parameters");
//     ml_prec_ = Teuchos::rcp(new ML_Epetra::MultiLevelPreconditioner(*P2f2f_, ml_plist_));
    ml_prec_ = Teuchos::rcp(new ML_Epetra::MultiLevelPreconditioner(*P2f2f_, ml_plist_, false));
    
    
//     ml_prec_ = Teuchos::rcp(new ML_Epetra::MultiLevelPreconditioner(*Pff));
    
    cout<<"Prec is created\n";
    

}

double MPCCoupledFlowEnergy::enorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> du){
        
        double norm ;
        
        norm = StrongMPC::enorm(u, du);
        
        cout<<"Enorm "<<norm<<endl;
        
        return norm;
}

} //namespace
