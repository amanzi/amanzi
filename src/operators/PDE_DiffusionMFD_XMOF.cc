/*
  Operators 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

*/

#include <vector>

// TPLs
#include "Epetra_Vector.h"

// Amanzi
#include "errors.hh"
#include "LinearOperator.hh"
#include "LinearOperatorFactory.hh"
#include "MatrixFE.hh"
#include "MFD3D_CrouzeixRaviart.hh"
#include "MFD3D_Diffusion.hh"
#include "PreconditionerFactory.hh"
#include "SuperMap.hh"
#include "WhetStoneDefs.hh"
#include "Mesh.hh"
// Operators
#include "Op.hh"
#include "Op_Cell_Edge.hh"
#include "Op_Cell_Node.hh"
#include "Op_Cell_FaceCell.hh"
#include "Op_Face_Cell.hh"
#include "Op_SurfaceFace_SurfaceCell.hh"

#include "OperatorDefs.hh"
#include "Operator_Edge.hh"
#include "Operator_FaceCell.hh"
#include "Operator_FaceCellScc.hh"
#include "Operator_FaceCellSff.hh"
#include "Operator_Node.hh"

#include "PDE_DiffusionMFD_XMOF.hh"
#include "GMVMesh.hh"
//#include "xmof2D.h"

namespace Amanzi {
namespace Operators {


  

  void PDE_DiffusionMFD_XMOF::ConstructMiniMesh(const Teuchos::Ptr<const CompositeVector>& vol_fraction,
                                                const Teuchos::Ptr<const CompositeVector>& centroids){

    const Epetra_MultiVector& vol_frac_vec = *vol_fraction->ViewComponent("cell", true);
    const Epetra_MultiVector& centroids_vec = *centroids->ViewComponent("cell", true);

    int num_mat = vol_frac_vec.NumVectors();
    double tol = 1e-15;
    std::vector<int> ncells_with_nmats(1, 0);

    for (int c=0; c != ncells_owned; ++c) {

      (*mat_data_)[c].cells_materials.resize(1);
      (*mat_data_)[c].cells_vfracs.resize(1);
      (*mat_data_)[c].cells_centroids.resize(1);

      int cell_mats = 0;
      for (int j=0;j<num_mat;j++){
        if (vol_frac_vec[j][c] > tol) cell_mats++;
      }
      for (int j=0;j<num_mat;j++){
        if (vol_frac_vec[j][c] > tol){
          (*mat_data_)[c].cells_materials[0].push_back(j);      
          (*mat_data_)[c].cells_vfracs[0].push_back(vol_frac_vec[j][c]);

          if (cell_mats==1) {
            (*mat_data_)[c].cells_centroids[0].push_back(XMOF2D::BAD_POINT);
          }else{
            double xc = centroids_vec[2*j][c];
            double yc = centroids_vec[2*j+1][c];
            (*mat_data_)[c].cells_centroids[0].push_back(XMOF2D::Point2D(xc,yc));
          }
        }
      }

      // std::cout<<c<<": "<<(*mat_data_)[c].cells_materials.size()<<" : ";
      // for (int i=0; i!= (*mat_data_)[c].cells_materials.size(); ++i){
      //   std::cout<<"material "<<i<<"\n";
      //   std::cout<<"mat id ";
      //   for (int j=0; j!= (*mat_data_)[c].cells_materials[i].size(); ++j) 
      //     std::cout<<(*mat_data_)[c].cells_materials[i][j]<<" ";
      //   std::cout<<"\n";
      // //std::cout<<c<<": ";
      //   std::cout<<"vrac ";        
      //   for (int j=0; j!= (*mat_data_)[c].cells_materials[i].size(); ++j) 
      //     std::cout<<(*mat_data_)[c].cells_vfracs[i][j]<<" ";
      //   std::cout<<"\n";
      // //      std::cout<<c<<": ";
      //   std::cout<<"centroid ";
      //   for (int j=0; j!= (*mat_data_)[c].cells_materials[i].size(); ++j) 
      //     std::cout<<(*mat_data_)[c].cells_centroids[i][j]<<" ";
      //   std::cout<<"\n";
      // }

      

      (*xmof_ir_)[c].set_materials_data((*mat_data_)[c]);

      //std::cout<<(*xmof_ir_)[c]<<"\n";

      std::cout << "\r Processing cell " << c << "/" << ncells_owned << "...   ";

      (*xmof_ir_)[c].construct_minimesh(0);

      // if ((*xmof_ir_)[c].get_base_mesh().get_cell(0).has_minimesh()) {
      //    int nmats = (*xmof_ir_)[c].get_base_mesh().get_cell(0).get_minimesh().ncells();
      //    if (nmats > ncells_with_nmats.size())
      //      ncells_with_nmats.resize(nmats, 0);
      //    ncells_with_nmats[nmats - 1]++;
      // }else{
      //   ncells_with_nmats[0]++;
      // }

    }
    

    // std::cout << "Total number of cells: " << ncells_owned << std::endl;
    // std::cout << "Number of single-material cells: " << ncells_with_nmats[0] << std::endl;
    // for (int inm = 0; inm < ncells_with_nmats.size() - 1; inm++)
    //   std::cout << "Number of cells with " << inm + 2 << " materials: " <<
    //     ncells_with_nmats[inm + 1] << std::endl;
    // std::cout << std::endl;


    CreateMassMatrices_();

  }

  void PDE_DiffusionMFD_XMOF::ConstructBaseMesh_(){

    ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
    mesh_cfg_ = Teuchos::rcp(new std::vector<XMOF2D::MeshConfig>(ncells_owned));
    mat_data_ = Teuchos::rcp(new std::vector<XMOF2D::CellsMatData>(ncells_owned));
    xmof_ir_ = Teuchos::rcp(new std::vector<XMOF2D::XMOF_Reconstructor>);


    AmanziMesh::Entity_ID_List nodes;   
    for (int c=0; c != ncells_owned; ++c) {
      mesh_->cell_get_nodes(c, &nodes);
      int nnodes = nodes.size();
      
      (*mesh_cfg_)[c].nodes_coords.resize(nnodes);
      (*mesh_cfg_)[c].ifaces_out_cells.resize(nnodes, 0);
      (*mesh_cfg_)[c].icells_faces.resize(1);
      (*mesh_cfg_)[c].cells_material.resize(1, -1);

      AmanziGeometry::Point ncoord;
      for (int i=0; i<nnodes;i++){
        mesh_->node_get_coordinates(nodes[i], &ncoord);
        (*mesh_cfg_)[c].nodes_coords[i].x = ncoord[0];
        (*mesh_cfg_)[c].nodes_coords[i].y = ncoord[1];       
        (*mesh_cfg_)[c].ifaces_nodes.push_back({i, (i+1)%nnodes});
        (*mesh_cfg_)[c].icells_faces[0].push_back(i);
      }

      // std::cout<<"nodes\n";
      // for (int i=0; i<nnodes;i++){
      //   std::cout<< (*mesh_cfg_)[c].nodes_coords[i].x<<" "<<(*mesh_cfg_)[c].nodes_coords[i].y<<"\n";
      // }
      // std::cout<<"face\n";
      // for (int i=0; i<nnodes;i++){
      //   std::cout<< (*mesh_cfg_)[c].ifaces_nodes[i][0]<<" "<<(*mesh_cfg_)[c].ifaces_nodes[i][1]<<"\n";
      // }
      // std::cout<<"cell2face\n";
      // for (int i=0; i<nnodes;i++) std::cout<< (*mesh_cfg_)[c].icells_faces[0][i]<<" ";
      // std::cout<<"\n";
           

      XMOF2D::XMOF_Reconstructor* ir = new XMOF2D::XMOF_Reconstructor( (*mesh_cfg_)[c],  ir_tolerances_);
      
      (*xmof_ir_).push_back(*ir);
    }

  }

  void PDE_DiffusionMFD_XMOF::InitDiffusion_(Teuchos::ParameterList& plist){

    //PDE_DiffusionMFD::InitDiffusion_(plist);

//  XMOF2D::IRTolerances ir_tolerances;
    ir_tolerances_.dist_eps = 1.0e-15;
    ir_tolerances_.div_eps = 1.0e-6;
    ir_tolerances_.ddot_eps = 1.0e-12;
    ir_tolerances_.vfrac_eps = 1.0e-13;
    ir_tolerances_.ang_eps = 1.0e-12;
    ir_tolerances_.mof_max_iter = 1e3;

    ConstructBaseMesh_();

    // const Epetra_BlockMap& cmap_owned = mesh_->cell_map(false); 
    // rhs_cells_ = Teuchos::rcp(new Epetra_Vector(cmap_owned));

    Abb_.resize(ncells_owned);
    Pb_.resize(ncells_owned);
    Pp_.resize(ncells_owned);
    T3_.resize(ncells_owned);


    for (int c=0; c<ncells_owned;c++) {
      Pp_[c].Reshape(1,1);
      Abb_[c].Reshape(1,1);
      Pb_[c].Reshape(1,1);
      T3_[c].Reshape(1,1);
    }


  }

  
  void PDE_DiffusionMFD_XMOF::CreateMassMatrices_(){

    WhetStone::MFD3D_Diffusion mfd(mesh_);
    int ok = WhetStone::WHETSTONE_ELEMENTAL_MATRIX_FAILED;
    mfd.ModifyStabilityScalingFactor(factor_);
    
    Wff_cells_mm_.resize(ncells_owned);
    WhetStone::Tensor Kc(mesh_->space_dimension(), 1);

    AmanziGeometry::Point cm;
    std::vector< AmanziGeometry::Point > fm;
    std::vector< AmanziGeometry::Point > fnor;
    std::vector< double > face_area;

    for (int c=0; c != ncells_owned; ++c) {
      WhetStone::DenseMatrix Wff;
      if ((*xmof_ir_)[c].get_base_mesh().get_cell(0).has_minimesh()){
        const XMOF2D::MiniMesh& mini_mesh = (*xmof_ir_)[c].get_base_mesh().get_cell(0).get_minimesh();        
        int num_cells=mini_mesh.ncells();
        Wff_cells_mm_[c].resize(num_cells);
        
        for (int k=0; k<num_cells; ++k){

          int mat_id = mini_mesh.get_cell(k).get_material_index();
          Kc = (*KMultiMat_)[c][mat_id];

          double vol = mini_mesh.get_cell(k).size();
          cm[0] = mini_mesh.get_cell(k).center().x;
          cm[1] = mini_mesh.get_cell(k).center().y;
          cm[2] = 0.;
          // std::cout<<"Vol "<< vol <<"\n";                   
          int nfaces = mini_mesh.get_cell(k).nfaces();
          fm.resize(nfaces);
          fnor.resize(nfaces);
          face_area.resize(nfaces);
          // std::cout<<"cell center "<<cm<<"\n";
          for (int j=0; j<nfaces; ++j){
            const XMOF2D::Face fc = mini_mesh.get_cell(k).get_face(j);
            face_area[j] = fc.size();
            fm[j][0] = fc.center().x;
            fm[j][1] = fc.center().y;
            fm[j][2] = 0.;
            fnor[j][0] = fc.normal()[0];
            fnor[j][1] = fc.normal()[1];
            fnor[j][2] = 0.;
            if (!fc.normal_is_out(k)){
              fnor[j][0] *= -1.;
              fnor[j][1] *= -1.;
            }
            // std::cout<<"Face\n"<<fc<<"\n";
            // std::cout<<"nor "<<fnor[j]<<"\n";
            // std::cout<<"area "<<face_area[j]<<"\n";
            // std::cout<<"face center "<<fm[j]<<"\n";
            // std::cout<<"\n";
          }
          Wff=0.;
          ok = mfd.MassMatrixInverse(cm, vol, fm, fnor, face_area, Kc, Wff);
          if (ok == WhetStone::WHETSTONE_ELEMENTAL_MATRIX_FAILED) {
            Errors::Message msg("PDE_DiffusionMFD: unexpected failure in WhetStone.");
            Exceptions::amanzi_throw(msg);
          } 
          Wff_cells_mm_[c][k] = Wff;

          // if ((c==16)||(c==24)) std::cout<<"Wc\n"<<Wff<<"\n";
          // if ((c==16)||(c==24)) std::cout<<"center "<<cm<<"\n";
          // if ((c==16)||(c==24)) std::cout<<"subcell size "<<vol<<" Kc\n"<<Kc<<"\n";
          //std::cout<<"wff\n"<<Wff<<"\n";
          //exit(0);
            
        }
      }else{
        // std::cout<<"pure cell "<<c<<"\n";
        Wff_cells_mm_[c].resize(1);

        Kc = (*KMultiMat_)[c][0];
        ok = mfd.MassMatrixInverse(c, Kc, Wff);
        Wff_cells_mm_[c][0] =  Wff;
        //std::cout<<"wff\n"<<Wff<<"\n";

        if (ok == WhetStone::WHETSTONE_ELEMENTAL_MATRIX_FAILED) {
          Errors::Message msg("PDE_DiffusionMFD: unexpected failure in WhetStone.");
          Exceptions::amanzi_throw(msg);
        }
      }


    }

  }

  void PDE_DiffusionMFD_XMOF::UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& flux,
                                             const Teuchos::Ptr<const CompositeVector>& u, 
                                             double dt){

    //CreateMassMatrices_();


    
    // update matrix blocks
    AmanziMesh::Entity_ID_List faces, cells;
    std::vector<int> dirs;
    int ierr;
    double mixed_mat(0.), pure_mat(0.);
    

    global_op_->rhs()->PutScalarGhosted(0.0);
    Epetra_MultiVector& rhs_faces = *global_op_->rhs()->ViewComponent("face", true);
    //rhs_cells_->PutScalar(0.);
    
    for (int c = 0; c < ncells_owned; c++) {
      // if ((c==16)||(c==24)) std::cout<<"\n"<<"cell "<<c<<"\n\n";
      mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
      int nfaces = faces.size();

      WhetStone::DenseMatrix Scell(nfaces,nfaces);
      
      if ((*xmof_ir_)[c].get_base_mesh().get_cell(0).has_minimesh()){  // mixed cell
        // if ((c==16)||(c==24)) std::cout<<"MIXED\n";
        
        const XMOF2D::MiniMesh& mini_mesh = (*xmof_ir_)[c].get_base_mesh().get_cell(0).get_minimesh();        
        int num_cells = mini_mesh.ncells();
        int num_faces = mini_mesh.nfaces();
        int n_bnd = 0, n_int = 0;
        std::vector<int> bnd_faces(num_faces, 0);
        for (int i=0;i<num_faces;i++) {
          const XMOF2D::Face fc = mini_mesh.get_face(i);         
          //std::cout<<"face "<<i<<" "<<fc.is_boundary()<<"\n";
          if (fc.is_boundary()) {
            n_bnd++;
            bnd_faces[i] = n_bnd;
          }else{
            n_int++;
            bnd_faces[i] = -n_int;
          }
        }

        WhetStone::DenseMatrix Aii(n_int,n_int), Aib(n_int, n_bnd);
        WhetStone::DenseMatrix Pi(n_int, num_cells);
        WhetStone::DenseMatrix M1(n_bnd, n_int), M2(num_cells, n_int);
        WhetStone::DenseMatrix rhs_local(n_bnd, 1), rhs_cell(num_cells, 1);

        Abb_[c].Reshape(n_bnd,n_bnd);
        Pb_[c].Reshape(n_bnd, num_cells);
        Pp_[c].Reshape(num_cells,num_cells);
        T3_[c].Reshape(n_bnd, num_cells);

        Aii = 0.; Abb_[c] = 0.; Aib = 0.;
        Pp_[c] = 0.; Pi = 0.; Pb_[c] = 0.;
        M1 = 0.; M2 = 0.; T3_[c] = 0.;
        rhs_local = 0.; rhs_cell=0.;
        
        for (int i=0;i<num_cells;i++){

          const XMOF2D::Cell cell = mini_mesh.get_cell(i);
          int faces_c = cell.nfaces();
          WhetStone::DenseMatrix Acell(faces_c + 1, faces_c + 1);
          WhetStone::DenseMatrix& Wc = Wff_cells_mm_[c][i];
          //std::vector<double> area(faces_c), tmp(faces_c);
        

          Acell = 0.;
          // Flux elimination
          double matsum = 0.0; 
          for (int n = 0; n < faces_c; n++) {
            double rowsum = 0.0;
            for (int m = 0; m < faces_c; m++) {
              double tmp = Wc(n, m);
              rowsum += tmp;
              Acell(n, m) = tmp;
            }
            Acell(n, faces_c) = -rowsum;
            matsum += rowsum;
          }
          Acell(faces_c, faces_c) = matsum;
        
          for (int n = 0; n < faces_c; n++) {
            double colsum = 0.0;
            for (int m = 0; m < faces_c; m++) colsum += Acell(m, n);
            Acell(faces_c, n) = -colsum;
          }
          

          if (std::isnan(Acell.Norm2())){
            std::cout<<"acell \n"<<Acell<<"\n";
            std::cout<<"c "<<c<<" i "<<i<<"\n";
            std::cout<<Wff_cells_mm_[c][i]<<"\n";
            exit(0);
          }

          for (int i1=0; i1<faces_c; i1++){
            for (int j1=0; j1<faces_c; j1++){             
        
              if (cell.get_face(i1).is_boundary() && cell.get_face(j1).is_boundary()){
                int i2 = bnd_faces[cell.get_face_index(i1)] - 1;
                int j2 = bnd_faces[cell.get_face_index(j1)] - 1;
                Abb_[c](i2,j2) += Acell(i1, j1);
              }else if ((!cell.get_face(i1).is_boundary()) && (!cell.get_face(j1).is_boundary())){
                int i2 = -bnd_faces[cell.get_face_index(i1)] - 1;
                int j2 = -bnd_faces[cell.get_face_index(j1)] - 1;
                Aii(i2,j2) += Acell(i1, j1);
              }else if ((cell.get_face(i1).is_boundary()) && (!cell.get_face(j1).is_boundary())){
                int i2 = bnd_faces[cell.get_face_index(i1)] - 1;
                int j2 = -bnd_faces[cell.get_face_index(j1)] - 1;
                Aib(j2,i2) += Acell(i1, j1);
              }              
            }
            
            if (cell.get_face(i1).is_boundary()) {
              int i2 = bnd_faces[cell.get_face_index(i1)] - 1;
              Pb_[c](i2, i) = Acell(i1, faces_c);
            }else{
              int i2 = -bnd_faces[cell.get_face_index(i1)] - 1;
              Pi(i2, i) = Acell(i1, faces_c);
            }                                       
          }
          //std::cout<<"Acell "<<i<<" "<<Acell(faces_c, faces_c)<<"\n";
          
          Pp_[c](i,i) = Acell(faces_c, faces_c);
          // Accumulation term
          if ((u != Teuchos::null)&&dt > 0) {
            const Epetra_MultiVector& u_cell = *u->ViewComponent("cell", false);
            Pp_[c](i,i) += cell.size()/dt;
            //rhs_cell(i,0) = u_cell[0][c]*cell.size()*D/dt;
            int mat_id = cell.get_material_index();
            rhs_cell(i,0) = u_cell[mat_id + 1][c]*cell.size()/dt;

          }

        }

        Aii.Inverse();

        ierr = M1.Multiply(Aib, Aii, true);
        ierr = M2.Multiply(Pi, Aii, true);

        for (int i1=0; i1<n_bnd; i1++){
          for (int j1=i1; j1<n_bnd; j1++){
            double s = 0.;
            for (int k1=0; k1<n_int; k1++) s += M1(i1, k1)*Aib(k1, j1);
            Abb_[c](i1,j1) -= s; // Abb_[c] - Aib^T*Aii^-1*Aib
            Abb_[c](j1,i1) = Abb_[c](i1,j1);
          }
          for (int j1=0; j1<num_cells; j1++){
            double s = 0.;
            for (int k1=0; k1<n_int; k1++) s += M1(i1, k1)*Pi(k1, j1);
            Pb_[c](i1,j1) -= s; // Pb - Aib^T*Aii^-1*Pi
          }
        }

        for (int i1=0; i1<num_cells; i1++){
          for (int j1=0; j1<num_cells; j1++){
            double s = 0.;
            for (int k1=0; k1<n_int; k1++) s += M2(i1, k1)*Pi(k1, j1);
            Pp_[c](i1,j1) -= s;      //Pp - Pi^T*Aii^-1*Pi
          }
        }                                            

        // std::cout<<"Abb_[c]\n"<<Abb_[c]<<"\n";
        // std::cout<<"Pb\n"<<Pb_[c]<<"\n";
        Pp_[c].Inverse();
        
        T3_[c].Multiply(Pb_[c], Pp_[c], false);
        
        for (int i1=0; i1<n_bnd; i1++){
          for (int j1=i1; j1<n_bnd; j1++){
            double s = 0.;
            for (int k1=0; k1<num_cells; k1++){
              s +=  T3_[c](i1, k1) * Pb_[c](j1, k1);   // Abb - Pb*Pp^-1*Pb^T
            }
            Abb_[c](i1,j1) -= s;
            Abb_[c](j1,i1)  = Abb_[c](i1,j1);
          }
        }

        // std::cout<<"T3\n"<<T3_[c]<<"\n";
        // std::cout<<"rhs_cell\n"<<rhs_cell<<"\n";
        
        rhs_local.Multiply(T3_[c], rhs_cell, false);
        rhs_local *= -1;
        
        // std::cout<<"*************\n";
        // std::cout<<"Abb_[c]\n"<<Abb_[c]<<"\n";

        
        // for (int i=0;i<num_faces;i++) {
        //   std::cout<<bnd_faces[i]<<" ";
        //   if (bnd_faces[i] > 0) {
        //     int par_id = mini_mesh.get_face(i).iparent();
        //     std::cout << par_id  <<" "<< mini_mesh.get_face(i).center()<<" *** "<<
        //       (*xmof_ir_)[c].get_base_mesh().get_face(par_id).center()<<"\n";
        //   }
        //   std::cout<<"\n";
        // }
        Scell = 0.;
        for (int i=0; i<num_faces; i++) {
          if (bnd_faces[i] > 0) {
           
            int i1 = bnd_faces[i] - 1;
            int par_i = mini_mesh.get_face(i).iparent();
            for (int j=0; j<num_faces; j++) {
              if (bnd_faces[j] > 0) {
                int j1 = bnd_faces[j] - 1;
                int par_j = mini_mesh.get_face(j).iparent();
                //std::cout<<i1<<" "<<j1<<" "<<par_i<<" "<<par_j<<"\n";
                Scell(par_i, par_j) += Abb_[c](i1, j1);
              }
            }
            // Add to RHS
            int f = faces[par_i];
            rhs_faces[0][f] += rhs_local(i1,0);            
          }          
        }

        // if ((c==16)||(c==24)) std::cout<<rhs_local<<"\n";
        //std::cout<<"Scell\n"<<Scell<<"\n";
        //mixed_mat += Scell.Norm2();
        
        
      }else{
        // std::vector<double> area(nfaces), tmp(nfaces);
        // for (int j=0; j<nfaces; j++) area[j] = mesh_->face_area(faces[j]);
        T3_[c].Reshape(nfaces,1);
        Pp_[c].Reshape(1,1);
        WhetStone::DenseMatrix& Wc = Wff_cells_mm_[c][0];
        WhetStone::DenseMatrix rhs_local(nfaces,1);
        double rhs_c = 0.;
        
        double alp = 0.;
        for (int i1=0; i1<nfaces; i1++){
          T3_[c](i1,0) = 0.;
          for (int j1=0; j1<nfaces; j1++){
            T3_[c](i1,0) += Wc(i1,j1);//*area[j1];
            alp += Wc(i1,j1);//*area[i1]*area[j1];
          }
        }
        
        // Accumulation term
        if ((u != Teuchos::null)&&dt > 0) {
          alp += mesh_->cell_volume(c)/dt;
          const Epetra_MultiVector& u_cell = *u->ViewComponent("cell", false);
          rhs_c = u_cell[0][c]*mesh_->cell_volume(c)/dt;
        }
        
        
        T3_[c] *= -(1./alp);
        Pp_[c] = (1/alp);

        rhs_local = T3_[c];
        rhs_local *= -rhs_c;
        
        for (int i1=0; i1<nfaces; i1++){
          for (int j1=0; j1<nfaces; j1++){
            Scell(i1, j1) = (Wc(i1,j1) - alp*T3_[c](i1,0)*T3_[c](j1,0));
          }
          // Add to RHS
          int f = faces[i1];
          rhs_faces[0][f] += rhs_local(i1,0);            
        }
        
        //pure_mat += Scell.Norm2();
        //std::cout<<rhs_local<<"\n";
      }
      

      // if ((c==16)||(c==24)) std::cout<<"Scell\n"<<Scell<<"\n";
      
      local_op_->matrices[c] = Scell;

    }

    std::cout<<"mixed "<<mixed_mat<<" pure "<<pure_mat<<"\n";

  }

 
  void PDE_DiffusionMFD_XMOF::ApplyBCs(bool primary, bool eliminate){

    ApplyBCs_Mixed_(*bcs_trial_[0], *bcs_test_[0], primary, eliminate);

  }
  
/* ******************************************************************
* Apply BCs on face values.
****************************************************************** */

void PDE_DiffusionMFD_XMOF::ApplyBCs_Mixed_(BCs& bc_trial, BCs& bc_test,
                                            bool primary, bool eliminate){


// apply diffusion type BCs to FACE-CELL system
  AmanziMesh::Entity_ID_List faces;

  const std::vector<int>& bc_model_trial = bc_trial.bc_model();
  const std::vector<int>& bc_model_test = bc_test.bc_model();

  const std::vector<double>& bc_value = bc_trial.bc_value();
  const std::vector<double>& bc_mixed = bc_trial.bc_mixed();

  // ASSERT(bc_model_trial.size() == nfaces_wghost);
  // ASSERT(bc_value.size() == nfaces_wghost);

  //global_op_->rhs()->PutScalarGhosted(0.0);
  Epetra_MultiVector& rhs_face = *global_op_->rhs()->ViewComponent("face", true);

  for (int c = 0; c != ncells_owned; ++c) {
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();
    
    bool flag(true);
    WhetStone::DenseMatrix& Scell = local_op_->matrices[c];

    // essential conditions for test functions
    for (int n = 0; n != nfaces; ++n) {
      int f = faces[n];
      if (bc_model_test[f] == OPERATOR_BC_DIRICHLET) {
        if (flag) {  // make a copy of elemental matrix
          local_op_->matrices_shadow[c] = Scell;
          flag = false;
        }
        for (int m = 0; m < nfaces; m++) Scell(n, m) = 0.0;
      }
    }

    // conditions for trial functions
    for (int n = 0; n != nfaces; ++n) {
      int f = faces[n];
      double value = bc_value[f];

      if (bc_model_trial[f] == OPERATOR_BC_DIRICHLET) {
        // make a copy of elemental matrix for post-processing
        if (flag) {
          local_op_->matrices_shadow[c] = Scell;
          flag = false;
        }

        if (eliminate) { 
          for (int m = 0; m < nfaces; m++) {
            rhs_face[0][faces[m]] -= Scell(m, n) * value;
            Scell(m, n) = 0.0;
          }
        }

        if (primary) {
          rhs_face[0][f] = value;
          Scell(n,n) = 1.0;
        }

      } else if (bc_model_trial[f] == OPERATOR_BC_NEUMANN && primary) {
        
          rhs_face[0][f] -= value * mesh_->face_area(f);
          
      }
    }
    //std::cout<<Scell<<"\n";
  }

  global_op_->rhs()->GatherGhostedToMaster("face", Add);    

}

  void PDE_DiffusionMFD_XMOF::UpdateFlux(const Teuchos::Ptr<CompositeVector>& p,
                                         const Teuchos::Ptr<CompositeVector>& u,
                                         double dt){

    ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
    AmanziMesh::Entity_ID_List faces, cells;
    std::vector<int> dirs;
    int ierr;

    Epetra_MultiVector& p_sol = *p->ViewComponent("cell", true);
    Epetra_MultiVector& lmd_sol = *p->ViewComponent("face", true); 


    double D = 1;

    for (int c=0; c != ncells_owned; ++c) {
      mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
      int nfaces = faces.size();

      WhetStone::DenseMatrix pnew_c, lmd_f;
      
      if ((*xmof_ir_)[c].get_base_mesh().get_cell(0).has_minimesh()){  // mixed cell
        const XMOF2D::MiniMesh& mini_mesh = (*xmof_ir_)[c].get_base_mesh().get_cell(0).get_minimesh();        
        int num_cells = mini_mesh.ncells();
        int num_faces = mini_mesh.nfaces();        

        int n_bnd = 0, n_int = 0;
        std::vector<int> bnd_faces(num_faces, 0);
        WhetStone::DenseMatrix rhs_c(num_cells, 1), tmp(num_cells, 1);

        for (int i=0;i<num_faces;i++) {
          const XMOF2D::Face fc = mini_mesh.get_face(i);         
          //std::cout<<"face "<<i<<" "<<fc.is_boundary()<<"\n";
          if (fc.is_boundary()) {
            n_bnd++;
            bnd_faces[i] = n_bnd;
          }else{
            n_int++;
            bnd_faces[i] = -n_int;
          }
        }

        lmd_f.Reshape(n_bnd, 1);
        pnew_c.Reshape(num_cells, 1);

        rhs_c = 0.;
        if (dt > 0) {
          for (int i=0;i<num_cells;i++) {
            int mat_id = mini_mesh.get_cell(i).get_material_index();
            rhs_c(i,0) = p_sol[mat_id + 1][c] * mini_mesh.get_cell(i).size()/dt;
          }
        }


        for (int i=0;i<num_faces;i++) {
          const XMOF2D::Face fc = mini_mesh.get_face(i);         
          //std::cout<<"face "<<i<<" "<<fc.is_boundary()<<"\n";
          if (fc.is_boundary()) {
            int f = fc.iparent();
            lmd_f(bnd_faces[i] - 1, 0) = lmd_sol[0][faces[f]];
          }
        }

        //if ((c==16)||(c==24)) std::cout << "lmd_f\n"<<lmd_f<<"\n";
        
        ierr = pnew_c.Multiply(T3_[c], lmd_f, true);
        pnew_c *= -1.;

        ierr = tmp.Multiply(Pp_[c], rhs_c, false);
        pnew_c += tmp;


        
        
        p_sol[0][c] = 0.;
        double sum_area = 0.;
        for (int i = 0; i < num_cells; i++){
          int mat_id = mini_mesh.get_cell(i).get_material_index();
          p_sol[mat_id + 1][c] = pnew_c(i,0);

          sum_area += mini_mesh.get_cell(i).size();
          p_sol[0][c] += pnew_c(i,0)*mini_mesh.get_cell(i).size();
        }
        p_sol[0][c] /= sum_area;

      }else{
     
        lmd_f.Reshape(nfaces, 1);
        pnew_c.Reshape(1, 1);

        for (int i=0;i<nfaces; i++) lmd_f(i,0) = lmd_sol[0][faces[i]];

        if (dt > 0){
          p_sol[0][c] = Pp_[c](0,0) * p_sol[0][c]*mesh_->cell_volume(c)/dt;
        }else{
          p_sol[0][c] = 0.;
        }
        
        // std::cout<<lmd_f<<"\n";
        // std::cout<<"T3 \n"<<T3_[c]<<"\n";
        ierr = pnew_c.Multiply(T3_[c], lmd_f, true);
        p_sol[0][c] -=  pnew_c(0,0);
            
      }
    }

  }

  void PDE_DiffusionMFD_XMOF::WriteSpecialGMV(std::string filename,
                                              const Epetra_MultiVector& vol_frac_vec,
                                              const Epetra_MultiVector& solution){

    int dim = mesh_->space_dimension();
    int dim_cell = mesh_->manifold_dimension();
    gmvwrite_openfile_ir_ascii((char*)filename.c_str(), 4, 8);

    unsigned int num_nodes = mesh_->num_entities(AmanziMesh::NODE, AmanziMesh::OWNED);
    std::vector<double> x,y,z;
    double val;
    int num_mat = vol_frac_vec.NumVectors();
    
    AmanziGeometry::Point xc(dim);
    for (int i=0; i<num_nodes; i++) {
      mesh_->node_get_coordinates(i, &xc);
      x.push_back(xc[0]);
      y.push_back(xc[1]);
      if (dim == 3)  
        z.push_back(xc[2]);
      else         
        z.push_back(0.0);
    }

    int *offset_c = new int[ncells_owned+1];
    offset_c[0] = num_nodes;

    int num_mixcell = 0; 
    int ncell_all = ncells_owned;

    for (int c=0; c != ncells_owned; ++c){
      if ((*xmof_ir_)[c].get_base_mesh().get_cell(0).has_minimesh()){ 
        num_mixcell++;
        const XMOF2D::MiniMesh& mini_mesh = (*xmof_ir_)[c].get_base_mesh().get_cell(0).get_minimesh();        
        int nnodes = mini_mesh.nnodes();
        ncell_all += mini_mesh.ncells();

        offset_c[c+1] = nnodes + offset_c[c];

        for (int i=0;i<nnodes;i++){
          const XMOF2D::Node node = mini_mesh.get_node(i);
          const XMOF2D::Point2D xv = node.get_crd();
          x.push_back(xv.x);
          y.push_back(xv.y);
          z.push_back(0.0);
        }
      }else{
        offset_c[c+1] = offset_c[c];
      }
    }

    int nnode_all = x.size();
    double *xd = new double[nnode_all];
    double *yd = new double[nnode_all];
    double *zd = new double[nnode_all];

    for (int i=0; i<nnode_all; i++){
      xd[i] = x[i];
      yd[i] = y[i];
      zd[i] = z[i];
    }

    double *num_mat_type = new double[ncell_all];
    double *solution_data = new double[ncell_all];
    double *materials_id = new double[ncell_all];


    gmvwrite_node_data(&nnode_all, xd, yd, zd);

    gmvwrite_cell_header(&ncell_all);
    unsigned int *xh = new unsigned int[20];

    AmanziMesh::Entity_ID_List nodes;   
    for (int c=0; c != ncells_owned; ++c) {
      mesh_->cell_get_nodes(c, &nodes);
      int nnodes = nodes.size();
      for (int j=0; j<nnodes; j++) xh[j] = nodes[j] + 1;
      gmvwrite_cell_type((char*) "general 1", nnodes, xh);
      solution_data[c] = solution[0][c];
    }
    
    int k = ncells_owned;
    for (int c=0; c != ncells_owned; ++c){
      num_mat_type[c] = 1;
      if ((*xmof_ir_)[c].get_base_mesh().get_cell(0).has_minimesh()){ 
        num_mat_type[c] = 2;
        materials_id[c] = 0;
        const XMOF2D::MiniMesh& mini_mesh = (*xmof_ir_)[c].get_base_mesh().get_cell(0).get_minimesh();        
        int ncells = mini_mesh.ncells();                                     
        for (int j=0;j<ncells;j++){
          num_mat_type[k] = 1;
          const std::vector<int> mini_cell_nodes = mini_mesh.get_cell(j).get_nodes();
          for (int l=0; l<mini_cell_nodes.size(); l++) xh[l] = offset_c[c] + mini_cell_nodes[l] + 1;
          gmvwrite_cell_type((char*) "general 1", mini_cell_nodes.size(), xh);
          int mat_id = mini_mesh.get_cell(j).get_material_index();
          solution_data[k] = solution[mat_id + 1][c];
          materials_id[k] = mat_id + 1;
          k++;
        }
      }else{
        for (int j=0;j<num_mat;j++){
          if (vol_frac_vec[j][c] > 0.9)  {
            materials_id[c] = j+1;
            break;
          }
        }
      }
    }


    gmvwrite_variable_header();
    std::string varname="num_mat";
    gmvwrite_variable_name_data(0, (char*) varname.c_str(), num_mat_type);
    varname = "mat";
    gmvwrite_variable_name_data(0, (char*) varname.c_str(), materials_id);
    varname="solution";
    gmvwrite_variable_name_data(0, (char*) varname.c_str(), solution_data);
    gmvwrite_variable_endvars();


    delete [] xd;
    delete [] yd;
    delete [] zd;
    delete [] xh;
    delete [] num_mat_type;
    delete [] solution_data;
    delete [] materials_id;

    gmvwrite_closefile();

  }
  
}  // namespace Operators
}  // namespace Amanzi
