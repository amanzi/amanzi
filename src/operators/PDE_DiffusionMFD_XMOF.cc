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

      // std::cout<<c<<": ";
      // for (int j=0; j!= (*mat_data_)[c].cells_materials[0].size(); ++j) 
      //   std::cout<<(*mat_data_)[c].cells_vfracs[0][j]<<" ";
      // std::cout<<"\n";

      (*xmof_ir_)[c].set_materials_data((*mat_data_)[c]);

      //std::cout << "\rProcessing cell " << c << "/" << ncells_owned << "...     ";

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

    PDE_DiffusionMFD::InitDiffusion_(plist);

//  XMOF2D::IRTolerances ir_tolerances;
    ir_tolerances_.dist_eps = 1.0e-15;
    ir_tolerances_.div_eps = 1.0e-6;
    ir_tolerances_.ddot_eps = 1.0e-12;
    ir_tolerances_.vfrac_eps = 1.0e-13;
    ir_tolerances_.ang_eps = 1.0e-12;
    ir_tolerances_.mof_max_iter = 1e3;

    ConstructBaseMesh_();

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
        // std::cout<<"tensor "<<Kc<<"\n";
        
        for (int k=0; k<num_cells; ++k){
          double vol = mini_mesh.get_cell(k).size();
          cm[0] = mini_mesh.get_cell(k).center().x;
          cm[1] = mini_mesh.get_cell(k).center().y;
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
            fnor[j][0] = fc.normal()[0];
            fnor[j][1] = fc.normal()[1];
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
          mfd.MassMatrixInverse(cm, vol, fm, fnor, face_area, Kc, Wff);
          
          Wff_cells_mm_[c][k] = Wff;         
          std::cout<<"wff\n"<<Wff<<"\n";
          //exit(0);
            
        }
      }else{
        std::cout<<"pure cell "<<c<<"\n";
        Wff_cells_mm_[c].resize(1);
        Kc(0, 0) = 1.0;
        ok = mfd.MassMatrixInverse(c, Kc, Wff);
        Wff_cells_mm_[c][0] =  Wff;
        std::cout<<"wff\n"<<Wff<<"\n";
      }
      if (ok == WhetStone::WHETSTONE_ELEMENTAL_MATRIX_FAILED) {
        Errors::Message msg("PDE_DiffusionMFD: unexpected failure in WhetStone.");
        Exceptions::amanzi_throw(msg);
      }

    }

  }

  void PDE_DiffusionMFD_XMOF::UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& flux,
                                             const Teuchos::Ptr<const CompositeVector>& u){

    CreateMassMatrices_();

    // update matrix blocks
    AmanziMesh::Entity_ID_List faces, cells;
    std::vector<int> dirs;

    for (int c = 0; c < ncells_owned; c++) {
      std::cout<<"cell "<<c<<"\n\n";
      mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
      int nfaces = faces.size();
      WhetStone::DenseMatrix Acell(nfaces, nfaces);
      
      if ((*xmof_ir_)[c].get_base_mesh().get_cell(0).has_minimesh()){  // mixed cell
        
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

        WhetStone::DenseMatrix Sii(n_int,n_int), Sbb(n_bnd,n_bnd), Sbi(n_bnd, n_int);

        for (int i=0;i<num_cells;i++){
          const XMOF2D::Cell cell = mini_mesh.get_cell(i);
          int faces_c = cell.nfaces();
          WhetStone::DenseMatrix Sc(faces_c, faces_c);
          WhetStone::DenseMatrix& Wc = Wff_cells_mm_[c][i];
          std::vector<double> area(faces_c), tmp(faces_c);
          for (int j=0;j<faces_c;j++) area[j] = cell.get_face(j).size();

          double alp = 0.;
          for (int i1=0; i1<faces_c; i1++){
            tmp[i1] = 0.;
            for (int j1=0; j1<faces_c; j1++){
              tmp[i1] += Wc(i1,j1)*area[j1];
              alp += Wc(i1,j1)*area[i1]*area[j1];
            }
          }

          for (int i1=0; i1<faces_c; i1++){
            for (int j1=0; j1<faces_c; j1++){
              Sc(i1, j1) = (Wc(i1,j1) - tmp[i1]*tmp[j1]/alp)*area[i1]*area[j1];

              if (cell.get_face(i1).is_boundary() && cell.get_face(j1).is_boundary()){
                int i2 = bnd_faces[cell.get_face_index(i1)] - 1;
                int j2 = bnd_faces[cell.get_face_index(j1)] - 1;
                Sbb(i2,j2) += Sc(i1, j1);
              }else if ((!cell.get_face(i1).is_boundary()) && (!cell.get_face(j1).is_boundary())){
                int i2 = -bnd_faces[cell.get_face_index(i1)] - 1;
                int j2 = -bnd_faces[cell.get_face_index(j1)] - 1;
                Sii(i2,j2) += Sc(i1, j1);
              }else if ((cell.get_face(i1).is_boundary()) && (!cell.get_face(j1).is_boundary())){
                int i2 = bnd_faces[cell.get_face_index(i1)] - 1;
                int j2 = -bnd_faces[cell.get_face_index(j1)] - 1;
                Sbi(i2,j2) += Sc(i1, j1);
              }              
            }
          }
        }
      }else{
        std::vector<double> area(nfaces), tmp(nfaces);
        for (int j=0; j<nfaces; j++) area[j] = mesh_->face_area(faces[j]);
        WhetStone::DenseMatrix& Wc = Wff_cells_mm_[c][1];        
        double alp = 0.;
        for (int i1=0; i1<nfaces; i1++){
          tmp[i1] = 0.;
          for (int j1=0; j1<nfaces; j1++){
            tmp[i1] += Wc(i1,j1)*area[j1];
            alp += Wc(i1,j1)*area[i1]*area[j1];
          }
        }
        for (int i1=0; i1<nfaces; i1++){
          for (int j1=0; j1<nfaces; j1++){
            Acell(i1, j1) = (Wc(i1,j1) - tmp[i1]*tmp[j1]/alp)*area[i1]*area[j1];
          }
        }
      }
        
      local_op_->matrices[c] = Acell;

    }
  }
  

}  // namespace Operators
}  // namespace Amanzi
