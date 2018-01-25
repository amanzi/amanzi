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

#include "PDE_DiffusionMFD_Tracer.hh"
//#include "xmof2D.h"

namespace Amanzi {
namespace Operators {

  void PDE_DiffusionMFD_Tracer::InitDiffusion_(Teuchos::ParameterList& plist){

    PDE_DiffusionMFD::InitDiffusion_(plist);

    cell_surfaces_.resize(ncells_owned);

  }

  void PDE_DiffusionMFD_Tracer::UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& flux,
                                               const Teuchos::Ptr<const CompositeVector>& u,
                                               const Teuchos::Ptr<const CompositeVector>& surf_presence,
                                               const Teuchos::Ptr<const CompositeVector>& surface_param){

    // update matrix blocks
    WhetStone::MFD3D_Diffusion mfd(mesh_);
    mfd.ModifyStabilityScalingFactor(factor_);

    AmanziMesh::Entity_ID_List nodes;

    const Epetra_MultiVector& surf_presence_vec = *surf_presence->ViewComponent("cell", true);
    const Epetra_MultiVector& surface_param_vec = *surface_param->ViewComponent("cell", true);

    int num_surfaces = surf_presence_vec.NumVectors();
    
    WhetStone::Tensor Kc(mesh_->space_dimension(), 1);
    Kc(0, 0) = 1.0;

    for (int c = 0; c < ncells_owned; c++) {
    //for (int c = 0; c < 1; c++) {      
      if (K_.get()) Kc = (*K_)[c];
      
      mesh_->cell_get_nodes(c, &nodes);
      int nnodes = nodes.size();

      WhetStone::DenseMatrix Acell(nnodes, nnodes);

      Acell = 0.;

      for (int s_id=0; s_id<num_surfaces;s_id++){
        if (surf_presence_vec[s_id][c] > 0){
        
          //AmanziGeometry::Point sur_cntr(3);//(0.5, 0.5, 0.5);
          AmanziGeometry::Point sur_norm(3);//(1.0, 0.6, 1.0);
          double surf_d = surface_param_vec[s_id*4+3][c];

          sur_norm[0] = surface_param_vec[s_id*4][c];
          sur_norm[1] = surface_param_vec[s_id*4+1][c];
          sur_norm[2] = surface_param_vec[s_id*4+2][c];


          //sur_norm[0]=1.; sur_norm[1]=0.6; sur_norm[2]=1.;  surf_d=-1.3;
          
          CellSurfaceInterception_(c, s_id, sur_norm, surf_d);
          //cell_surfaces_[c][0]->print();          

          mfd.StiffnessMatrixTracer(c, Kc, sur_norm,
                                    cell_surfaces_[c][0]->surface_pnt(),
                                    cell_surfaces_[c][0]->v_ids(),
                                    cell_surfaces_[c][0]->inter_coef(),
                                    Acell);
          break;          
        }else{
          //for (int i=0;i<nnodes;i++) Acell(i,i)=1.;
        }
      }

      // std::cout<<"cell "<<c<<" nodes ";
      // for (int i=0;i<nnodes;i++)std::cout<<nodes[i]<<" "; std::cout<<"\n";
      // std::cout<<Acell<<"\n";      
      local_op_->matrices[c] = Acell;     
    }
    
  }

  void PDE_DiffusionMFD_Tracer::ApplyBCs(bool primary, bool eliminate){

     Teuchos::RCP<BCs> bc_f, bc_n;
      for (std::vector<Teuchos::RCP<BCs> >::iterator bc = bcs_trial_.begin();
          bc != bcs_trial_.end(); ++bc) {
        if ((*bc)->kind() == AmanziMesh::FACE) {
          bc_f = *bc;
        } else if ((*bc)->kind() == AmanziMesh::NODE) {
          bc_n = *bc;
        }
      }
      ApplyBCs_Nodal_(bc_f.ptr(), bc_n.ptr(), primary, eliminate);
      
  }

  void PDE_DiffusionMFD_Tracer::CellSurfaceInterception_(int c, int surf_id,
                                                         const AmanziGeometry::Point& sur_norm, double surf_d){

    AmanziMesh::Entity_ID_List cells, faces, vertices;
    std::vector<int> edges_dirs;
    int dir;
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    Teuchos::RCP<Local_Surface> local_surf = Teuchos::rcp(new Local_Surface(surf_id, sur_norm, surf_d));
    //surface_pnt_[c].resize(nfaces);

    for (int i=0; i<nfaces; i++){
      int f = faces[i];
      AmanziGeometry::Point face_normal = mesh_->face_normal(f, false, c, &dir);
      AmanziGeometry::Point face_centr = mesh_->face_centroid(f);
      AmanziGeometry::Point line_vec = sur_norm^face_normal;
      line_vec *= 1./norm(line_vec);
      int coord0 = 2;
      double eps = 1e-12;
      for (coord0=2; coord0>=0; --coord0){
        if (std::abs(line_vec[coord0]) > eps) break;
      }
      AmanziGeometry::Point line_point(3);

      int coord1, coord2;
      coord1 = (coord0+1)%3;
      coord2 = (coord0+2)%3;
      double a1,a2;
      a1 = sur_norm[coord1]; a2=face_normal[coord1];     
      double b1,b2;
      b1 = sur_norm[coord2]; b2=face_normal[coord2];     
      double d1,d2;
      d1 = surf_d;
      d2 = -face_normal*face_centr;

      double det = a1*b2 - a2*b1;
      line_point[coord0] = 0.;
      line_point[coord1] = (b1*d2 - d1*b2)/ det;
      line_point[coord2] = (a2*d1 - d2*a1)/ det;

      mesh_->face_get_nodes(f, &vertices);

      int nedges = vertices.size();
      int num_inter = 0;
      AmanziGeometry::Point p_tmp(line_point + line_vec);
      
      std::vector<AmanziGeometry::Point> pnts(2);
      std::vector<int> vrt_ids(4);
      std::vector<double> wgt(2);
      
      AmanziGeometry::Point pa(3), pb(3);
                 
      for (int j=0; j<nedges; j++){
        int v_id[2];
        v_id[0] = vertices[j];
        v_id[1] = vertices[(j+1)%nedges];
        AmanziGeometry::Point vs[2];
        for (int k=0; k<2; k++) mesh_->node_get_coordinates(v_id[k], &vs[k]);
        
        double t1, t2, eps=1e-10;
        bool intersection = LineLineIntersect_(vs[0], vs[1], line_point, p_tmp, eps, pa, pb, &t1, &t2);
        if (intersection){   // if possible to find shortest distance
          if ((t1 > 1e-10)&&(t1<1 - 1e-10)){
            if (norm(pa - pb)<eps){
              //std::cout << std::setw(5)<< pa.x()<<" "<<std::setw(5)<<pa.y()<<" "<<std::setw(5)<<pa.z()<<"\n";
              pnts[num_inter] = pa;
              vrt_ids[2*num_inter] = v_id[0];
              vrt_ids[2*num_inter + 1] = v_id[1];
              wgt[num_inter] = t1;
              num_inter++;
            }
          }
        }
        if (num_inter > 2) {
          Errors::Message msg;
          msg << "Surface intersectes a face in more then 2 points\n";
          Exceptions::amanzi_throw(msg);
        }
      }
      if (num_inter == 1) {
        Errors::Message msg;
        msg << "Surface  intersectes only in one point\n Cell "<<c<<"\n";
        Exceptions::amanzi_throw(msg);
      }else if ( num_inter == 2){
        local_surf->Add_Face(pnts, vrt_ids, wgt);
      }
      //std::cout<<"\n";
    }

    cell_surfaces_[c].push_back(local_surf);
      
  }

  /*
   Calculate the shortest segment  between
   two lines P1P2 and P3P4. Calculate also the values of mu_1 and mu_2 where
       int_p1 = P1 + mu_1 (P2 - P1)
       int_p2 = P3 + mu_2 (P4 - P3)
   Return FALSE if no solution exists.
*/

  bool PDE_DiffusionMFD_Tracer::LineLineIntersect_(const AmanziGeometry::Point& p1,
                                                  const AmanziGeometry::Point& p2,
                                                  const AmanziGeometry::Point& p3,
                                                  const AmanziGeometry::Point& p4,
                                                  double EPS,
                                                  AmanziGeometry::Point& int_p1,
                                                  AmanziGeometry::Point& int_p2,
                                                  double *mu_1, double *mu_2){

    AmanziGeometry::Point p13,p43,p21;
    double d1343,d4321,d1321,d4343,d2121;
    double numer,denom;

    p13[0] = p1[0] - p3[0];
    p13[1] = p1[1] - p3[1];
    p13[2] = p1[2] - p3[2];
    p43[0] = p4[0] - p3[0];
    p43[1] = p4[1] - p3[1];
    p43[2] = p4[2] - p3[2];
    if (std::abs(p43[0]) < EPS && std::abs(p43[1]) < EPS && std::abs(p43[2]) < EPS)
      return(false);
    p21[0] = p2[0] - p1[0];
    p21[1] = p2[1] - p1[1];
    p21[2] = p2[2] - p1[2];
    if (std::abs(p21[0]) < EPS && std::abs(p21[1]) < EPS && std::abs(p21[2]) < EPS)
      return(false);

    d1343 = p13[0] * p43[0] + p13[1] * p43[1] + p13[2] * p43[2];
    d4321 = p43[0] * p21[0] + p43[1] * p21[1] + p43[2] * p21[2];
    d1321 = p13[0] * p21[0] + p13[1] * p21[1] + p13[2] * p21[2];
    d4343 = p43[0] * p43[0] + p43[1] * p43[1] + p43[2] * p43[2];
    d2121 = p21[0] * p21[0] + p21[1] * p21[1] + p21[2] * p21[2];

    denom = d2121 * d4343 - d4321 * d4321;
    if (std::abs(denom) < EPS)
      return(false);
    numer = d1343 * d4321 - d1321 * d4343;

    *mu_1 = numer / denom;
    *mu_2 = (d1343 + d4321 * (*mu_1)) / d4343;

    int_p1[0] = p1[0] + *mu_1 * p21[0];
    int_p1[1] = p1[1] + *mu_1 * p21[1];
    int_p1[2] = p1[2] + *mu_1 * p21[2];
    int_p2[0] = p3[0] + *mu_2 * p43[0];
    int_p2[1] = p3[1] + *mu_2 * p43[1];
    int_p2[2] = p3[2] + *mu_2 * p43[2];

    return(true);
    
    
  }

}
}

