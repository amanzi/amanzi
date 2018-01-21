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

  }

  void PDE_DiffusionMFD_Tracer::UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& flux,
                                               const Teuchos::Ptr<const CompositeVector>& u){

    AmanziGeometry::Point sur_cntr(0.5, 0.5, 0.5);
    AmanziGeometry::Point sur_norm(1.0, 0.6, 1.0);

    CellSurfaceInterception_(0, sur_cntr, sur_norm);
    
  }

  void PDE_DiffusionMFD_Tracer::CellSurfaceInterception_(int c, 
                                                         const AmanziGeometry::Point& sur_cntr, 
                                                         const AmanziGeometry::Point& sur_norm){

    AmanziMesh::Entity_ID_List cells, faces, vertices;
    std::vector<int> edges_dirs;
    int dir;
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();
    
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
      d1 = -sur_norm*sur_cntr;
      d2 = -face_normal*face_centr;

      double det = a1*b2 - a2*b1;
      line_point[coord0] = 0.;
      line_point[coord1] = (b1*d2 - d1*b2)/ det;
      line_point[coord2] = (a2*d1 - d2*a1)/ det;

      mesh_->face_get_nodes(f, &vertices);

      int nedges = vertices.size();
      int num_inter = 0;
      AmanziGeometry::Point p_tmp(line_point + line_vec);
      
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
              std::cout << std::setw(5)<< pa.x()<<" "<<std::setw(5)<<pa.y()<<" "<<std::setw(5)<<pa.z()<<"\n";
              num_inter++;
            }
          }
        }
        if (num_inter > 2) {
          Errors::Message msg;
          msg << "Surface intercent a face in more then 2 points\n";
          Exceptions::amanzi_throw(msg);
        }
      }
      std::cout<<"\n";
    }    
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

