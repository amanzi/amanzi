/*
  WhetStone, version 2.1
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Miscaleneous discretization methods for diffusion.
*/

#include <cmath>
#include <vector>

#include "Mesh.hh"
#include "Point.hh"
#include "errors.hh"

#include "MFD3D_Diffusion.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* The conventional FV scheme for a general mesh.
****************************************************************** */
int MFD3D_Diffusion::MassMatrixInverseTPFA(int c, const Tensor& K, DenseMatrix& W)
{
  Entity_ID_List faces;
  std::vector<int> dirs;

  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
  int nfaces = faces.size();

  W.Reshape(nfaces, nfaces);
  W.PutScalar(0.0);

  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
  AmanziGeometry::Point a(d_);

  for (int n = 0; n < nfaces; n++) {
    int f = faces[n];
    const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);

    a = xf - xc;
    double s = mesh_->face_area(f) * dirs[n] / norm(a);
    double Knn = ((K * a) * normal) * s;
    double dxn = a * normal;
    W(n, n) = Knn / fabs(dxn);
  }

  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* The one-sided transmissibility coefficient. Any change to this 
* routine must be consistent with the above routine.
****************************************************************** */
double MFD3D_Diffusion::Transmissibility(int f, int c, const Tensor& K)
{
  int dir;
  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
  AmanziGeometry::Point a(d_);

  const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
  const AmanziGeometry::Point& normal = mesh_->face_normal(f, false, c, &dir);

  a = xf - xc;
  double s = mesh_->face_area(f) * dir / norm(a);
  double Knn = ((K * a) * normal) * s;
  double dxn = a * normal;
  double W = Knn / fabs(dxn);

  return W;
}


/* ******************************************************************
* The debug version of the above FV scheme for a scalar tensor and
* an orthogonal brick element.
****************************************************************** */
int MFD3D_Diffusion::MassMatrixInverseDiagonal(int c, const Tensor& K, DenseMatrix& W)
{
  double volume = mesh_->cell_volume(c);

  Entity_ID_List faces;
  mesh_->cell_get_faces(c, &faces);
  int nfaces = faces.size();

  W.Reshape(nfaces, nfaces);
  W.PutScalar(0.0);

  for (int n = 0; n < nfaces; n++) {
    int f = faces[n];
    double area = mesh_->face_area(f);
    W(n, n) = nfaces * K(0, 0) * area * area / (d_ * volume);
  }
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Second-generation MFD method as inlemented in RC1.
****************************************************************** */
int MFD3D_Diffusion::MassMatrixInverseSO(int c, const Tensor& K, DenseMatrix& W)
{
  Entity_ID_List faces;
  std::vector<int> fdirs;

  mesh_->cell_get_faces_and_dirs(c, &faces, &fdirs);
  int num_faces = faces.size();

  Entity_ID_List nodes, corner_faces;
  mesh_->cell_get_nodes(c, &nodes);
  int nnodes = nodes.size();

  Tensor Kinv(K);
  Kinv.Inverse();

  // collect all corner matrices
  std::vector<Tensor> Mv;
  std::vector<double> cwgt;

  Tensor N(d_, 2), NK(d_, 2), Mv_tmp(d_, 2);

  for (int n = 0; n < nnodes; n++) {
    int v = nodes[n];
    mesh_->node_get_cell_faces(v, c, (ParallelTypeCast)WhetStone::USED, &corner_faces);
    int nfaces = corner_faces.size();
    if (nfaces < d_) {
      Errors::Message msg;
      msg << "WhetStone MFD3D_Diffusion: number of faces forming a corner is small.";
      Exceptions::amanzi_throw(msg);
    }

    for (int i = 0; i < d_; i++) {
      int f = corner_faces[i];
      N.SetColumn(i, mesh_->face_normal(f));
    }
    double cwgt_tmp = fabs(N.Det());

    N.Inverse();
    NK = N * Kinv;

    N.Transpose();
    Mv_tmp = NK * N;
    Mv.push_back(Mv_tmp);

    for (int i = 0; i < d_; i++) {
      int f = corner_faces[i];
      cwgt_tmp /= mesh_->face_area(f);
    }
    cwgt.push_back(cwgt_tmp);
  }

  // rescale corner weights
  double factor = 0.0;
  for (int n = 0; n < nnodes; n++) factor += cwgt[n];
  factor = mesh_->cell_volume(c) / factor;

  for (int n = 0; n < nnodes; n++) cwgt[n] *= factor;

  // assemble corner matrices
  W.Reshape(num_faces, num_faces);
  W.PutScalar(0.0);
  for (int n = 0; n < nnodes; n++) {
    int v = nodes[n];
    mesh_->node_get_cell_faces(v, c, (ParallelTypeCast)WhetStone::USED, &corner_faces);

    Tensor& Mv_tmp = Mv[n];
    for (int i = 0; i < d_; i++) {
      int k = FindPosition_(corner_faces[i], faces);
      for (int j = i; j < d_; j++) {
        int l = FindPosition_(corner_faces[j], faces);
        W(k, l) += Mv_tmp(i, j) * cwgt[n] * fdirs[k] * fdirs[l];
        W(l, k) = W(k, l);
      }
    }
  }
 
  // invert matrix W
  int ierr = W.Inverse();
  if (ierr != 0) {
    Errors::Message msg;
    msg << "WhetStone MFD3D_Diffusion: support operator generated bad elemental mass matrix.";
    Exceptions::amanzi_throw(msg);
  }

  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


int MFD3D_Diffusion::StiffnessMatrixTracer(int c, const Tensor& K, const AmanziGeometry::Point& nG,
                                           std::vector< std::vector<AmanziGeometry::Point> >& surface_pnt,
                                           std::vector< std::vector<int> >& v_ids,
                                           std::vector< std::vector<double> >& inter_coef,
                                           DenseMatrix& A){

  DenseMatrix N, R, A_tmp, q1, N4;
  Entity_ID_List vertices;
  double area;

  mesh_->cell_get_nodes(c, &vertices);
  int nvertices = vertices.size();
  
  N.Reshape(nvertices, 2);
  R.Reshape(nvertices, 2);
  N4.Reshape(nvertices, 4);
  A_tmp.Reshape(nvertices, 2);
  q1.Reshape(nvertices, 1);

  Tensor Kinv;
  Kinv.Init(2,2);

  AmanziGeometry::Point m2(3), m3(3), tmp_p(3);
 
  for (int i=0; i<3; i++){
    if (std::abs(nG[i]) > 1e-10){
      m2[i] = (-nG[(i+1)%3]-nG[(i+2)%3])/nG[i];
      m2[(i+1)%3] = 1.;
      m2[(i+2)%3] = 1.;
      break;
    }
  }
  
  m2 *= 1./norm(m2);
  m3 = m2^nG;

  tmp_p = K*m2;
  Kinv(0,0) = tmp_p*m2;
  Kinv(0,1) = tmp_p*m3;
  tmp_p = K*m3;
  Kinv(1,1) = tmp_p*m3;
  Kinv(1,0) = Kinv(0,1);

  // std::cout<<nG<<"\n"<<m2<<"\n"<<m3<<"\n"<<nG*m2<<" "<<nG*m3<<"\n";
  // std::cout<<"Kinv\n"<<Kinv<<"\n";
  
  Kinv.Inverse();
  
  

  AmanziGeometry::Point vs;  
  for (int i=0; i<nvertices; i++){
    mesh_->node_get_coordinates(vertices[i], &vs);
    N(i,0) = m2*vs; N(i,1) = m3*vs;
    q1(i,0) = nG*vs;
  }

  int nfaces = surface_pnt.size();
  std::vector<AmanziGeometry::Point> nrm(nfaces);
  AmanziGeometry::Point e, cntr(3);
  for (int i=0; i<surface_pnt.size(); i++){
    e = surface_pnt[i][1] - surface_pnt[i][0];
    cntr += 0.5*(surface_pnt[i][1] + surface_pnt[i][0]);
    nrm[i] = e^nG;
    nrm[i] *= 1./norm(nrm[i]); 
  }
  cntr *= 1./nfaces;
  //std::cout << "cntr "<<cntr<<"\n";

  SurfaceElementArea(cntr, surface_pnt, &area);
  
  for (int i=0; i<surface_pnt.size(); i++){
    e = 0.5*(surface_pnt[i][1] + surface_pnt[i][0]);
    if ((e - cntr)*nrm[i] < 0) nrm[i] *= -1.;
    //std::cout<<"nrm "<<nrm[i]<<"\n";
  }


  AmanziGeometry::Point Kn, vrt(3), pa(3);

  R=0.;
  
  for (int i=0; i<nfaces; i++){
    double e_sum = 0.;
    double nKm2, nKm3, e_len;

    Kn = K*nrm[i];
    nKm2 = Kn*m2;
    nKm3 = Kn*m3;
    e_len = norm(surface_pnt[i][1] - surface_pnt[i][0]);
    
    for (int j=0; j<4; j++){
      int k=0;
      while (k < nvertices){
        if (v_ids[i][j]!=vertices[k]){
          k++;
        }else{
          break;
        }
      }
      double wgt;
      if (j%2==0) wgt = 0.5*(1 - inter_coef[i][j/2])*e_len;
      else wgt = 0.5*inter_coef[i][j/2]*e_len;

      e_sum += wgt;

      // Internal verification 
      // mesh_->node_get_coordinates(v_ids[i][j], &vrt);                    
      // pa += wgt*vrt;
      // if ((j==1)||(j==3)){
      //   std::cout << "pa "<<pa<<"  surface_pnt "<<surface_pnt[i][j/2]<<"\n";
      //   pa[0] = 0.;pa[1] = 0.;pa[2] = 0.;
      // }

      R(k,0) += wgt*nKm2;
      R(k,1) += wgt*nKm3;

    }
  }

  // std::cout<<"area "<<area<<"\n";
  // std::cout<<"N\n"<<N<<"\n";
  // std::cout<<"R\n"<<R<<"\n";
  // std::cout<<"q1\n"<<q1<<"\n";

  for (int i=0; i<nvertices; i++){
    for (int j=0; j<2; j++){
      A_tmp(i,j) = 0.;
      for (int k=0; k<2; k++) A_tmp(i,j) += R(i,k)*Kinv(k,j);
    }
  }

  for (int i=0; i<nvertices; i++){
      for (int j=0; j<nvertices; j++){
        A(i,j) = 0.;
        for (int k=0; k<2; k++) A(i,j) += A_tmp(i,k)*R(j,k);
      }
  }

  A *= 1./area;
  
  std::cout<<std::setprecision(10)<<A<<"\n";
  //A = R*Kinv*;

  for (int i=0; i<nvertices; i++){
    N4(i,0) = 1.;
    N4(i,1) = q1(i,0);
    N4(i,2) = N(i,0);
    N4(i,3) = N(i,1);
  }
  //std::cout<<"N4\n"<<N4<<"\n";
  
  StabilityScalar_(N4, A);

  /// Internal check
  // std::cout<<std::setprecision(10)<<A<<"\n";
 

  // AmanziGeometry::Point t1(3), t2(3);
  // std::srand (time(NULL));
  // for (int i=0;i<3;i++) {
  //   t1[i] = 1.*sqrt(c+1)*std::rand()/RAND_MAX;
  //   t2[i] = 1.*sqrt(c+1)*std::rand()/RAND_MAX;
  // }
  // std::cout<<t1<<" "<<t2<<"\n";
  // std::vector<double> r1(nvertices), r2(nvertices);
  // for (int i=0; i<nvertices; i++){
  //   double sum =0.;
  //   for (int j=0; j<nvertices; j++) sum += A(i,j);
  //   std::cout<<i<<": "<<sum<<"\n";
  //   mesh_->node_get_coordinates(vertices[i], &vs);
  //   r1[i] = t1*vs; r2[i] = t2*vs;    
  // }  

  // double val = 0;
  // for (int i=0; i<nvertices; i++){
  //   for (int j=0; j<nvertices; j++) val += A(i,j)*r1[i]*r2[j];
  // }

  // std::cout<<"area "<<area<<"\n";
  // std::cout<<"val1 "<<val<<" val2 "<<t1*t2*area<<" diff "<<val - t1*t2*area<<"\n";
  
  
  //exit(0);
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
  
}


void MFD3D_Diffusion::SurfaceElementArea(AmanziGeometry::Point& cntr,
                                         std::vector< std::vector<AmanziGeometry::Point> >& surface_pnt,
                                         double* area){

  int nedges = surface_pnt.size();
  int dim = cntr.dim();

  *area = 0.;
  for (int i=0;i<nedges;i++){
    AmanziGeometry::Point v1 = surface_pnt[i][0] - cntr;
    AmanziGeometry::Point v2 = surface_pnt[i][1] - cntr;
    (*area) += 0.5*norm(v1^v2);
  } 

}
  
}  // namespace WhetStone
}  // namespace Amanzi



