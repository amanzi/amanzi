/*
 * Elements.hh
 *
 *  Created on: Oct 15, 2020
 *      Author: l328622
 */

#ifndef SRC_PKS_SHALLOW_WATER_ELEMENTS_HH_
#define SRC_PKS_SHALLOW_WATER_ELEMENTS_HH_

namespace Amanzi {
namespace ShallowWater {

class Element2DQuad {

public:

  void Elements(){};
  void ~Elements(){};

  double basis(int,std::vector<double>&);
  double basis_grad(int,std::vector<double>&);
  void quadrature();

  int nquad_vol, nquad_edge;
  std::vector<std::vector<double> > quad_nodes_vol;
  std::vector<double> weight_vol;
  std::vector<std::vector<double> > quad_nodes_edge;
  std::vector<double> weight_edge;

};

double Element2DQuad::basis(int k, std::vector<double>& x){
  switch (k) {
    case 0:
      return x[0]*x[1];
    case 1:
      return (1.-x[0])*x[1];
    case 2:
      return (1.-x[0])*(1.-x[1]);
    case 3:
      return x[0]*(1.-x[1]);
    default:
      std::cout << "P1, basis function " << k << std::endl;
      exit(1);
  }
}

std::vector<double> Element2DQuad::basis_grad(int k, std::vector<double>& x){
  std::vector<double> grad;
  grad.resize(2);

  switch (k) {
    case 0:
     grad[0] = x[1]; grad[1] = x[0];
    case 1:
     grad[0] = -x[1]; grad[1] = 1.-x[0];
    case 2:
     grad[0] = -(1.-x[1]); grad[1] = -(1.-x[0]);
    case 3:
     grad[0] = 1.-x[1]; fy = -x[0];
    default:
      std::cout << "P1, gradient of basis function " << k << std::endl;
      exit(1);
  }
}

std::vector<double> Element2DQuad::EvalSol(std::vector<std::vector<double> >& f, std::vector<double>& x) {

  sol = 0;
  for (int i = 0; i < 3; i++) {
    for (int k = 0; k < nDOFs; k++) {
      sol += f[i][k]*basis[k];
    }
  }

  return sol;
}

void quadrature(int c){

  // quadrature in the cell
  //---- 9-point, 2d product of 3-point Gauss ----!
  nquad_vol = 9;
  int nvertex = 4;
  quad_nodes_vol.resize(nvertex);
  for int (n = 0; n < nquad_vol; n++) quad_nodes_vol[n].resize(nquad_vol);
  weight_vol.resize(nquad_vol);

  std::vector<double> qp, qw;
  qp.resize(3); qw.resize(3);

  double s=sqrt(0.6);
  qp[0] = 0.5*(1.0 - s);
  qp[1] = 0.5;
  qp[2] = 0.5*(1.0 + s);
  qw[0] = 5.0/18.0;
  qw[1] = 8.0/18.0;
  qw[2] = 5.0/18.0;

  int iq = -1;
  for (int iq2 = 0; iq2 < 3; iq2++) {
     for (int iq1 = 0; iq1 < 3; iq1++) {
       iq += 1
       quad_nodes_vol[0][iq] = qp(iq1);
       quad_nodes_vol[1][iq] = qp(iq2);
       weight_vol[iq] = qw(iq1)*qw(iq2);
     }
  }

  for (iq = 0; iq < nquad_vol; iq++) {
      quad_nodes_vol[2][iq] = 1.0 - quad_nodes_vol[0,iq];
      quad_nodes_vol[3][iq] = 1.0 - quad_nodes_vol[1,iq];
  }

  // quadrature on edges
  //---- 3-point Gauss formula ----!
  nquad_edge = 3;

  std::vector<std::vector<double> > qp_edge;
  std::vector<double> qw_edge;
  qp_edge.resize(nquad_edge); qw_edge.resize(nquad_edge);

  qp_edge[0].resize(nquad_edge);
  qp_edge[0][0] = 0.5*(1.0 - s);
  qp_edge[1][0] = 0.5*(1.0 + s);
  qw_edge[0] = 5.0/18.0;
  qp_edge[1].resize(nquad_edge);
  qp_edge[0][1] = 0.5*(1.0 + s);
  qp_edge[1][1] = 0.5*(1.0 - s);
  qw_edge[1] = 5.0/18.0;
  qp_edge[2].resize(nquad_edge);
  qp_edge[0][2] = 0.5;
  qp_edge[1][2] = 0.5;
  qw_edge[2] = 8.0/18.0;

  // triangles
//  mesh_->cell_get_faces(c,&cfaces);

//
//  for (int edgenum = 0; edgenum < cfaces.size(); edgenum++) {
//    for (iq = 0; iq < nquad; iq++) {
//       quad_nodes_edge(edgenum,edgenum,iq)         = qp_edge(1,iq)
//       quad_nodes_edge(edgenum,e%next(edgenum),iq) = qp_edge(2,iq)
//       weight_edge(edgenum,iq)                     = qw_edge(iq)
//    }
//  }

  // quads

  int edgenum;

  edgenum = 0
  for (int iq = 0; iq < nquad_edge; iq++) {
      quad_nodes_edge(edgenum,0,iq)   = quad_edge(0,iq);
      quad_nodes_edge(edgenum,1,iq)   = 0.0;
      quad_nodes_edge(edgenum,2,iq) = 1.0 - e%quad_edge(edgenum,0,iq);
      quad_nodes_edge(edgenum,3,iq) = 1.0 - e%quad_edge(edgenum,1,iq);
  }
  edgenum = 1
  for (int iq = 0; iq < nquad_edge; iq++) {
    quad_nodes_edge(edgenum,0,iq)   = 1.0;
    quad_nodes_edge(edgenum,1,iq)   = quad_edge(0,iq);
    quad_nodes_edge(edgenum,2,iq) = 1.0 - e%quad_edge(edgenum,0,iq);
    quad_nodes_edge(edgenum,3,iq) = 1.0 - e%quad_edge(edgenum,1,iq);
  }
  edgenum = 2
  for (int iq = 0; iq < nquad_edge; iq++) {
    quad_nodes_edge(edgenum,0,iq)   = quad_edge(0,iq);
    quad_nodes_edge(edgenum,1,iq)   = 1.0;
    quad_nodes_edge(edgenum,2,iq) = 1.0 - e%quad_edge(edgenum,0,iq);
    quad_nodes_edge(edgenum,3,iq) = 1.0 - e%quad_edge(edgenum,1,iq);
  }
  edgenum = 3
  for (int iq = 0; iq < nquad_edge; iq++) {
    quad_nodes_edge(edgenum,0,iq)   = 0.0;
    quad_nodes_edge(edgenum,1,iq)   = quad_edge(0,iq);
    quad_nodes_edge(edgenum,2,iq) = 1.0 - e%quad_edge(edgenum,0,iq);
    quad_nodes_edge(edgenum,3,iq) = 1.0 - e%quad_edge(edgenum,1,iq);
  }

  for (edgenum = 0; edgenum < 3; edgenum++) {
    for (int iq = 0; iq < nquad_edge; iq++) {
        weight_edge(edgenum,iq)   = weight_edge(iq);
    }
  }

}



}  // namespace ShallowWater
}  // namespace Amanziss

#endif /* SRC_PKS_SHALLOW_WATER_ELEMENTS_HH_ */
