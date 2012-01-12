/*
This is the mimetic discretization component of the Amanzi code. 
License: BSD
Author: Konstantin Lipnikov (lipnikov@lanl.gov)
Usage:  
*/

#ifndef __WHETSTONE_hpp__
#define __WHETSTONE_hpp__

/*
This is the discretization package, release alpha.

The package uses the formula M = Mc + Ms, where matrix Mc is build from a 
consistency condition (Mc N = R) and matrix Ms is build from a stability 
condition (Ms N = 0), to generate mass and stiffness matrices for a variety 
of physics packages: flow, transport, thermal, and geomechanics. 
The material properties are imbedded into the the matrix Mc. 

Notation used below: M (mass), W (inverse of M), A (stiffness).

IMPORTANT: all matrices must be reshaped before calling member functions.
*/

#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

#include "Mesh.hh"
#include "Point.hh"

#include "tensor.hpp"


namespace Amanzi {
namespace WhetStone {

const int WHETSTONE_ELEMENTAL_MATRIX_OK = 0;
const int WHETSTONE_ELEMENTAL_MATRIX_WRONG = 1;

class MFD3D { 
 public:
  MFD3D(Teuchos::RCP<AmanziMesh::Mesh> mesh) { mesh_ = mesh; };
  ~MFD3D() {};

  // primary methods
  int darcy_mass(int cell,
                 Tensor& permeability,
                 Teuchos::SerialDenseMatrix<int, double>& M);
  int darcy_mass_inverse(int cell,
                         Tensor& permeability,
                         Teuchos::SerialDenseMatrix<int, double>& W);

  int dispersion_corner_fluxes(int node,
                               int cell,
                               Tensor& dispersion,
                               std::vector<AmanziGeometry::Point>& corner_points,
                               double cell_value,
                               std::vector<double>& corner_values,
                               std::vector<double>& corner_fluxes);

  int elasticity_stiffness(int cell,
                           Tensor& deformation,
                           Teuchos::SerialDenseMatrix<int, double>& A); 

  // suppporting primary methods
  int L2_consistency(int cell,
                     Tensor& T,
                     Teuchos::SerialDenseMatrix<int, double>& N,
                     Teuchos::SerialDenseMatrix<int, double>& Mc);
  int L2_consistency_inverse(int cell,
                             Tensor& permeability,
                             Teuchos::SerialDenseMatrix<int, double>& R,
                             Teuchos::SerialDenseMatrix<int, double>& Wc);
  int H1_consistency(int cell,
                     Tensor& T,
                     Teuchos::SerialDenseMatrix<int, double>& N,
                     Teuchos::SerialDenseMatrix<int, double>& Mc);
  int H1_consistency_elasticity(int cell,
                                Tensor& T,
                                Teuchos::SerialDenseMatrix<int, double>& N,
                                Teuchos::SerialDenseMatrix<int, double>& Ac);

  void gramm_schmidt(Teuchos::SerialDenseMatrix<int, double>& N);

  void stability_scalar(int cell,
                        Teuchos::SerialDenseMatrix<int, double>& N,  // use R, Wc, and W for the inverse matrix
                        Teuchos::SerialDenseMatrix<int, double>& Mc,
                        Teuchos::SerialDenseMatrix<int, double>& M);
  void stability_monotone(int cell,
                          Teuchos::SerialDenseMatrix<int, double>& N,
                          Teuchos::SerialDenseMatrix<int, double>& Mc,
                          Teuchos::SerialDenseMatrix<int, double>& M);

  void calculate_harmonic_points(int face, 
                                 std::vector<Tensor>& T, 
                                 AmanziGeometry::Point& harmonic_point,
                                 double& harmonic_point_weight);

   // extension of mesh API
   int cell_get_face_adj_cell(const int cell, const int face);

   // debug methods
   int darcy_mass_inverse_diagonal(int cell,
                                   Tensor& permeability,
                                   Teuchos::SerialDenseMatrix<int, double>& W);

 private:
  int find_position(int v, AmanziMesh::Entity_ID_List nodes);
  Teuchos::RCP<AmanziMesh::Mesh> mesh_;
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

