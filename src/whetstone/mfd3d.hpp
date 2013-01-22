/*
This is the mimetic discretization component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided Reconstruction.cppin the top-level COPYRIGHT file.

Release name: aka-to.
Author: Konstantin Lipnikov (lipnikov@lanl.gov)
Usage: 
*/

#ifndef __MFD3D_HPP__
#define __MFD3D_HPP__

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
const int WHETSTONE_ELEMENTAL_MATRIX_PASSED = 2;
const int WHETSTONE_ELEMENTAL_MATRIX_FAILED = 4;  // only for unexpected situations

const int WHETSTONE_STABILITY_GENERIC = 1;
const int WHETSTONE_STABILITY_GENERIC_SCALED = 2;
const int WHETSTONE_STABILITY_OPTIMIZED_DMP = 3;
const int WHETSTONE_STABILITY_OPTIMIZED_GEOMETRY = 4;

class MFD3D { 
 public:
  explicit MFD3D(Teuchos::RCP<AmanziMesh::Mesh> mesh);
  ~MFD3D() {};

  // primary methods
  int DarcyMass(int cell, const Tensor& permeability,
               Teuchos::SerialDenseMatrix<int, double>& M);
  int DarcyMassInverse(int cell, const Tensor& permeability,
                       Teuchos::SerialDenseMatrix<int, double>& W);
  int DarcyMassInverseSO(int cell, const Tensor& permeability,
                         Teuchos::SerialDenseMatrix<int, double>& W);
  int DarcyMassInverseHex(int cell, const Tensor& permeability,
                          Teuchos::SerialDenseMatrix<int, double>& W);
  int DarcyMassInverseDiagonal(int cell, const Tensor& permeability,
                               Teuchos::SerialDenseMatrix<int, double>& W);
  int DarcyMassInverseOptimized(int cell, const Tensor& permeability,
                                Teuchos::SerialDenseMatrix<int, double>& W);

  int ElasticityStiffness(int cell, const Tensor& deformation,
                          Teuchos::SerialDenseMatrix<int, double>& A); 

  // suppporting primary methods
  int L2consistency(int cell, const Tensor& T,
                    Teuchos::SerialDenseMatrix<int, double>& N,
                    Teuchos::SerialDenseMatrix<int, double>& Mc);
  int L2consistencyInverse(int cell, const Tensor& permeability,
                           Teuchos::SerialDenseMatrix<int, double>& R,
                           Teuchos::SerialDenseMatrix<int, double>& Wc);
  int H1consistency(int cell, const Tensor& T,
                    Teuchos::SerialDenseMatrix<int, double>& N,
                    Teuchos::SerialDenseMatrix<int, double>& Mc);
  int H1consistencyElasticity(int cell, const Tensor& T,
                              Teuchos::SerialDenseMatrix<int, double>& N,
                              Teuchos::SerialDenseMatrix<int, double>& Ac);

  void GrammSchmidt(Teuchos::SerialDenseMatrix<int, double>& N);

  double CalculateStabilityScalar(Teuchos::SerialDenseMatrix<int, double>& Mc);
  void ModifyStabilityScalingFactor(double factor);

  void StabilityScalar(int cell,
                       Teuchos::SerialDenseMatrix<int, double>& N,  // use R, Wc, and W for the inverse matrix
                       Teuchos::SerialDenseMatrix<int, double>& Mc,
                       Teuchos::SerialDenseMatrix<int, double>& M);
  int StabilityMonotoneHex(int cell, const Tensor& T,
                           Teuchos::SerialDenseMatrix<int, double>& Mc,
                           Teuchos::SerialDenseMatrix<int, double>& M);
  int StabilityOptimized(const Tensor& T,
                         Teuchos::SerialDenseMatrix<int, double>& N,
                         Teuchos::SerialDenseMatrix<int, double>& Mc,
                         Teuchos::SerialDenseMatrix<int, double>& M);

  // MFD extension of VAG scheme
  void CalculateHarmonicPoints(int face, std::vector<Tensor>& T, 
                               AmanziGeometry::Point& harmonic_point,
                               double& harmonic_point_weight);

  int DispersionCornerFluxes(int node, int cell, Tensor& dispersion,
                             std::vector<AmanziGeometry::Point>& corner_points,
                             double cell_value,
                             std::vector<double>& corner_values,
                             std::vector<double>& corner_fluxes);

  // extension of mesh API
  int cell_get_face_adj_cell(const int cell, const int face);

  // access members
  double get_scaling_factor() { return scaling_factor_; }
  double get_scalar_stability() { return scalar_stability_; }
  Teuchos::SerialDenseMatrix<int, double>& get_matrix_stability() { return matrix_stability_; }

 private:
  int FindPosition(int v, AmanziMesh::Entity_ID_List nodes);
  Teuchos::RCP<AmanziMesh::Mesh> mesh_;

  int stability_method_;  
  double scalar_stability_, scaling_factor_;
  Teuchos::SerialDenseMatrix<int, double> matrix_stability_;
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

