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

#ifndef __MFD3D_VAG_HH__
#define __MFD3D_VAG_HH__

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

#include "tensor.hh"


namespace Amanzi {
namespace WhetStone {

class MFD3D_VAG { 
 public:
  explicit MFD3D_VAG(Teuchos::RCP<const AmanziMesh::Mesh> mesh) : mesh_(mesh) {};
  ~MFD3D_VAG() {};

  void CalculateHarmonicPoints(int face, std::vector<Tensor>& T, 
                               AmanziGeometry::Point& harmonic_point,
                               double& harmonic_point_weight);

  int DispersionCornerFluxes(int node, int cell, Tensor& dispersion,
                             std::vector<AmanziGeometry::Point>& corner_points,
                             double cell_value,
                             std::vector<double>& corner_values,
                             std::vector<double>& corner_fluxes);

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

