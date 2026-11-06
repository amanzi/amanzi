/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

  Patch-based linear reconstruction of a vector field v from fluxes
  (2-cochains) f = int_f v . n. In cell c, v(x) = v0 + G (x - xc) where
  v0 and G (d + d^2 unknowns) minimize the flux mismatch on faces of the
  cell and of all cells sharing a node with it, subject to exact
  reproduction of fluxes on faces of c. The reconstruction is exact for
  linear fields on meshes with planar faces. An exception is thrown if
  the patch does not determine the unknowns (e.g. too thin a domain).

  Data layout (component "cell"): d components of the mean value v0
  followed by d^2 components of G stored row-wise, G[i][j] = dv_i/dx_j.
*/

#ifndef AMANZI_RECONSTRUCTION_CELL_COCHAIN_HH_
#define AMANZI_RECONSTRUCTION_CELL_COCHAIN_HH_

#include "Teuchos_RCP.hpp"

#include "CompositeVector.hh"
#include "Mesh.hh"
#include "Point.hh"

namespace Amanzi {
namespace Operators {

class ReconstructionCellCochain {
 public:
  ReconstructionCellCochain(Teuchos::RCP<const AmanziMesh::Mesh> mesh);
  ~ReconstructionCellCochain() {};

  // flux data in the "face" component
  void Compute(const CompositeVector& flux);

  // reconstructed vector at point p using the polynomial of cell c
  AmanziGeometry::Point getValue(int c, const AmanziGeometry::Point& p) const;
  AmanziGeometry::Point getMean(int c) const;
  double getGradient(int c, int i, int j) const { return (*data_c_)[d_ + i * d_ + j][c]; }

  Teuchos::RCP<CompositeVector> data() { return data_; }

 private:
  int d_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<CompositeVector> data_;
  Teuchos::RCP<Epetra_MultiVector> data_c_;
};

} // namespace Operators
} // namespace Amanzi

#endif
