/*
  WhetStone, version 2.1
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Algorithms underpinning elliptic projectors.
*/

#ifndef AMANZI_WHETSTONE_PROJECTOR_H1_HH_
#define AMANZI_WHETSTONE_PROJECTOR_H1_HH_

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "Point.hh"

#include "DenseMatrix.hh"
#include "Tensor.hh"
#include "VectorPolynomial.hh"
#include "WhetStone_typedefs.hh"

namespace Amanzi {
namespace WhetStone {

class ProjectorH1 { 
 public:
  ProjectorH1(Teuchos::RCP<const AmanziMesh::Mesh> mesh) :
      mesh_(mesh),
      d_(mesh_->space_dimension()) {};
  ~ProjectorH1() {};

  // -- energy-based projector for non-conforming virtual space
  void HarmonicCell_CR1(
      int c, const std::vector<VectorPolynomial>& vf, VectorPolynomial& uc) const;

  void HarmonicCell_CRk(
      int c, int order,
      const std::vector<VectorPolynomial>& vf, VectorPolynomial& uc) const;

  void HarmonicFace_CR1(
      int f, const AmanziGeometry::Point& p0,
      const std::vector<VectorPolynomial>& ve, VectorPolynomial& uf) const;

  // -- energy-based projector for conforming virtual space
  void HarmonicCell_Pk(
      int c, int order,
      const std::vector<VectorPolynomial>& vf, VectorPolynomial& uc) const;

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  int d_;
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

