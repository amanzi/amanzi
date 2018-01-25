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
#include "WhetStone_typedefs.hh"

namespace Amanzi {
namespace WhetStone {

class VectorPolynomial;

class Projector { 
 public:
  enum class Type {
    L2,
    H1,
    SERENDIPITY_L2
  };

 public:
  Projector(Teuchos::RCP<const AmanziMesh::Mesh> mesh) :
      mesh_(mesh),
      d_(mesh_->space_dimension()) {};
  ~Projector() {};

  // projectors for non-conforming (Crouzier-Raviart) virtual space
  // -- simple implementation of low-order elliptic projectors
  void HarmonicCell_CR1(
      int c, const std::vector<VectorPolynomial>& vf, VectorPolynomial& uc) const;

  void HarmonicFace_CR1(
      int f, const AmanziGeometry::Point& p0,
      const std::vector<VectorPolynomial>& ve, VectorPolynomial& uf) const;

  // -- elliptic projector takes cell-moments on input
  void EllipticCell_CRk(
      int c, int order, const std::vector<VectorPolynomial>& vf,
      const std::shared_ptr<DenseVector>& moments, VectorPolynomial& uc) const {
    GenericCell_CRk_(c, order, vf, Type::H1, false, moments, uc);
  }

  // -- harmonic projector calculates and returns cell-moments
  void HarmonicCell_CRk(
      int c, int order, const std::vector<VectorPolynomial>& vf,
      const std::shared_ptr<DenseVector>& moments, VectorPolynomial& uc) const {
    GenericCell_CRk_(c, order, vf, Type::H1, true, moments, uc);
  }

  // projectors for conforming (Lagrange) virtual space
  // -- L2 projector that calculates and returns cell moments
  void L2Cell_Pk(
      int c, int order, const std::vector<VectorPolynomial>& vf,
      const std::shared_ptr<DenseVector>& moments, VectorPolynomial& uc) const {
    GenericCell_Pk_(c, order, vf, Type::L2, true, moments, uc);
  }

  // -- elliptic projector takes cell-moments on input
  void EllipticCell_Pk(
      int c, int order, const std::vector<VectorPolynomial>& vf,
      const std::shared_ptr<DenseVector>& moments, VectorPolynomial& uc) const {
    GenericCell_Pk_(c, order, vf, Type::H1, false, moments, uc);
  }

  // -- harmonic projector calculates and returns cell-moments
  void HarmonicCell_Pk(
      int c, int order, const std::vector<VectorPolynomial>& vf,
      const std::shared_ptr<DenseVector>& moments, VectorPolynomial& uc) const {
    GenericCell_Pk_(c, order, vf, Type::H1, true, moments, uc);
  }

  // -- L2 projector from serendipidy virtual element
  void L2Cell_SerendipityPk(
      int c, int order, const std::vector<VectorPolynomial>& vf,
      const std::shared_ptr<DenseVector>& moments, VectorPolynomial& uc) const;

 private:
  void GenericCell_CRk_(
      int c, int order, const std::vector<VectorPolynomial>& vf,
      const Type type, bool is_harmonic,
      const std::shared_ptr<DenseVector>& moments, VectorPolynomial& uc) const;

  void GenericCell_Pk_(
      int c, int order, const std::vector<VectorPolynomial>& vf,
      const Type type, bool is_harmonic,
      const std::shared_ptr<DenseVector>& moments, VectorPolynomial& uc) const;

  void CalculateDOFsOnBoundary_Pk_(
      int c, int order, const std::vector<VectorPolynomial>& vf,
      DenseVector& vdof, int i) const;

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  int d_;
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

