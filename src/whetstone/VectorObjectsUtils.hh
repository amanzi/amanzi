/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Miscalleneous tools for working with vector objects.
*/ 

#ifndef AMANZI_WHETSTONE_VECTOR_OBJECTS_UTILS_HH_
#define AMANZI_WHETSTONE_VECTOR_OBJECTS_UTILS_HH_

#include "DenseMatrix.hh"
#include "Polynomial.hh"
#include "SpaceTimePolynomial.hh"
#include "VectorObjects.hh"

namespace Amanzi {
namespace WhetStone {

// differential operators
VectorPolynomial Gradient(const Polynomial& p);
VectorSpaceTimePolynomial Gradient(const SpaceTimePolynomial& p);

VectorPolynomial Curl3D(const VectorPolynomial& p);
Polynomial Curl2D(const VectorPolynomial& p);
VectorPolynomial Rot2D(const Polynomial& p);

Polynomial Divergence(const VectorObjects<Polynomial>& vp);

// matrix form of differential operators
DenseMatrix Curl3DMatrix(int d, int order);

// vector decompositions
// -- q_k = curl(p_k ^ x) + x . p_{k-1}
void VectorDecomposition3DCurl(const Monomial& q, int component,
                               VectorPolynomial& p1, Polynomial& p2);
// -- q_k = rot(p_{k+1}) + x . p_{k-1}
void VectorDecomposition2DRot(const VectorPolynomial& q,
                              Polynomial& p1, Polynomial& p2);

//  3D: q_k = grad(p_{k+1}) + x ^ p_{k-1}
//  2D: q_k = grad(p_{k+1}) + x*. p_{k-1}, x* = (-y, x)
// void VectorDecomposition3DGrad(T& p1, VectorObjects<T>& p2) {};
// void VectorDecomposition2DGrad(T& p1, T& p2) {};

// projector
VectorPolynomial ProjectVectorPolynomialOnManifold(
   const VectorPolynomial& vpoly, const AmanziGeometry::Point& x0,
   const std::vector<AmanziGeometry::Point>& tau);

// miscalleneous
DenseVector ExpandCoefficients(VectorPolynomial& vp);

// project gradient of the given polynomial on unit sphere using
// the Taylor expansion with k terms
VectorPolynomial GradientOnUnitSphere(const Polynomial& poly, int k);

}  // namespace WhetStone
}  // namespace Amanzi

#endif

