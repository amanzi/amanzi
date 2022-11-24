/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Natural basis with basis with monomials of form (x-x0)^k. The
  transformation matrix R is identity. This class is introduced 
  for factory to work uniformly.
*/

#ifndef AMANZI_DG_BASIS_NATURAL_HH_
#define AMANZI_DG_BASIS_NATURAL_HH_

#include "Teuchos_RCP.hpp"

#include "MeshLight.hh"

#include "Basis.hh"
#include "DenseMatrix.hh"
#include "DenseVector.hh"
#include "Polynomial.hh"
#include "WhetStoneDefs.hh"

namespace Amanzi {
namespace WhetStone {

class Basis_Natural : public Basis {
 public:
  Basis_Natural() { id_ = TAYLOR_BASIS_NATURAL; }
  ~Basis_Natural(){};

  // initialization
  virtual void Init(const Teuchos::RCP<const AmanziMesh::MeshLight>& mesh,
                    int c,
                    int order,
                    Polynomial& integrals)
  {
    order_ = order;
  }

  // transformation of bilinear form
  virtual void BilinearFormNaturalToMy(DenseMatrix& A) const {};
  virtual void BilinearFormNaturalToMy(std::shared_ptr<Basis> bl,
                                       std::shared_ptr<Basis> br,
                                       DenseMatrix& A) const {};

  // transformation of linear form
  virtual void LinearFormNaturalToMy(DenseVector& v) const {};

  // transformation of vector
  virtual void ChangeBasisMyToNatural(DenseVector& v) const {};
  virtual void ChangeBasisNaturalToMy(DenseVector& v) const {};

  // recover polynomial in natural basis
  virtual Polynomial CalculatePolynomial(const Teuchos::RCP<const AmanziMesh::MeshLight>& mesh,
                                         int c,
                                         int order,
                                         DenseVector& coefs) const
  {
    int d = mesh->space_dimension();
    Polynomial poly(d, order, coefs);
    poly.set_origin(mesh->cell_centroid(c));
    return poly;
  }

  // assess
  int id() { return id_; };

 protected:
  int id_;

 private:
  int order_;
};

} // namespace WhetStone
} // namespace Amanzi

#endif
