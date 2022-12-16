/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  WhetStone, Version 2.2
  Release name: naka-to.

  The base class for dG basis.

  Let each column of R represent new basis vector in the old basis.
  Then, transformation of the bilinear form matrix A from old to new
  basis is calculated as follows: R^T A R. Also, to change vector v
  from the new to old basis, we compute v_old = R v_new. Inverse of
  R defines the backward transformation. Finally, transformation of
  the linear form, represented by vector f, is given by R^T f.
*/

#ifndef AMANZI_DG_BASIS_HH_
#define AMANZI_DG_BASIS_HH_

#include <memory>

#include "Teuchos_RCP.hpp"

#include "MeshLight.hh"

#include "DenseMatrix.hh"
#include "DenseVector.hh"
#include "Polynomial.hh"

namespace Amanzi {
namespace WhetStone {

class Basis {
 public:
  Basis(){};
  virtual ~Basis(){};

  // initialization
  virtual void Init(const Teuchos::RCP<const AmanziMesh::MeshLight>& mesh,
                    int c,
                    int order,
                    Polynomial& integrals) = 0;

  // transformation of bilinear form
  virtual void BilinearFormNaturalToMy(DenseMatrix& A) const = 0;
  virtual void BilinearFormNaturalToMy(std::shared_ptr<Basis> bl,
                                       std::shared_ptr<Basis> br,
                                       DenseMatrix& A) const = 0;

  // transformation of a linear form
  virtual void LinearFormNaturalToMy(DenseVector& v) const = 0;

  // transformation of vector
  virtual void ChangeBasisMyToNatural(DenseVector& v) const = 0;
  virtual void ChangeBasisNaturalToMy(DenseVector& v) const = 0;

  // recover polynomial in the natural basis
  virtual Polynomial CalculatePolynomial(const Teuchos::RCP<const AmanziMesh::MeshLight>& mymesh,
                                         int c,
                                         int order,
                                         DenseVector& coefs) const = 0;

  // assess
  int id() const { return id_; };

 protected:
  int id_;
};

} // namespace WhetStone
} // namespace Amanzi

#endif
