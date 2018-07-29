/*
  WhetStone, version 2.1
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Virtual framework for L2 and H1 projectors. The projectors
  may require degrees of freedom on cell boundary (for low-order 
  and high-order schemes) and inside it (for high-order schemes).
  Degrees of freedom on faces are typically function moments or,
  equivalently, polynomial approximation of the function. Degrees
  of freedom in cell are high-order moments. 

  NOTE. For a scheme of order k, polynomial of order k on each 
  face, provides uniform interface for large variety of schemes.

  NOTE. For a scheme of order k, internal moment of order k-s are
  needed, where s=2 for most implemented schemes and s > 2 for 
  serendipity schemes. At the moment, we assume s=2 and return 
  moments that are derived quantities. 
*/

#ifndef AMANZI_WHETSTONE_PROJECTORS_HH_
#define AMANZI_WHETSTONE_PROJECTORS_HH_

#include "Teuchos_RCP.hpp"

namespace Amanzi {
namespace WhetStone {

class DenseVector;
class VectorPolynomial;

class Projectors { 
 public:
  enum class Type {
    L2,
    H1,
    LS   // least square
  };

 public:
  explicit Projectors() {};
  ~Projectors() {};

  // elliptic projector
  virtual void H1Cell(
      int c, const std::vector<VectorPolynomial>& vf,
      VectorPolynomial& moments, VectorPolynomial& uc) {
    Errors::Message msg("Elliptic projector is not supported for this scheme.");
    Exceptions::amanzi_throw(msg);
  }

  // L2 projector 
  virtual void L2Cell(
      int c, const std::vector<VectorPolynomial>& vf,
      VectorPolynomial& moments, VectorPolynomial& uc) {
    Errors::Message msg("L2 projector is not supported for this scheme.");
    Exceptions::amanzi_throw(msg);
  }
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

