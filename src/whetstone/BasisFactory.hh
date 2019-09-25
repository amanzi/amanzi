/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Factory of polynomial bases for dG methods.
*/

#ifndef AMANZI_DG_BASIS_FACTORY_HH_
#define AMANZI_DG_BASIS_FACTORY_HH_

#include <memory>

#include "Basis_Natural.hh"
#include "Basis_Normalized.hh"
#include "Basis_Orthonormalized.hh"
#include "Basis_Regularized.hh"

namespace Amanzi {
namespace WhetStone {

template<class MyMesh>
class BasisFactory {
 public:
  explicit BasisFactory() {};
  ~BasisFactory() {};

  // select numerical scheme using its name and order 
  std::shared_ptr<Basis<MyMesh> > Create(const std::string& name) {
    if (name == "regularized") {
      auto basis = std::make_shared<Basis_Regularized<MyMesh> >();
      return basis;
    }
    else if (name == "normalized") {
      auto basis = std::make_shared<Basis_Normalized<MyMesh> >();
      return basis;
    }
    else if (name == "orthonormalized") {
      auto basis = std::make_shared<Basis_Orthonormalized<MyMesh> >();
      return basis;
    }
    else if (name == "natural") {
      auto basis = std::make_shared<Basis_Natural<MyMesh> >();
      return basis;
    }
    return NULL;
  }
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

