/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Process Kernels

  A collection of helper functions.
*/

#ifndef AMANZI_PK_STATE_HELPERS_HH_
#define AMANZI_PK_STATE_HELPERS_HH_

#include <string>

#include "Teuchos_ParameterList.hpp"

#include "CompositeVector.hh"
#include "Key.hh"
#include "State.hh"
#include "Tag.hh"

namespace Amanzi {

inline Teuchos::ParameterList
RequireFieldForEvaluator(State& S, const Key& key)
{
  Key domain = Keys::getDomain(key);
  S.Require<CompositeVector, CompositeVectorSpace>(key, Tags::DEFAULT, key)
    .SetMesh(S.GetMesh(domain))
    ->SetGhosted(true)
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  Teuchos::ParameterList elist(key);
  elist.set<std::string>("my key", key).set<std::string>("tag", Tags::DEFAULT.get());

  return elist;
}

} // namespace Amanzi

#endif
