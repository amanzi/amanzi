/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

  This smoothness indicator compares a reconstructed polynomial
  with that in neighbooring cells.
*/

#ifndef AMANZI_SMOOTHNESS_INDICATOR_SHU_HH_
#define AMANZI_SMOOTHNESS_INDICATOR_SHU_HH_

#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Mesh.hh"

#include "Reconstruction.hh"
#include "SmoothnessIndicator.hh"

namespace Amanzi {
namespace Operators {

class SmoothnessIndicatorShu : public SmoothnessIndicator {
 public:
  SmoothnessIndicatorShu(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh)
    : SmoothnessIndicator(mesh){};
  ~SmoothnessIndicatorShu(){};

  virtual void Init(Teuchos::ParameterList& plist);

  virtual void Compute(const Teuchos::RCP<Reconstruction>& lifting);
};

} // namespace Operators
} // namespace Amanzi

#endif
