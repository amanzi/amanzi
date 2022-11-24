/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  The base class for various smoothness indicators.
*/

#ifndef AMANZI_SMOOTHNESS_INDICATOR_HH_
#define AMANZI_SMOOTHNESS_INDICATOR_HH_

#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Mesh.hh"

#include "Reconstruction.hh"

namespace Amanzi {
namespace Operators {

class SmoothnessIndicator {
 public:
  SmoothnessIndicator(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh) : mesh_(mesh){};
  virtual ~SmoothnessIndicator() = default;

  virtual void Init(Teuchos::ParameterList& plist) = 0;

  virtual void Compute(const Teuchos::RCP<Reconstruction>& lifting) = 0;

  std::vector<double> get_measure() const { return measure_; }

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  std::vector<double> measure_;
};

} // namespace Operators
} // namespace Amanzi

#endif
