/*
  This is the Transport component of Amanzi.

  License: BSD
  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_DISPERSION_MATRIX_FACTORY_HH_
#define AMANZI_DISPERSION_MATRIX_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Dispersion.hh"

namespace Amanzi {
namespace AmanziTransport {

class DispersionMatrixFactory {
 public:
  DispersionMatrixFactory() {};
  ~DispersionMatrixFactory() {};

  Teuchos::RCP<Dispersion> Create(
     const string& matrix_name, std::vector<Teuchos::RCP<DispersionModel> >* specs,
     Teuchos::RCP<const AmanziMesh::Mesh> mesh, Teuchos::RCP<State> S);
};

}  // namespace AmanziTransport
}  // namespace Amanzi

#endif
