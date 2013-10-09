/*
  This is the Transport component of Amanzi.

  License: BSD
  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "Teuchos_RCP.hpp"

#include "errors.hh"

#include "DispersionMatrixFactory.hh"
#include "Dispersion.hh"
#include "Dispersion_TPFA.hh"

namespace Amanzi {
namespace AmanziTransport {

/* ******************************************************************
 * Initialization of the dispersion matrix
 ****************************************************************** */
Teuchos::RCP<Dispersion> DispersionMatrixFactory::Create(
    const string& matrix_name, std::vector<Teuchos::RCP<DispersionModel> >* specs,
    Teuchos::RCP<const AmanziMesh::Mesh> mesh, Teuchos::RCP<Transport_State> TS)
{
  if (matrix_name == "tpfa") {
    Teuchos::RCP<Dispersion_TPFA> prec = Teuchos::rcp(new Dispersion_TPFA(specs, mesh, TS));
    return prec;
  } else {
    std::stringstream msgstream;
    msgstream << "DispersionMatrixFactory: Unknown matrix name " << matrix_name;
    Errors::Message message(msgstream.str());
    Exceptions::amanzi_throw(message);
  }
}

}  // namespace AmanziTransport
}  // namespace Amanzi
