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
#include "Dispersion_MFD.hh"
#include "Dispersion_NLFV.hh"

namespace Amanzi {
namespace AmanziTransport {

/* ******************************************************************
 * Initialization of the dispersion matrix
 ****************************************************************** */
Teuchos::RCP<Dispersion> DispersionMatrixFactory::Create(
    const std::string& matrix_name, std::vector<Teuchos::RCP<DispersionModel> >* specs,
    Teuchos::RCP<const AmanziMesh::Mesh> mesh, Teuchos::RCP<State> S)
{
  if (matrix_name == "tpfa") {
    Teuchos::RCP<Dispersion_TPFA> matrix = Teuchos::rcp(new Dispersion_TPFA(specs, mesh, S));
    matrix->Init();
    return matrix;
  } else if (matrix_name == "mfd") {
    Teuchos::RCP<Dispersion_MFD> matrix = Teuchos::rcp(new Dispersion_MFD(specs, mesh, S));
    matrix->Init();
    return matrix;
  } else if (matrix_name == "nlfv") {
    Teuchos::RCP<Dispersion_NLFV> matrix = Teuchos::rcp(new Dispersion_NLFV(specs, mesh, S));
    matrix->Init();
    matrix->InitNLFV();  // additional initialization
    return matrix;
  } else {
    std::stringstream msgstream;
    msgstream << "DispersionMatrixFactory: Unknown matrix name " << matrix_name;
    Errors::Message message(msgstream.str());
    Exceptions::amanzi_throw(message);
  }
  return Teuchos::null;  // This line cannot be reached.
}

}  // namespace AmanziTransport
}  // namespace Amanzi
