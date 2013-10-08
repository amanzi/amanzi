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
    const std::string& matrix_name, const Teuchos::ParameterList& matrix_list)
{
  if (matrix_name == "tpfa") {
    Teuchos::RCP<Dispersion_TPFA> prec = Teuchos::rcp(new Dispersion_TPFA());
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
