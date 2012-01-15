/*
This is the flow component of the Amanzi code. 
License: BSD
Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <vector>

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "errors.hh"
#include "tabular-function.hh"

#include "Flow_PK.hpp"

namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* Routine processes parameter list. It needs to be called only once
* on each processor.                                                     
****************************************************************** */
void Flow_PK::processParameterList()
{
  Teuchos::ParameterList flow_list;
  flow_list = parameter_list.get<Teuchos::ParameterList>("Flow");
}


/* ************************************************************* */
/* Printing information about Transport status                   */
/* ************************************************************* */
void Flow_PK::print_statistics() const
{
  if (!MyPID && verbosity_level > 0) {
    cout << "Flow PK:" << endl;
    cout << "    Execution mode = " << (standalone_mode ? "standalone" : "MPC") << endl;
    cout << "    Verbosity level = " << verbosity_level << endl;
    cout << "    Enable internal tests = " << (internal_tests ? "yes" : "no")  << endl;
  }
}
 
}  // namespace AmanziFlow
}  // namespace Amanzi

