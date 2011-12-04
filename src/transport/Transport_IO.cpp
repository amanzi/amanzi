/*
This is the transport component of the Amanzi code. 
License: BSD
Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <vector>

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "errors.hh"
#include "tabular-function.hh"

#include "Transport_PK.hpp"

namespace Amanzi {
namespace AmanziTransport {

/* ******************************************************************
* Routine processes parameter list. It needs to be called only once
* on each processor.                                                     
****************************************************************** */
void Transport_PK::process_parameter_list()
{
  Teuchos::RCP<AmanziMesh::Mesh> mesh = TS->get_mesh_maps();

  Teuchos::ParameterList transport_list;
  transport_list = parameter_list.get<Teuchos::ParameterList>("Transport");

  // global transport parameters
  cfl = transport_list.get<double>("CFL", 1.0);

  spatial_disc_order = transport_list.get<int>("spatial discretization order", 1);
  if (spatial_disc_order < 1 || spatial_disc_order > 2) spatial_disc_order = 1;
  temporal_disc_order = transport_list.get<int>("temporal discretization order", 1);
  if (temporal_disc_order < 1 || temporal_disc_order > 2) temporal_disc_order = 1;

  string dispersivity_name = transport_list.get<string>("dispersivity model", "none");
  if (dispersivity_name == "isotropic") { 
    dispersivity_model = TRANSPORT_DISPERSIVITY_MODEL_ISOTROPIC;
  }
  else if (dispersivity_name == "Bear") {
    dispersivity_model = TRANSPORT_DISPERSIVITY_MODEL_BEAR;
  }
  else if (dispersivity_name == "Lichtner") {
    dispersivity_model = TRANSPORT_DISPERSIVITY_MODEL_LICHTNER;
  }

  dispersivity_longitudinal = transport_list.get<double>("dispersivity longitudinal", 0.0);
  dispersivity_transverse = transport_list.get<double>("dispersivity transverse", 0.0);

  string advection_limiter_name = transport_list.get<string>("advection limiter", "Tensorial");
  if (advection_limiter_name == "BarthJespersen") {
    advection_limiter = TRANSPORT_LIMITER_BARTH_JESPERSEN;
  }
  else if (advection_limiter_name == "Tensorial") {
    advection_limiter = TRANSPORT_LIMITER_TENSORIAL;
  }

  flow_mode = TRANSPORT_FLOW_STEADYSTATE;
  string flow_mode_name = transport_list.get<string>("flow mode", "steady-state");
  if (flow_mode_name == "transient") {
    flow_mode = TRANSPORT_FLOW_TRANSIENT;
  }  
   
  // control parameter
  verbosity_level = transport_list.get<int>("verbosity level", 0);
  internal_tests = transport_list.get<string>("enable internal tests", "no") == "yes";
  tests_tolerance = transport_list.get<double>("internal tests tolerance", TRANSPORT_CONCENTRATION_OVERSHOOT);
  dT_debug = transport_list.get<double>("maximal time step", TRANSPORT_LARGE_TIME_STEP);
 
  // extract list of lists of boundary conditions
  Teuchos::ParameterList BCs_list;
  BCs_list = transport_list.get<Teuchos::ParameterList>("Transport BCs");
 
  // populate the list of boundary influx functions
  int nBCs = BCs_list.get<int>("number of BCs");
  bcs.clear();
  bcs_tcc_index.clear();

  for (int n=0; n<nBCs; n++) {
    char bc_char_name[10];
    
    sprintf(bc_char_name, "BC %d", n);
    string bc_name(bc_char_name);

    if (!BCs_list.isSublist(bc_name)) {
      Errors::Message msg;
      msg << "Boundary condition with name " << bc_char_name << " does not exist" << "\n";
      Exceptions::amanzi_throw(msg);
    }
    Teuchos::ParameterList BC_list = BCs_list.sublist(bc_name);  // A single sublist. 
 
    bool flag_BCX = false;
    for (int i=0; i<number_components; i++) {
      char tcc_char_name[20];

      sprintf(tcc_char_name, "Component %d", i);
      string tcc_name(tcc_char_name);

      if (BC_list.isParameter(tcc_name)) {
        flag_BCX = true; 
        std::vector<std::string> regions, functions;
        std::vector<double> times, values;

        regions = BC_list.get<Teuchos::Array<std::string> >("Regions").toVector();
        times = BC_list.get<Teuchos::Array<double> >("Times").toVector();
        values = BC_list.get<Teuchos::Array<double> >(tcc_name).toVector();
        functions = BC_list.get<Teuchos::Array<std::string> >("Time functions").toVector();

        int nfunctions = functions.size();  // convert strings to forms
        std::vector<TabularFunction::Form> forms(functions.size());
        for (int k=0; k<nfunctions; k++) {
          forms[k] = (functions[k] == "Constant") ? TabularFunction::CONSTANT : TabularFunction::LINEAR;
        }

        Teuchos::RCP<Function> f;
        f = Teuchos::rcp(new TabularFunction(times, values, forms));

        BoundaryFunction* bnd_fun = new BoundaryFunction(mesh_);
        bnd_fun->Define(regions, f);
        bcs.push_back(bnd_fun);
        bcs_tcc_index.push_back(i);
        break;
      }
    }
    if (!flag_BCX) {
      Errors::Message msg;
      msg << "Sublist BC X was not found.\n";
      Exceptions::amanzi_throw(msg);
    }
  }
}


/* ************************************************************* */
/* Printing information about Transport status                   */
/* ************************************************************* */
void Transport_PK::print_statistics() const
{
  if (!MyPID && verbosity_level > 0) {
    cout << "Transport PK: CFL = " << cfl << endl;
    cout << "    Execution mode = " << (standalone_mode ? "standalone" : "MPC") << endl;
    cout << "    Total number of components = " << number_components << endl;
    cout << "    Verbosity level = " << verbosity_level << endl;
    cout << "    Spatial/temporal discretication orders = " << spatial_disc_order 
         << " " << temporal_disc_order << endl;
    cout << "    Enable internal tests = " << (internal_tests ? "yes" : "no")  << endl;
    cout << "    Advection limiter = " << (advection_limiter == TRANSPORT_LIMITER_TENSORIAL ? "Tensorial" : "BarthJespersen") << endl;
  }
}
 
}  // namespace AmanziTransport
}  // namespace Amanzi

