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
#include "gmv_mesh.hh"

#include "Transport_PK.hpp"

namespace Amanzi {
namespace AmanziTransport {

/* ******************************************************************
* Routine processes parameter list. It needs to be called only once
* on each processor.                                                     
****************************************************************** */
void Transport_PK::processParameterList()
{
  Teuchos::ParameterList transport_list;
  transport_list = parameter_list;

  // create verbosity list if it does not exist
  if (!transport_list.isSublist("VerboseObject")) {
    Teuchos::ParameterList verbosity_list;
    verbosity_list.set<std::string>("Verbosity Level", "none");
    transport_list.set("VerboseObject", verbosity_list);
  }

  // extract verbosity level
  Teuchos::ParameterList verbosity_list = transport_list.get<Teuchos::ParameterList>("VerboseObject");
  std::string verbosity_name = verbosity_list.get<std::string>("Verbosity Level");
  if (verbosity_name == "none") {
    verbosity = TRANSPORT_VERBOSITY_NONE;
  } else if (verbosity_name == "low") {
    verbosity = TRANSPORT_VERBOSITY_LOW;
  } else if (verbosity_name == "medium") {
    verbosity = TRANSPORT_VERBOSITY_MEDIUM;
  } else if (verbosity_name == "high") {
    verbosity = TRANSPORT_VERBOSITY_HIGH;
  } else if (verbosity_name == "extreme") {
    verbosity = TRANSPORT_VERBOSITY_EXTREME;
  }

  Teuchos::RCP<AmanziMesh::Mesh> mesh = TS->mesh();

  // global transport parameters
  cfl_ = transport_list.get<double>("CFL", 1.0);

  spatial_disc_order = transport_list.get<int>("spatial discretization order", 1);
  if (spatial_disc_order < 1 || spatial_disc_order > 2) spatial_disc_order = 1;
  temporal_disc_order = transport_list.get<int>("temporal discretization order", 1);
  if (temporal_disc_order < 1 || temporal_disc_order > 2) temporal_disc_order = 1;

  string dispersivity_name = transport_list.get<string>("dispersivity model", "isotropic");
  if (dispersivity_name == "isotropic") { 
    dispersivity_model = TRANSPORT_DISPERSIVITY_MODEL_ISOTROPIC;
  } else if (dispersivity_name == "Bear") {
    dispersivity_model = TRANSPORT_DISPERSIVITY_MODEL_BEAR;
  } else if (dispersivity_name == "Lichtner") {
    dispersivity_model = TRANSPORT_DISPERSIVITY_MODEL_LICHTNER;
  } else {
    Errors::Message msg;
    msg << "Dispersivity model is wrong (isotropic, Bear, Lichtner)." << "\n";
    Exceptions::amanzi_throw(msg);    
  }

  dispersivity_longitudinal = transport_list.get<double>("dispersivity longitudinal", 0.0);
  dispersivity_transverse = transport_list.get<double>("dispersivity transverse", 0.0);

  string advection_limiter_name = transport_list.get<string>("advection limiter", "Tensorial");
  if (advection_limiter_name == "BarthJespersen") {
    advection_limiter = TRANSPORT_LIMITER_BARTH_JESPERSEN;
  } else if (advection_limiter_name == "Tensorial") {
    advection_limiter = TRANSPORT_LIMITER_TENSORIAL;
  } else if (advection_limiter_name == "Kuzmin") {
    advection_limiter = TRANSPORT_LIMITER_KUZMIN;
  } else {
    Errors::Message msg;
    msg << "Advection limiter is wrong (BarthJespersen, Tensorial, Kuzmin)." << "\n";
    Exceptions::amanzi_throw(msg);    
  }

  flow_mode = TRANSPORT_FLOW_TRANSIENT;
  string flow_mode_name = transport_list.get<string>("flow mode", "transient");
  if (flow_mode_name == "steady-state") {
    flow_mode = TRANSPORT_FLOW_STEADYSTATE;
  } else if (flow_mode_name == "transient") {
    flow_mode = TRANSPORT_FLOW_TRANSIENT;
  } else {
    Errors::Message msg;
    msg << "Flow mode is wrong (steady-state, transient)." << "\n";
    Exceptions::amanzi_throw(msg);
  }  
   
  // control parameter
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
        functions = BC_list.get<Teuchos::Array<std::string> >("Time Functions").toVector();

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
  //print_statistics();
}


/* ************************************************************* */
/* Printing information about Transport status                   */
/* ************************************************************* */
void Transport_PK::printStatistics() const
{
  if (!MyPID && verbosity > TRANSPORT_VERBOSITY_NONE) {
    cout << "Transport PK: CFL = " << cfl_ << endl;
    cout << "    Execution mode = " << (standalone_mode ? "standalone" : "MPC") << endl;
    cout << "    Total number of components = " << number_components << endl;
    cout << "    Verbosity level = " << verbosity << endl;
    cout << "    Spatial/temporal discretication orders = " << spatial_disc_order 
         << " " << temporal_disc_order << endl;
    cout << "    Enable internal tests = " << (internal_tests ? "yes" : "no")  << endl;
    cout << "    Advection limiter = " << (advection_limiter == TRANSPORT_LIMITER_TENSORIAL ? "Tensorial" : "BarthJespersen or Kuzmin(experimental)") << endl;
  }
}


/* ****************************************************************
* DEBUG: creating GMV file 
**************************************************************** */
void Transport_PK::writeGMVfile(Teuchos::RCP<Transport_State> TS) const
{
  Teuchos::RCP<AmanziMesh::Mesh> mesh = TS->mesh();
  Epetra_MultiVector& tcc = TS->ref_total_component_concentration();

  GMV::open_data_file(*mesh, (std::string)"transport.gmv");
  GMV::start_data();
  GMV::write_cell_data(tcc, 0, "component0");
  GMV::write_cell_data(TS->ref_water_saturation(), "saturation");
  GMV::close_data_file();
}
 
}  // namespace AmanziTransport
}  // namespace Amanzi

