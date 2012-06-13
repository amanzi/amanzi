/*
This is the transport component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided Reconstruction.cppin the top-level COPYRIGHT file.

Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <set>
#include <string>
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
void Transport_PK::ProcessParameterList()
{
  Teuchos::ParameterList transport_list;
  transport_list = parameter_list;

  // create verbosity list if it does not exist
  if (! transport_list.isSublist("VerboseObject")) {
    Teuchos::ParameterList verbosity_list;
    verbosity_list.set<std::string>("Verbosity Level", "none");
    transport_list.set("VerboseObject", verbosity_list);
  }

  // extract verbosity level
  Teuchos::ParameterList verbosity_list = transport_list.get<Teuchos::ParameterList>("VerboseObject");
  std::string verbosity_name = verbosity_list.get<std::string>("Verbosity Level");
  ProcessStringVerbosity(verbosity_name, &verbosity);

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

  string advection_limiter_name = transport_list.get<string>("advection limiter");
  ProcessStringAdvectionLimiter(advection_limiter_name, &advection_limiter);

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
  bcs.clear();
  bcs_tcc_index.clear();

  for (Teuchos::ParameterList::ConstIterator it = BCs_list.begin(); it != BCs_list.end(); ++it) {
    if (BCs_list.isSublist(it->first)) {
      Teuchos::ParameterList BC_list = BCs_list.sublist(it->first); 
      
      bool flag_BCX = false;
      for (int i = 0; i < number_components; i++) {
	char tcc_char_name[20];
	
	sprintf(tcc_char_name, "Component %d", i);
	string tcc_name(tcc_char_name);
	string tcc_name_alt(TS->get_component_name(i));
		
	if (BC_list.isParameter(tcc_name) || BC_list.isParameter(tcc_name_alt)) {
	  flag_BCX = true;
	  std::vector<std::string> regions, functions;
	  std::vector<double> times, values;
	  
	  regions = BC_list.get<Teuchos::Array<std::string> >("Regions").toVector();
	  times = BC_list.get<Teuchos::Array<double> >("Times").toVector();
	  if (BC_list.isParameter(tcc_name)) {
	    values = BC_list.get<Teuchos::Array<double> >(tcc_name).toVector();
	  } else if ( BC_list.isParameter(tcc_name_alt)) {
	    values = BC_list.get<Teuchos::Array<double> >(tcc_name_alt).toVector();
	  }
	  functions = BC_list.get<Teuchos::Array<std::string> >("Time Functions").toVector();
	  
	  int nfunctions = functions.size();  // convert strings to forms
	  std::vector<TabularFunction::Form> forms(functions.size());
	  for (int k = 0; k < nfunctions; k++) {
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
      if (! flag_BCX) {
	Errors::Message msg;
	msg << "Sublist BC X was not found.\n";
	Exceptions::amanzi_throw(msg);
      }
    }
  }
}

/* ****************************************************************
* Process string for the discretization method.
**************************************************************** */
void Transport_PK::ProcessStringVerbosity(const std::string name, int* verbosity)
{
  Errors::Message msg;
  if (name == "none") {
    *verbosity = TRANSPORT_VERBOSITY_NONE;
  } else if (name == "low") {
    *verbosity = TRANSPORT_VERBOSITY_LOW;
  } else if (name == "medium") {
    *verbosity = TRANSPORT_VERBOSITY_MEDIUM;
  } else if (name == "high") {
    *verbosity = TRANSPORT_VERBOSITY_HIGH;
  } else if (name == "extreme") {
    *verbosity = TRANSPORT_VERBOSITY_EXTREME;
  } else {
    msg << "Transport PK: unknown verbosity level.\n";
    Exceptions::amanzi_throw(msg);
  }
}


/* ****************************************************************
* Process time integration sublist.
**************************************************************** */
void Transport_PK::ProcessStringAdvectionLimiter(const std::string name, int* method)
{
  Errors::Message msg;
  if (name == "BarthJespersen") {
    advection_limiter = TRANSPORT_LIMITER_BARTH_JESPERSEN;
  } else if (name == "Tensorial") {
    advection_limiter = TRANSPORT_LIMITER_TENSORIAL;
  } else if (name == "Kuzmin") {
    advection_limiter = TRANSPORT_LIMITER_KUZMIN;
  } else {
    msg << "Transport PK: unknown advection limiter (BarthJespersen, Tensorial, Kuzmin).\n";
    Exceptions::amanzi_throw(msg);
  }
}


/* ************************************************************* */
/* Printing information about Transport status                   */
/* ************************************************************* */
void Transport_PK::PrintStatistics() const
{
  if (!MyPID && verbosity > TRANSPORT_VERBOSITY_NONE) {
    cout << "Transport PK: CFL = " << cfl_ << endl;
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
void Transport_PK::WriteGMVfile(Teuchos::RCP<Transport_State> TS) const
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

