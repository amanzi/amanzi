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

#include "errors.hh"
#include "MultiFunction.hh"
#include "GMVMesh.hh"

#include "Mesh.hh"
#include "Transport_PK.hh"
#include "Transport_BC_Factory.hh"

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

  Teuchos::RCP<const AmanziMesh::Mesh> mesh = TS->mesh();

  // global transport parameters
  cfl_ = transport_list.get<double>("CFL", 1.0);

  spatial_disc_order = transport_list.get<int>("spatial discretization order", 1);
  if (spatial_disc_order < 1 || spatial_disc_order > 2) spatial_disc_order = 1;
  temporal_disc_order = transport_list.get<int>("temporal discretization order", 1);
  if (temporal_disc_order < 1 || temporal_disc_order > 2) temporal_disc_order = 1;

  string dispersivity_name = transport_list.get<string>("dispersivity model", "isotropic");
  ProcessStringDispersionModel(dispersivity_name, &dispersivity_model);

  dispersivity_longitudinal = transport_list.get<double>("dispersivity longitudinal", 0.0);
  dispersivity_transverse = transport_list.get<double>("dispersivity transverse", 0.0);

  string advection_limiter_name = transport_list.get<string>("advection limiter");
  ProcessStringAdvectionLimiter(advection_limiter_name, &advection_limiter);

  std::string flow_mode_name = transport_list.get<string>("flow mode", "transient");
  ProcessStringFlowMode(flow_mode_name, &flow_mode);

  // control parameter
  internal_tests = transport_list.get<string>("enable internal tests", "no") == "yes";
  tests_tolerance = transport_list.get<double>("internal tests tolerance", TRANSPORT_CONCENTRATION_OVERSHOOT);
  dT_debug = transport_list.get<double>("maximum time step", TRANSPORT_LARGE_TIME_STEP);

  // populate the list of boundary influx functions
  bcs.clear();
  bcs_tcc_index.clear();

  if (transport_list.isSublist("boundary conditions")) {  // New flexible format.
    std::vector<std::string> bcs_tcc_name;
    Teuchos::RCP<Teuchos::ParameterList>
       bcs_list = Teuchos::rcp(new Teuchos::ParameterList(transport_list.get<Teuchos::ParameterList>("boundary conditions")));
    TransportBCFactory bc_factory(mesh_, bcs_list);
    bc_factory.CreateConcentration(bcs, bcs_tcc_name);
    for (int i = 0; i < bcs_tcc_name.size(); i++) {
      bcs_tcc_index.push_back(TS->get_component_number(bcs_tcc_name[i]));
    }
  } else {
    Errors::Message msg;
    msg << "Transport PK: does not have boundary conditions.\n";
    Exceptions::amanzi_throw(msg);  
  }

  // Create the source object if any
  if (transport_list.isSublist("source terms")) {
    Teuchos::RCP<Teuchos::ParameterList> src_list = Teuchos::rcpFromRef(transport_list.sublist("source terms", true));
    TransportSourceFactory src_factory(mesh_, src_list);
    src_sink = src_factory.CreateSource();
    
    // revisit the code below (lipnikov@lanl.gov)
    src_sink_distribution = src_sink->CollectActionsList();
    if (src_sink_distribution & Functions::TransportActions::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_PERMEABILITY) {
      Errors::Message msg;
      msg << "Transport PK: support of permeability weighted source distribution is pending.\n";
      Exceptions::amanzi_throw(msg);  
    }
  } else {
    src_sink = NULL;
    src_sink_distribution = 0;
  }
}


/* ****************************************************************
* Process string for the discretization method.
**************************************************************** */
void Transport_PK::ProcessStringFlowMode(const std::string name, int* method)
{
  Errors::Message msg;
  *method = TRANSPORT_FLOW_TRANSIENT;
  if (name == "steady-state") {
    *method = TRANSPORT_FLOW_STEADYSTATE;
  } else if (name == "transient") {
    *method = TRANSPORT_FLOW_TRANSIENT;
  } else {
    msg << "Trasnport PK: flow mode is wrong (steady-state, transient)." << "\n";
    Exceptions::amanzi_throw(msg);
  }

}


/* ****************************************************************
* Process string for the discretization method.
**************************************************************** */
void Transport_PK::ProcessStringDispersionModel(const std::string name, int* method)
{
  Errors::Message msg;
  if (name == "isotropic") {
    *method = TRANSPORT_DISPERSIVITY_MODEL_ISOTROPIC;
  } else if (name == "Bear") {
    *method = TRANSPORT_DISPERSIVITY_MODEL_BEAR;
  } else if (name == "Lichtner") {
    *method = TRANSPORT_DISPERSIVITY_MODEL_LICHTNER;
  } else {
    Errors::Message msg;
    msg << "Dispersivity model is wrong (isotropic, Bear, Lichtner).\n";
    Exceptions::amanzi_throw(msg);
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
  Teuchos::RCP<const AmanziMesh::Mesh> mesh = TS->mesh();
  Epetra_MultiVector& tcc = TS->ref_total_component_concentration();

  GMV::open_data_file(*mesh, (std::string)"transport.gmv");
  GMV::start_data();
  GMV::write_cell_data(tcc, 0, "component0");
  GMV::write_cell_data(TS->ref_water_saturation(), "saturation");
  GMV::close_data_file();
}

}  // namespace AmanziTransport
}  // namespace Amanzi

