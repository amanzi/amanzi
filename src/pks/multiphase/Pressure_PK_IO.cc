/*
  MultiPhase

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Quan Bui (mquanbui@math.umd.edu)

  This routine implements various I/O components of Pressure_PK which
  are common to Pressure_PK and Saturation_PK.
  Copied from Flow_IO.cc in Flow_PK and adapt it to multiphase.
*/

#include <vector>
#include "errors.hh"

// Amanzi
#include "FunctionTabular.hh"
#include "Mesh.hh"

// Multipahase
#include "Multiphase_BC_Factory.hh"
#include "Pressure_PK.hh"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* Routine processes parameter list. It needs to be called only once
* on each processor.                                                     
****************************************************************** */
void Pressure_PK::ProcessParameterList(Teuchos::ParameterList& plist)
{
  // Process main one-line options (not sublists)
  //atm_pressure_ = plist.get<double>("atmospheric pressure", FLOW_PRESSURE_ATMOSPHERIC);
 
  // Create the BC objects.
  if (!plist.isSublist("boundary conditions")) {
    Errors::Message msg;
    msg << "Pressure_PK: problem does not have <boundary conditions> list\n";
    Exceptions::amanzi_throw(msg);
  }
  Teuchos::RCP<Teuchos::ParameterList>
        bc_list = Teuchos::rcp(new Teuchos::ParameterList(plist.sublist("boundary conditions", true)));
  MultiphaseBCFactory bc_factory(mesh_, bc_list);

  bc_pressure = bc_factory.CreatePressure(bc_submodel);
  bc_flux = bc_factory.CreateMassFlux(bc_submodel, 0);
  
  // Create the source object if any
  /*
  if (!plist.isSublist("source terms")) {
    Errors::Message msg;
    msg << "Pressure_PK: problem does not have <source terms> list\n";
    Exceptions::amanzi_throw(msg);
  } else {
    std::string distribution_method_name = plist.get<std::string>("source and sink distribution method", "none");
    ProcessStringSourceDistribution(distribution_method_name, &src_sink_distribution_); 

    Teuchos::RCP<Teuchos::ParameterList> src_list = Teuchos::rcpFromRef(plist.sublist("source terms", true));
    Flow::FlowSourceFactory src_factory(mesh_, src_list);
    src_sink_ = src_factory.createSource();
    src_sink_distribution_ = src_sink_->CollectActionsList();
  }
  */
}


/* ****************************************************************
* Process string for the linear solver.
**************************************************************** */
void Pressure_PK::ProcessStringSourceDistribution(const std::string name, int* method)
{
  if (name != "none") {
    Errors::Message msg;
    msg << "\nPressure_PK: \"source and sink distribution method\" is obsolete.\n"
        << "         see desription of sublist \"source terms\" in the native spec.\n";
    Exceptions::amanzi_throw(msg);
  }
}

}  // namespace Flow
}  // namespace Amanzi

