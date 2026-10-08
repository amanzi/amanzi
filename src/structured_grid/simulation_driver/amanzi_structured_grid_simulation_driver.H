/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#ifndef _amanzi_structured_grid_simulation_driver_H
#define _amanzi_structured_grid_simulation_driver_H

#include "BoxLib.H"
#include "ParallelDescriptor.H"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_TimeMonitor.hpp"

#include "ObservationData.hh"
#include "Simulator.hh"

XERCES_CPP_NAMESPACE_USE

struct AmanziStructuredGridSimulationDriver
  : public Amanzi::Simulator,
    public Teuchos::VerboseObject<AmanziStructuredGridSimulationDriver> {

 public:
  AmanziStructuredGridSimulationDriver(const std::string& input_file, DOMDocument* input,
                                       const std::string& ouput_prefix);
  ~AmanziStructuredGridSimulationDriver();

  virtual ReturnType
  Run(const Amanzi::Comm_ptr_type& comm, Amanzi::ObservationData& output_observations);

  virtual void
  Summarize() { if (getVerbLevel() > Teuchos::VERB_LOW) Teuchos::TimeMonitor::summarize(); }
};
#endif
