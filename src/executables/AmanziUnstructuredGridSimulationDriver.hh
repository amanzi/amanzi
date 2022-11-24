/*
  Simulator

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Markus Brendt (brendt@lanl.gov)
*/

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "ObservationData.hh"
#include "Simulator.hh"


struct AmanziUnstructuredGridSimulationDriver
  : Amanzi::Simulator,
    public Teuchos::VerboseObject<AmanziUnstructuredGridSimulationDriver> {
 public:
  // constructor for native XML
  explicit AmanziUnstructuredGridSimulationDriver(const std::string& xmlInFileName);

  // constructor for v2 XML
  AmanziUnstructuredGridSimulationDriver(const std::string& xmlInFileName,
                                         xercesc::DOMDocument* input,
                                         const std::string& output_prefix);

  Amanzi::Simulator::ReturnType
  Run(const Amanzi::Comm_ptr_type& comm, Amanzi::ObservationData& observations_data) override;

 private:
  // Read our parameter list.
  void ReadParameterList();
  Teuchos::RCP<Teuchos::ParameterList> plist_;
};
