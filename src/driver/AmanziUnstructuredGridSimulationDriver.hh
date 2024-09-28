/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Markus Brendt (brendt@lanl.gov)
*/

/*
  Simulator

*/

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_TimeMonitor.hpp"

#include "GeometricModel.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "ObservationData.hh"
#include "Simulator.hh"


struct AmanziUnstructuredGridSimulationDriver
  : public Amanzi::Simulator,
    public Teuchos::VerboseObject<AmanziUnstructuredGridSimulationDriver> {
 public:
  // constructor for native XML
  explicit AmanziUnstructuredGridSimulationDriver(const std::string& xmlInFileName);

  // constructor for v2 XML
  AmanziUnstructuredGridSimulationDriver(const std::string& xmlInFileName,
                                         xercesc::DOMDocument* input,
                                         const std::string& output_prefix);

  virtual Amanzi::Simulator::ReturnType
  Run(const Amanzi::Comm_ptr_type& comm, Amanzi::ObservationData& observations_data);

  virtual void
  Summarize() { if (getVerbLevel() > Teuchos::VERB_LOW) Teuchos::TimeMonitor::summarize(); }

  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> InitGeometricModel();
  int InitMesh(Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel>& gm,
               Teuchos::RCP<Amanzi::AmanziMesh::Mesh>& mesh,
               std::string& domain,
               Teuchos::RCP<Amanzi::AmanziMesh::Mesh>& submesh);

  // access
  void set_comm(Amanzi::Comm_ptr_type comm) { comm_ = comm; }
  Teuchos::RCP<Teuchos::ParameterList> get_plist() { return plist_; }

 private:
  // Read our parameter list.
  void ReadParameterList();
  Teuchos::RCP<Teuchos::ParameterList> plist_;
  Amanzi::Comm_ptr_type comm_;
};
