#include <iostream>
#include "stdlib.h"
#include "math.h"
#include "UnitTest++.h"


#include <Epetra_MpiComm.h>
#include "Epetra_SerialComm.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"

#include "CycleDriver.hh"
#include "IO.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "mpc_pks_registration.hh"
#include "PK_Factory.hh"
#include "PK.hh"
#include "State.hh"
#include "Tag.hh"

#include "ats_flow_pks_registration.hh"
#include "ats_flow_relations_registration.hh"
#include "ats_relations_registration.hh"
#include "pks_transport_registration.hh"
#include "pks_chemistry_registration.hh"

// Temporarily
namespace Amanzi {
namespace Flow {

class ATS_Richards : public Richards {
 public:
  ATS_Richards(Teuchos::ParameterList& pk_tree,
               const Teuchos::RCP<Teuchos::ParameterList>& glist,
               const Teuchos::RCP<State>& S,
               const Teuchos::RCP<TreeVector>& soln)
    : PK(pk_tree, glist, S, soln),
      Richards(pk_tree, glist, S, soln) {};

  virtual void Setup() override {
    Richards::Setup();

    auto owner = S_->GetRecord("saturation_liquid", Tags::DEFAULT).owner();
    S_->GetRecordW("saturation_liquid", Tags::DEFAULT, owner).set_owner("saturation_liquid");
    S_->RequireEvaluator("saturation_liquid", Tags::DEFAULT); 

    S_->Require<double>("time", Tags::CURRENT, "time");
    S_->Require<double>("atmospheric_pressure", Tags::DEFAULT, "coordinator");
    S_->Require<Amanzi::AmanziGeometry::Point>("gravity", Tags::DEFAULT, "coordinator");
  }

  virtual void Initialize() override {
    Richards::Initialize();

    S_->GetRecordW("darcy_flux", Tags::DEFAULT, "state").set_initialized();
  }

  virtual bool AdvanceStep(double t_old, double t_new, bool reinit) override {
    S_->GetW<double>("time", Tags::CURRENT, "time") = t_old;
    S_->GetW<double>("time", Tags::NEXT, "time") = t_new;

    S_->GetW<CompositeVector>("water_content", Tags::CURRENT, "flow") =
      S_->Get<CompositeVector>("water_content", Tags::NEXT);

    S_->GetW<CompositeVector>("saturation_liquid", Tags::CURRENT, "flow") =
      S_->Get<CompositeVector>("saturation_liquid", Tags::NEXT);

    Richards::AdvanceStep(t_old, t_new, reinit);

    return false;
  }

  virtual void CommitStep(double t_old, double t_new, const Tag& tag) override {
    Richards::CommitStep(t_old, t_new, Tags::NEXT);

    S_->GetW<CompositeVector>("darcy_flux", Tags::DEFAULT, "state") =
      S_->Get<CompositeVector>("water_flux", Tags::NEXT);

    // reset time to beginning of time step as expected by Amanzi FIXME
    S_->GetW<double>("time", Tags::NEXT, "time") = t_old;
  }

  virtual void CalculateDiagnostics(const Tag& tag) override {
    Richards::CalculateDiagnostics(Tags::NEXT);
  }

 private:
  static RegisteredPKFactory<ATS_Richards> reg_;
};

RegisteredPKFactory<ATS_Richards> ATS_Richards::reg_("ats richards flow");

}  // namespace Flow
}  // namespace Amanzi


TEST(INTEROPERABILITY_FLOW_TRANSPORT) {

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;

  auto comm = Amanzi::getDefaultComm();
  
  // read the main parameter list
  std::string xmlInFileName = "test/interoperability_flow_transport.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlInFileName);
  Teuchos::ParameterList plist = xmlreader.getParameters();
  
  // For now create one geometric model from all the regions in the spec
  Teuchos::ParameterList region_list = plist.get<Teuchos::ParameterList>("regions");
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
      Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, region_list, *comm));
   
  // create mesh
  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  Teuchos::RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 100, 1, 1);
  AMANZI_ASSERT(!mesh.is_null());

  // create dummy observation data object
  double avg1;
  Amanzi::ObservationData obs_data;    
  Teuchos::RCP<Teuchos::ParameterList> glist = Teuchos::rcp(new Teuchos::ParameterList(plist));

  Teuchos::ParameterList state_plist = glist->sublist("state");
  Teuchos::RCP<Amanzi::State> S = Teuchos::rcp(new Amanzi::State(state_plist));
  S->RegisterMesh("domain", mesh);

  {
    Amanzi::CycleDriver cycle_driver(glist, S, comm, obs_data);
    try {
      cycle_driver.Go();
      S->Get<CompositeVector>("pressure").MeanValue(&avg1);
    } catch (std::exception& e) {
      std::cout << e.what() << std::endl;
    } catch (...) {
      CHECK(false);
    }
  }
}

