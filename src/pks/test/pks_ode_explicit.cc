/*
  State

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "Epetra_MpiComm.h"
#include "Epetra_Vector.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "UnitTest++.h"

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "State.hh"
#include "TreeVector.hh"

#include "PK.hh"
#include "PK_Default.hh"
#include "PK_MixinLeaf.hh"
#include "PK_MixinExplicit.hh"
#include "PK_Adaptors.hh"

using namespace Amanzi;

template<class Base_t>
class PK_ODE_Explicit : public Base_t {
 public:
  PK_ODE_Explicit(const Teuchos::RCP<Teuchos::ParameterList>& pk_tree,
                  const Teuchos::RCP<Teuchos::ParameterList>& global_plist,
                  const Teuchos::RCP<State>& S,
                  const Teuchos::RCP<TreeVector>& solution) :
      Base_t(pk_tree, global_plist, S, solution) {}

  void Setup() {
    Base_t::Setup();
    
    std::cout << "requiring: vector " << this->key_ << std::endl;
    this->S_->template Require<CompositeVector,CompositeVectorSpace>(this->key_)
        .SetComponent("cell",AmanziMesh::CELL,1)
        ->SetGhosted(false);
  }

  void Initialize() {
    Base_t::Initialize();
    this->S_->template GetW<CompositeVector>("primary", "", "primary").PutScalar(1.);
    this->S_->GetRecordW("primary", "primary").set_initialized();
  }
  
  void Functional(double t, const TreeVector& u, TreeVector& f) { f.PutScalar(1.); }
};




Teuchos::RCP<PK>
create(const Teuchos::RCP<State>& S) {
  Teuchos::RCP<Teuchos::ParameterList> pk_tree = Teuchos::rcp(new Teuchos::ParameterList("my pk"));
  Teuchos::RCP<Teuchos::ParameterList> global_list = Teuchos::rcp(new Teuchos::ParameterList("main"));
  auto sublist = Teuchos::sublist(Teuchos::sublist(global_list, "PKs"), "my pk");
  sublist->set<std::string>("domain name", "domain");
  sublist->set<std::string>("primary variable key", "primary");

  // intentionally leaks memory
  Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_WORLD);
  
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "Test: 2D transient Darcy, 2-layer model" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/flow_darcy_transient_2D.xml";

  // create mesh
  Teuchos::ParameterList regions_list;
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
      Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, regions_list, comm));

  Amanzi::AmanziMesh::FrameworkPreference pref;
  pref.clear();
  pref.push_back(Amanzi::AmanziMesh::MSTK);

  Amanzi::AmanziMesh::MeshFactory meshfactory(comm);
  meshfactory.preference(pref);
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh = meshfactory(0.0, -2.0, 1.0, 0.0, 18, 18, gm);
  S->RegisterDomainMesh(mesh);

  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
  return Teuchos::rcp(new PK_Adaptor<PK_ODE_Explicit<PK_MixinExplicit<PK_MixinLeaf<PK_Default> > > >(pk_tree, global_list, S, soln));
}


SUITE(PK_ODE_EXPLICIT) {

  TEST(Construct) {
    auto S = Teuchos::rcp(new State());
    auto pk = create(S);
    pk->Setup();
    S->Setup();
    pk->Initialize();
    S->Initialize();
  }


  TEST(Advance) {
    auto S = Teuchos::rcp(new State());
    auto pk = create(S);
    pk->Setup();
    S->Setup();
    pk->Initialize();
    S->Initialize();

    double dt = pk->get_dt();
    S->advance_time("next", dt);
    S->advance_cycle("next", 1);
    pk->AdvanceStep("", "next");

    CHECK_CLOSE(1.0, (*S->Get<CompositeVector>("primary", "").ViewComponent("cell",false))[0][0], 1.e-10);
    CHECK_CLOSE(2.0, (*S->Get<CompositeVector>("primary", "next").ViewComponent("cell",false))[0][0], 1.e-10);
    
  }

  TEST(Commit) {
    auto S = Teuchos::rcp(new State());
    auto pk = create(S);
    pk->Setup();
    S->Setup();
    pk->Initialize();
    S->Initialize();

    double dt = pk->get_dt();
    S->advance_time("next", dt);
    S->advance_cycle("next", 1);
    pk->AdvanceStep("", "next");

    CHECK_CLOSE(1.0, (*S->Get<CompositeVector>("primary", "").ViewComponent("cell",false))[0][0], 1.e-10);
    CHECK_CLOSE(2.0, (*S->Get<CompositeVector>("primary", "next").ViewComponent("cell",false))[0][0], 1.e-10);
  }

    
}

