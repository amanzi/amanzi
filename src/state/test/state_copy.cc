/*
  State

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt
           Ethan Coon (ecoon@lanl.gov)
*/

#include "Epetra_MpiComm.h"
#include "Epetra_MultiVector.h"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_RCP.hpp"

#include "MeshFactory.hh"
#include "UnitTest++.h"

#include "State.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;

struct test_state {
  Comm_ptr_type comm;
  Teuchos::RCP<Mesh> mesh;
  Teuchos::RCP<State> state;

  test_state() {
    comm = Amanzi::getDefaultComm();
    MeshFactory mesh_fact(comm);
    mesh = mesh_fact.create(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);
    state = Teuchos::rcp(new State());
    state->RegisterDomainMesh(mesh);

    std::vector<Entity_kind> locations(2);
    locations[0] = CELL;
    locations[1] = FACE;

    std::vector<std::string> names(2);
    names[0] = "cell";
    names[1] = "face";

    Tag tag1 = make_tag("tag1");
    Tag tag2 = make_tag("tag2");
    std::vector<int> num_dofs(2, 1);

    state->Require<CompositeVector, CompositeVectorSpace>("fieldname", tag1, "owner")
      .SetMesh(state->GetMesh())->SetComponents(names, locations, num_dofs);

    state->Require<CompositeVector, CompositeVectorSpace>("fieldname", tag2, "owner")
      .SetMesh(state->GetMesh())->SetComponents(names, locations, num_dofs);

    state->Setup();
    // state->GetRecordW("fieldname", "owner".set_initialized();
    // state->Initialize();
  }
  ~test_state() {};
};


SUITE(COPY) {
  TEST_FIXTURE(test_state, StateCopy) {
    Tag tag1 = make_tag("tag1");
    Tag tag2 = make_tag("tag2");
    state->GetW<CompositeVector>("fieldname", tag1, "owner").PutScalar(2.0);

    state->Set<CompositeVector>("fieldname", tag2, "owner",
                                state->Get<CompositeVector>("fieldname", tag1));

    auto& field = state->Get<CompositeVector>("fieldname", tag2);
    CHECK_CLOSE(2.0, (*field.ViewComponent("cell"))[0][0], 0.00001);
    std::cout << "State copy... passed\n";
  }
}
