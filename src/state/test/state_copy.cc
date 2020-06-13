/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Markus Berndt
      Ethan Coon (coonet@ornl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_RCP.hpp"

#include "MeshFactory.hh"
#include "UnitTest++.h"

#include "Field.hh"
#include "Field_CompositeVector.hh"
#include "State.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;

struct test_field {
  Epetra_MpiComm* comm;
  Teuchos::RCP<Mesh> mesh;

  Teuchos::RCP<Field_CompositeVector> field;
  test_field()
  {
    comm = new Epetra_MpiComm(MPI_COMM_WORLD);
    MeshFactory mesh_fact(comm);
    mesh = mesh_fact(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);

    std::vector<Entity_kind> locations(2);
    locations[0] = CELL;
    locations[0] = FACE;

    std::vector<std::string> names(2);
    names[0] = "cell";
    names[1] = "face";

    std::vector<int> num_dofs(2, 1);

    Teuchos::RCP<CompositeVectorSpace> data_sp =
      Teuchos::rcp(new CompositeVectorSpace());
    data_sp.SetMesh(mesh).SetGhosted(false).SetComponents(
      names, locations, num_dofs);
    Teuchos::RCP<CompositeVector> data =
      Teuchos::rcp(new CompositeVector(*data_sp));
    field = Teuchos::rcp(new Field_CompositeVector("fieldname", "owner", data));
  }
  ~test_field() { delete comm; }
};

double
get_value(Teuchos::RCP<Field>& field)
{
  Teuchos::RCP<Field_CompositeVector> field_ptr =
    Teuchos::rcp_static_cast<Field_CompositeVector>(field);
  Teuchos::RCP<const CompositeVector> data = field_ptr->GetFieldData();
  return (*(*data->ViewComponent("cell", true))(0))[0];
};

double
get_value(Teuchos::RCP<Field_CompositeVector>& field)
{
  Teuchos::RCP<const CompositeVector> data = field->GetFieldData();
  return (*(*data->ViewComponent("cell", true))(0))[0];
};

double
get_value(const State& state)
{
  Teuchos::RCP<const CompositeVector> data = state.GetFieldData("fieldname");
  return (*(*data->ViewComponent("cell", true))(0))[0];
};

struct test_state {
  Epetra_MpiComm* comm;
  Teuchos::RCP<Mesh> mesh;

  Teuchos::RCP<State> state;

  test_state()
  {
    comm = new Epetra_MpiComm(MPI_COMM_WORLD);
    MeshFactory mesh_fact(comm);
    mesh = mesh_fact(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);
    state = Teuchos::rcp(new State());
    state->RegisterDomainMesh(mesh);

    std::vector<Entity_kind> locations(2);
    locations[0] = CELL;
    locations[1] = FACE;

    std::vector<std::string> names(2);
    names[0] = "cell";
    names[1] = "face";

    std::vector<int> num_dofs(2, 1);

    Teuchos::RCP<CompositeVectorSpace> vec_factory =
      state->RequireField("fieldname", "owner");
    vec_factory.SetMesh(state->GetMesh());
    vec_factory.SetComponents(names, locations, num_dofs);
    state->Setup();
    state->GetField("fieldname", "owner")->set_initialized();
    state->Initialize();
  }
  ~test_state() { delete comm; }
};

SUITE(COPY)
{
  // test the field copy constructor
  TEST_FIXTURE(test_field, FieldCopy)
  {
    field->GetFieldData()->putScalar(2.0);
    Teuchos::RCP<Field_CompositeVector> newfield =
      Teuchos::rcp(new Field_CompositeVector(*field));
    newfield->GetFieldData()->putScalar(3.0);
    CHECK_CLOSE(3.0, get_value(newfield), 0.00001);
    CHECK_CLOSE(2.0, get_value(field), 0.00001);
  }

  // test the field Clone constructor
  TEST_FIXTURE(test_field, FieldClone)
  {
    field->GetFieldData()->putScalar(2.0);
    Teuchos::RCP<Field> newfield = field->Clone();
    newfield->GetFieldData()->putScalar(3.0);
    CHECK_CLOSE(3.0, get_value(newfield), 0.00001);
    CHECK_CLOSE(2.0, get_value(field), 0.00001);
  }

  // test the state copy constructor
  TEST_FIXTURE(test_state, StateCopy)
  {
    state->GetFieldData("fieldname", "owner")->putScalar(2.0);
    State newstate(*state);

    newstate.GetFieldData("fieldname", "owner")->putScalar(3.0);
    CHECK_CLOSE(3.0, get_value(newstate), 0.00001);
    CHECK_CLOSE(2.0, get_value(*state), 0.00001);
  }

  // test the state assignment operator
  TEST_FIXTURE(test_state, StateAssignment)
  {
    state->GetFieldData("fieldname", "owner")->putScalar(2.0);

    // copy construct the new state to get the same structure
    Teuchos::RCP<State> newstate = Teuchos::rcp(new State(*state));

    // reset the value to test the operator=
    newstate->GetFieldData("fieldname", "owner")->putScalar(0.0);
    CHECK_CLOSE(2.0, get_value(*state), 0.00001);
    CHECK_CLOSE(0.0, get_value(*newstate), 0.00001);

    // test operator=
    *newstate = *state;
    CHECK_CLOSE(get_value(*newstate), get_value(*state), 0.00001);
    CHECK_CLOSE(2.0, get_value(*newstate), 0.00001);

    // test copies are deep
    newstate->GetFieldData("fieldname", "owner")->putScalar(3.0);
    CHECK_CLOSE(3.0, get_value(*newstate), 0.00001);
    CHECK_CLOSE(2.0, get_value(*state), 0.00001);
  }
}
