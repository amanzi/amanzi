#include "UnitTest++.h"

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_RCP.hpp"

#include "Epetra_MpiComm.h"
#include "Epetra_MultiVector.h"

#include "MeshFactory.hh"
#include "Mesh_STK.hh"

#include "Field.hh"
#include "State.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;

int main (int argc, char *argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
  return UnitTest::RunAllTests();
}

struct test_field {
  Epetra_MpiComm *comm;
  Teuchos::RCP<Mesh> mesh;

  Teuchos::RCP<Field> field;
  test_field() {
    comm = new Epetra_MpiComm(MPI_COMM_WORLD);
    MeshFactory mesh_fact(*comm);
    mesh = mesh_fact(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);
    field = Teuchos::rcp(new Field("test_fieldname", FIELD_LOCATION_CELL,
            mesh, "test_owner"));
  }
  ~test_field() { delete comm; }
};

double get_value(Field& field) {
  Teuchos::RCP<Epetra_MultiVector> data = field.get_data("test_owner");
  return (*(*data)(0))[0];
};

double get_value(State& state) {
  Teuchos::RCP<Epetra_MultiVector> data = state.get_field("test_fieldname", "test_owner");
  return (*(*data)(0))[0];
};

struct test_state {
  Epetra_MpiComm *comm;
  Teuchos::RCP<Mesh> mesh;

  Teuchos::RCP<State> state;
  test_state() {
    comm = new Epetra_MpiComm(MPI_COMM_WORLD);
    MeshFactory mesh_fact(*comm);
    mesh = mesh_fact(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);
    state = Teuchos::rcp(new State(mesh));
    state->require_field("test_fieldname", FIELD_LOCATION_CELL, "test_owner");
  }
  ~test_state() { delete comm; }
};

SUITE(COPY) {
  // test the field copy constructor
  TEST_FIXTURE(test_field, FieldCopy) {
    field->set_data("test_owner", 2.0);
    Field newfield(*(field));
    CHECK_CLOSE(get_value(newfield), get_value(*field), 0.00001);
    CHECK_CLOSE(get_value(newfield), 2.0, 0.00001);
    newfield.set_data("test_owner", 3.0);
    CHECK_CLOSE(get_value(newfield), 3.0, 0.00001);
    CHECK_CLOSE(get_value(*field), 2.0, 0.00001);
  }

  // test the field assignment operator
  TEST_FIXTURE(test_field, FieldAssignment) {
    field->set_data("test_owner", 2.0);
    Teuchos::RCP<Field> newfield = Teuchos::rcp(new Field("new_fieldname",
            FIELD_LOCATION_CELL, mesh, "test_owner"));
    *newfield = *field;
    CHECK_CLOSE(get_value(*newfield), get_value(*field), 0.00001);
    CHECK_CLOSE(get_value(*newfield), 2.0, 0.00001);
    newfield->set_data("test_owner", 3.0);
    CHECK_CLOSE(get_value(*newfield), 3.0, 0.00001);
    CHECK_CLOSE(get_value(*field), 2.0, 0.00001);
  }

  // test the state copy constructor
  TEST_FIXTURE(test_state, StateCopy) {
    state->set_field("test_fieldname", "test_owner", 2.0);
    State newstate(*state);
    CHECK_CLOSE(get_value(newstate), get_value(*state), 0.00001);
    CHECK_CLOSE(get_value(newstate), 2.0, 0.00001);
    newstate.set_field("test_fieldname", "test_owner", 3.0);
    CHECK_CLOSE(get_value(newstate), 3.0, 0.00001);
    CHECK_CLOSE(get_value(*state), 2.0, 0.00001);
  }

  // test the state assignment operator
  TEST_FIXTURE(test_state, StateAssignment) {
    state->set_field("test_fieldname", "test_owner", 2.0);

    // copy construct the new state to get the same structure
    Teuchos::RCP<State> newstate = Teuchos::rcp(new State(*state));

    // reset the value to test the operator=
    newstate->set_field("test_fieldname", "test_owner", 0.0);

    // test operator=
    *newstate = *state;
    CHECK_CLOSE(get_value(*newstate), get_value(*state), 0.00001);
    CHECK_CLOSE(get_value(*newstate), 2.0, 0.00001);

    // test copies are deep
    newstate->set_field("test_fieldname", "test_owner", 3.0);
    CHECK_CLOSE(get_value(*newstate), 3.0, 0.00001);
    CHECK_CLOSE(get_value(*state), 2.0, 0.00001);
  }
}
