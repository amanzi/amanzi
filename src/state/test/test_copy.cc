#include "UnitTest++.h"

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_RCP.hpp"

#include "Epetra_MpiComm.h"
#include "Epetra_MultiVector.h"

#include "MeshFactory.hh"
#include "Mesh_STK.hh"

#include "Field.hh"
#include "Field_CV.hh"
#include "State.hh"

using namespace Amanzi;

SUITE(COPY) {
  struct test_field {
    Epetra_MpiComm *comm;
    Teuchos::RCP<AmanziMesh::Mesh> mesh;

    Teuchos::RCP<Field_CV> field;
    test_field() {
      comm = new Epetra_MpiComm(MPI_COMM_WORLD);
      AmanziMesh::MeshFactory mesh_fact(*comm);
      mesh = mesh_fact(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);

      std::vector<AmanziMesh::Entity_kind> locations(2);
      locations[0] = AmanziMesh::CELL;
      locations[0] = AmanziMesh::FACE;

      std::vector<std::string> names(2);
      names[0] = "cell";
      names[1] = "face";

      Teuchos::RCP<CompositeVector> data =
        Teuchos::rcp(new CompositeVector(mesh, names, locations));
      field = Teuchos::rcp(new Field_CV("test_fieldname", "test_owner", data));
    }
    ~test_field() {
      delete comm;
    }
  };

  double get_value(Teuchos::RCP<Field>& field) {
    Teuchos::RCP<Field_CV> field_ptr = Teuchos::rcp_static_cast<Field_CV>(field);
    Teuchos::RCP<const CompositeVector> data = field_ptr->GetFieldData();
    return (*(*data->ViewComponent("cell", true))(0))[0];
  };

  double get_value(Teuchos::RCP<Field_CV>& field) {
    Teuchos::RCP<const CompositeVector> data = field->GetFieldData();
    return (*(*data->ViewComponent("cell", true))(0))[0];
  };

  double get_value(const State& state) {
    Teuchos::RCP<const CompositeVector> data = state.GetFieldData("test_fieldname");
    return (*(*data->ViewComponent("cell", true))(0))[0];
  };

  struct test_state {
    Epetra_MpiComm *comm;
    Teuchos::RCP<AmanziMesh::Mesh> mesh;

    Teuchos::RCP<State> state;

    test_state() {
      comm = new Epetra_MpiComm(MPI_COMM_WORLD);
      AmanziMesh::MeshFactory mesh_fact(*comm);
      mesh = mesh_fact(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);
      state = Teuchos::rcp(new State(mesh));

      std::vector<AmanziMesh::Entity_kind> locations(2);
      locations[0] = AmanziMesh::CELL;
      locations[0] = AmanziMesh::FACE;

      std::vector<std::string> names(2);
      names[0] = "cell";
      names[1] = "face";

      state->RequireField("test_fieldname", "test_owner", names, locations);
    }
    ~test_state() {
      delete comm;
    }
  };


  // test the field copy constructor
  TEST_FIXTURE(test_field, FieldCopy) {
    field->GetFieldData("test_owner")->PutScalar(2.0);
    Teuchos::RCP<Field_CV> newfield = Teuchos::rcp(new Field_CV(*field));
    CHECK_CLOSE(get_value(newfield), get_value(field), 0.00001);
    CHECK_CLOSE(get_value(newfield), 2.0, 0.00001);
    newfield->GetFieldData("test_owner")->PutScalar(3.0);
    CHECK_CLOSE(get_value(newfield), 3.0, 0.00001);
    CHECK_CLOSE(get_value(field), 2.0, 0.00001);
  }

  // test the field Clone constructor
  TEST_FIXTURE(test_field, FieldClone) {
    field->GetFieldData("test_owner")->PutScalar(2.0);
    Teuchos::RCP<Field> newfield = field->Clone();
    CHECK_CLOSE(get_value(newfield), get_value(field), 0.00001);
    CHECK_CLOSE(get_value(newfield), 2.0, 0.00001);
    newfield->GetFieldData("test_owner")->PutScalar(3.0);
    CHECK_CLOSE(get_value(newfield), 3.0, 0.00001);
    CHECK_CLOSE(get_value(field), 2.0, 0.00001);
  }

  // // test the field assignment operator
  // TEST_FIXTURE(test_field, FieldAssignment) {
  //   field->GetFieldData("test_owner")->PutScalar(2.0);

  //   // copy construct the new field to get the same structure
  //   Teuchos::RCP<Field_CV> newfield = Teuchos::rcp(new Field_CV(*field));
  //   *newfield = *field;
  //   CHECK_CLOSE(get_value(*newfield), get_value(*field), 0.00001);
  //   CHECK_CLOSE(get_value(*newfield), 2.0, 0.00001);

  //   // reset and check the new value
  //   newfield->GetFieldData("test_owner")->PutScalar(3.0);
  //   CHECK_CLOSE(get_value(*newfield), 3.0, 0.00001);
  //   CHECK_CLOSE(get_value(*field), 2.0, 0.00001);

  //   // test operator=
  //   *newfield = *field;
  //   CHECK_CLOSE(get_value(*newfield), get_value(*field), 0.00001);
  //   CHECK_CLOSE(get_value(*newfield), 2.0, 0.00001);

  //   // ensure operator= did a deep copy, not a pointer copy
  //   newfield->GetFieldData("test_owner")->PutScalar(3.0);
  //   CHECK_CLOSE(get_value(*newfield), 3.0, 0.00001);
  //   CHECK_CLOSE(get_value(*field), 2.0, 0.00001);
  // }

  // test the state copy constructor
  TEST_FIXTURE(test_state, StateCopy) {
    state->GetFieldData("test_fieldname", "test_owner")->PutScalar(2.0);
    State newstate(*state);
    CHECK_CLOSE(get_value(newstate), get_value(*state), 0.00001);
    CHECK_CLOSE(get_value(newstate), 2.0, 0.00001);
    newstate.GetFieldData("test_fieldname", "test_owner")->PutScalar(3.0);
    CHECK_CLOSE(get_value(newstate), 3.0, 0.00001);
    CHECK_CLOSE(get_value(*state), 2.0, 0.00001);
  }

  // test the state assignment operator
  TEST_FIXTURE(test_state, StateAssignment) {
    state->GetFieldData("test_fieldname", "test_owner")->PutScalar(2.0);

    // copy construct the new state to get the same structure
    Teuchos::RCP<State> newstate = Teuchos::rcp(new State(*state));

    // reset the value to test the operator=
    newstate->GetFieldData("test_fieldname", "test_owner")->PutScalar(0.0);
    CHECK_CLOSE(get_value(*state), 2.0, 0.00001);
    CHECK_CLOSE(get_value(*newstate), 0.0, 0.00001);

    // test operator=
    *newstate = *state;
    CHECK_CLOSE(get_value(*newstate), get_value(*state), 0.00001);
    CHECK_CLOSE(get_value(*newstate), 2.0, 0.00001);

    // test copies are deep
    newstate->GetFieldData("test_fieldname", "test_owner")->PutScalar(3.0);
    CHECK_CLOSE(get_value(*newstate), 3.0, 0.00001);
    CHECK_CLOSE(get_value(*state), 2.0, 0.00001);
  }
}
