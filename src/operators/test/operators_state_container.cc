/*
 State

 Tests for state as a container of data
*/

// TPLs
#include "UnitTest++.h"

#include "IO.hh"
#include "MeshFactory.hh"
#include "Operator_DataHelpers.hh"
#include "State.hh"
#include "errors.hh"

#include "Op_Cell_Cell.hh"
#include "Op_Factory.hh"

TEST(STATE_VIRTUAL_DATA) {
  using namespace Amanzi;

  // create a mesh
  Comm_ptr_type comm = Amanzi::getDefaultComm();
  AmanziMesh::MeshFactory fac(comm);
  auto mesh = fac.create(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);

  // create a state
  State s;
  s.RegisterDomainMesh(mesh);

  // require some data
  auto& f = s.Require<Operators::Op, Operators::Op_Factory>("my_op", Tags::DEFAULT, "my_op_owner");
  f.set_mesh(mesh);
  f.set_name("cell");
  f.set_schema(Operators::Schema{Operators::OPERATOR_SCHEMA_BASE_CELL |
                                 Operators::OPERATOR_SCHEMA_DOFS_CELL});

  s.Setup();

  // existence
  CHECK(s.HasRecord("my_op"));
  CHECK_EQUAL(Operators::OPERATOR_SCHEMA_DOFS_CELL |
                  Operators::OPERATOR_SCHEMA_BASE_CELL,
              s.Get<Operators::Op>("my_op").schema_old());
}
