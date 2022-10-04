/*
 State

 Tests for state as a container of data
*/

// TPLs
#include "UnitTest++.h"

#include "IO.hh"
#include "MeshFactory.hh"
#include "TreeVector.hh"
#include "State.hh"
#include "errors.hh"

#include "Data_Helpers.hh"
#include "Vec.hh"

TEST(STATE_CREATION) {
  using namespace Amanzi;

  State s;
  s.Require<double>("my_double", Tags::DEFAULT, "my_double");
  s.Setup();

  CHECK(s.GetRecord("my_double", Tags::DEFAULT).ValidType<double>());
  CHECK(!s.GetRecord("my_double", Tags::DEFAULT).ValidType<int>());

  CHECK(s.GetRecordSet("my_double").ValidType<double>());
  CHECK(!s.GetRecordSet("my_double").ValidType<int>());
}


TEST(STATE_ASSIGNMENT) {
  using namespace Amanzi;

  State s;
  s.Require<double>("my_double", Tags::DEFAULT, "my_double");
  s.Setup();
  s.GetW<double>("my_double", "my_double") = 1.1;
  CHECK_EQUAL(1.1, s.Get<double>("my_double"));
}


TEST(STATE_FACTORIES_WITH_CREATE) {
  using namespace Amanzi;

  // create a mesh
  Comm_ptr_type comm = Amanzi::getDefaultComm();
  AmanziMesh::MeshFactory fac(comm);
  auto mesh = fac.create(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);

  // create a state
  State s;
  s.RegisterDomainMesh(mesh);

  // require data with factory
  s.Require<CompositeVector, CompositeVectorSpace>("my_vec", Tags::DEFAULT, "my_vec_owner")
      .SetMesh(s.GetMesh())
      ->SetGhosted();

  s.Require<CompositeVector, CompositeVectorSpace>("my_vec", Tags::DEFAULT).SetComponent(
      "cell", AmanziMesh::CELL, 1);

  s.Setup();
  Teuchos::RCP<CompositeVector> cv = s.GetPtrW<CompositeVector>("my_vec", Tags::DEFAULT, "my_vec_owner");
  // putting TreeVectors into state is a little tricky because they are
  // typically created from existing CompositeVectors.  We setup first, then
  // Require and stuff the pointer in.
  s.Require<TreeVector, TreeVectorSpace>("my_tree_vec", Tags::DEFAULT, "my_tree_vec");
  Teuchos::RCP<TreeVector> tv = Teuchos::rcp(new TreeVector());
  tv->SetData(cv);
  s.SetPtr<TreeVector>("my_tree_vec", Tags::DEFAULT, "my_tree_vec", tv);
}


TEST(STATE_FACTORIES_WITH_CONSTRUCTOR) {
  using namespace Amanzi;
  auto comm = Amanzi::getDefaultComm();
  State s;
  Epetra_Map my_map(-1, 3, 0, *comm);
  s.Require<Epetra_Vector, Epetra_Map>(my_map, "e_vec", Tags::DEFAULT, "e_vec");
  s.Setup();

  CHECK_EQUAL(comm->NumProc() * 3, s.Get<Epetra_Vector>("e_vec", Tags::DEFAULT).MyLength());
}


//
// NOTE: This cannot compile because there is no default constructor for Epetra_Map.  Yay!
//
// TEST(STATE_FACTORIES_WITH_CONSTRUCTOR2) {
//   using namespace Amanzi;
//   auto comm = Amanzi::getDefaultComm();
//   State s;
//   Epetra_Map my_map(-1, 3, 0, *comm);
//   s.Require<Epetra_Vector, Epetra_Map>("e_vec", Tags::DEFAULT, "e_vec");
//   s.Setup();

//   CHECK_EQUAL(comm->NumProc() * 3, s.Get<Epetra_Vector>("e_vec", Tags::DEFAULT).MyLength());
// }


TEST(STATE_HETEROGENEOUS_DATA) {
  using namespace Amanzi;

  // create a mesh
  Comm_ptr_type comm = Amanzi::getDefaultComm();
  AmanziMesh::MeshFactory fac(comm);
  auto mesh = fac.create(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);

  // create a state
  State s;
  s.RegisterDomainMesh(mesh);

  // require some data
  s.Require<double>("my_double", Tags::DEFAULT, "my_double_owner");

  // require a copy
  Tag tag_prev = make_tag("prev");
  s.Require<double>("my_double", tag_prev, "my_double_prev_owner");

  // require data with factory
  s.Require<CompositeVector, CompositeVectorSpace>("my_vec", Tags::DEFAULT, "my_vec_owner")
      .SetMesh(s.GetMesh())->SetComponent("cell", AmanziMesh::CELL, 1)->SetGhosted();

  s.Setup();

  // existence
  CHECK(s.HasRecord("my_double"));
  CHECK(s.HasRecord("my_vec"));
  CHECK(!s.HasRecord("my_nonexistent_data"));

  // defaults
  CHECK(!s.GetRecord("my_double").initialized());
  CHECK(s.GetRecord("my_double").io_checkpoint());
  CHECK(s.GetRecord("my_double").io_vis());

  // data access, construction
  CHECK(s.Get<CompositeVector>("my_vec").HasComponent("cell"));

  // incorrect type in Get
  CHECK_THROW(s.Get<double>("my_vec"), Errors::Message);

  // nonexistent data -- two checks ensure that previous HasRecord() call didn't
  // create the data!
  CHECK_THROW(s.Get<double>("my_nonexistent_data"), Errors::Message);
  CHECK_THROW(s.Get<double>("my_other_nonexistent_data"), Errors::Message);

  // setting data
  s.Assign("my_double", Tags::DEFAULT, "my_double_owner", 1.1);
  CHECK_EQUAL(1.1, s.Get<double>("my_double"));

  // copies
  CHECK(s.HasRecord("my_double", tag_prev));
  CHECK(!s.HasRecord("my_vec", tag_prev));

  s.Assign("my_double", tag_prev, "my_double_prev_owner", 2.2);
  CHECK_EQUAL(1.1, s.Get<double>("my_double"));
  CHECK_EQUAL(2.2, s.Get<double>("my_double", tag_prev));

  // set by reference
  s.GetW<double>("my_double", tag_prev, "my_double_prev_owner") = 3.3;
  CHECK_EQUAL(1.1, s.Get<double>("my_double"));
  CHECK_EQUAL(3.3, s.Get<double>("my_double", tag_prev));
}
