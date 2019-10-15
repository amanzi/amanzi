/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! <MISSING_ONELINE_DOCSTRING>

// Tests for state as a container of data

// TPLs
#include "UnitTest++.h"

#include "MeshFactory.hh"
#include "State.hh"
#include "errors.hh"

//#include "Op_Cell_Cell.hh"
//#include "Op_Factory.hh"
#include "Vec.hh"

TEST(STATE_CREATION)
{
  using namespace Amanzi;

  State s;
  s.Require<double>("my_double", "", "my_double");
  s.Setup();
}

TEST(STATE_ASSIGNMENT)
{
  using namespace Amanzi;

  State s;
  s.Require<double>("my_double", "", "my_double");
  s.Setup();
  s.GetW<double>("my_double", "my_double") = 1.1;
  CHECK_EQUAL(1.1, s.Get<double>("my_double"));
}

TEST(STATE_FACTORIES_PERSIST)
{
  using namespace Amanzi;

  // create a mesh
  auto comm = Comm_ptr_type(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
  AmanziMesh::MeshFactory fac(comm);
  auto mesh = fac.create(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);

  // create a state
  State s;
  s.RegisterDomainMesh(mesh);

  // require data with factory
  s.Require<CompositeVector, CompositeVectorSpace>("my_vec", "", "my_vec_owner")
    .SetMesh(s.GetMesh())
    ->SetGhosted();

  s.Require<CompositeVector, CompositeVectorSpace>("my_vec").SetComponent(
    "cell", AmanziMesh::CELL, 1);

  s.Setup();
}

TEST(STATE_HETEROGENEOUS_DATA)
{
  using namespace Amanzi;

  // create a mesh
  auto comm = Comm_ptr_type(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
  AmanziMesh::MeshFactory fac(comm);
  auto mesh = fac.create(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);

  // create a state
  State s;
  s.RegisterDomainMesh(mesh);

  // require some data
  s.Require<double>("my_double", "", "my_double_owner");

  // require a copy
  s.Require<double>("my_double", "prev", "my_double_prev_owner");

  // require data with factory
  s.Require<CompositeVector, CompositeVectorSpace>("my_vec", "", "my_vec_owner")
    .SetMesh(s.GetMesh())
    ->SetComponent("cell", AmanziMesh::CELL, 1)
    ->SetGhosted();

  s.Setup();

  // existence
  CHECK(s.HasData("my_double"));
  CHECK(s.HasData("my_vec"));
  CHECK(!s.HasData("my_nonexistent_data"));

  // defaults
  CHECK(!s.GetRecord("my_double").initialized());
  CHECK(s.GetRecord("my_double").io_checkpoint());
  CHECK(s.GetRecord("my_double").io_vis());

  // data access, construction
  CHECK(s.Get<CompositeVector>("my_vec").HasComponent("cell"));

  // incorrect type in Get
  CHECK_THROW(s.Get<double>("my_vec"), Errors::Message);

  // nonexistent data -- two checks ensure that previous HasData() call didn't
  // create the data!
  CHECK_THROW(s.Get<double>("my_nonexistent_data"), std::out_of_range);
  CHECK_THROW(s.Get<double>("my_other_nonexistent_data"), std::out_of_range);

  // setting data
  s.Set("my_double", "", "my_double_owner", 1.1);
  CHECK_EQUAL(1.1, s.Get<double>("my_double"));

  // copies
  CHECK(s.HasData("my_double", "prev"));
  CHECK(!s.HasData("my_vec", "prev"));

  s.Set("my_double", "prev", "my_double_prev_owner", 2.2);
  CHECK_EQUAL(1.1, s.Get<double>("my_double"));
  CHECK_EQUAL(2.2, s.Get<double>("my_double", "prev"));

  // set by reference
  s.GetW<double>("my_double", "prev", "my_double_prev_owner") = 3.3;
  CHECK_EQUAL(1.1, s.Get<double>("my_double"));
  CHECK_EQUAL(3.3, s.Get<double>("my_double", "prev"));
}
