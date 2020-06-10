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

#include "RegionBox.hh"
#include "Patch.hh"
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
  CHECK_THROW(s.Get<double>("my_nonexistent_data"), Errors::Message);
  CHECK_THROW(s.Get<double>("my_other_nonexistent_data"), Errors::Message);

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



TEST(STATE_PATCH_DATA)
{
  using namespace Amanzi;

  auto comm = Comm_ptr_type(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
  
  auto gm = Teuchos::rcp(new AmanziGeometry::GeometricModel(2));
  gm->AddRegion(Teuchos::rcp(new AmanziGeometry::RegionBox("box1", 0,
          AmanziGeometry::Point(0.,0.),
          AmanziGeometry::Point(2.0, 4.0))));
  gm->AddRegion(Teuchos::rcp(new AmanziGeometry::RegionBox("box2", 1,
          AmanziGeometry::Point(2.,0.),
          AmanziGeometry::Point(4.0, 4.0))));
  
  // create a mesh
  AmanziMesh::MeshFactory fac(comm, gm);
  auto mesh = fac.create(0.0, 0.0, 4.0, 4.0, 2, 2);

  // create a state
  State s;
  s.RegisterDomainMesh(mesh);

  auto& mps = s.Require<MultiPatch,MultiPatchSpace>("multipatch", "", "multipatch");
  mps.set_mesh(mesh);
  mps.ghosted = false;
  mps.AddPatch("box1", AmanziMesh::FACE, 1);
  mps.AddPatch("box2", AmanziMesh::FACE, 1);

  s.Setup();

  // existence
  CHECK(s.HasData("multipatch"));

  {
    auto& mp = s.GetW<MultiPatch>("multipatch","","multipatch");
    CHECK_EQUAL(7, mp[0].data.extent(0));
    CHECK_EQUAL(1, mp[0].data.extent(1));
    mp[0].data(3,0) = 1.0;
  }

  {
    const auto& mp = s.Get<MultiPatch>("multipatch", "");
    CHECK_EQUAL(1.0, mp[0].data(3,0));
  }  
}
                    
             
