/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella, others
*/

#include <iostream>
#include <UnitTest++.h>

#include <AmanziComm.hh>

#include "dbc.hh"
#include "GeometricModel.hh"
#include "../MeshFrameworkFactory.hh"
#include "../MeshException.hh"

// Check to see if we have some files to read
#define BOGUS_TEST_FILE "test/bogus.exo"
#define EXODUS_TEST_FILE "test/hex_3x3x3_ss.exo"
#define NEMESIS_TEST_FILE "test/hex_10x10x10_ss.par"
#define MOAB_TEST_FILE "test/hex_3x3x3_ss_4P.h5m"

// -------------------------------------------------------------
// check_preference
// -------------------------------------------------------------

SUITE (MeshFramework)
{


  TEST (Generate3D)
  {
    auto comm = Amanzi::getDefaultComm();
    bool parallel(comm->NumProc() > 1);

    Amanzi::AmanziMesh::Preference pref;
    Amanzi::AmanziMesh::MeshFrameworkFactory meshfactory(comm);

    double x0( 0.0), y0( 0.0), z0( 0.0);
    double x1(10.0), y1(10.0), z1(10.0);
    int nx(10), ny(10), nz(10);

    // The Simple framework is always available, but will only
    // generate in serial
    pref.clear(); pref.push_back(Amanzi::AmanziMesh::Framework::SIMPLE);
    meshfactory.set_preference(pref);

    if (parallel) {
      Teuchos::RCP<Amanzi::AmanziMesh::MeshFramework> mesh;
      CHECK_THROW(mesh = meshfactory.create(x0, y0, z0,
                                      x1, y1, z1,
                                      nx, ny, nz),
                  Amanzi::AmanziMesh::Message);
    } else {
      Teuchos::RCP<Amanzi::AmanziMesh::MeshFramework> mesh;
      mesh = meshfactory.create(x0, y0, z0,
                          x1, y1, z1,
                          nx, ny, nz);
    }


    // The MSTK framework, if available, will always generate
    if (Amanzi::AmanziMesh::framework_enabled(
            Amanzi::AmanziMesh::Framework::MSTK)) {
      pref.clear(); pref.push_back(Amanzi::AmanziMesh::Framework::MSTK);
      meshfactory.set_preference(pref);
      Teuchos::RCP<Amanzi::AmanziMesh::MeshFramework> mesh =
          meshfactory.create(x0, y0, z0,
                             x1, y1, z1,
                             nx, ny, nz);
      CHECK(!mesh.is_null());
    }

    // The MOAB framework cannot generate
    if (Amanzi::AmanziMesh::framework_enabled(
            Amanzi::AmanziMesh::Framework::MOAB)) {
      pref.clear(); pref.push_back(Amanzi::AmanziMesh::Framework::MOAB);
      meshfactory.set_preference(pref);
      CHECK_THROW(auto mesh = meshfactory.create(x0, y0, z0,
              x1, y1, z1,
              nx, ny, nz),
                  Amanzi::AmanziMesh::Message);
    }
  }

  TEST (Generate2D)
  {
    auto comm = Amanzi::getDefaultComm();

    Amanzi::AmanziMesh::Preference pref;
    Teuchos::RCP<Amanzi::AmanziMesh::MeshFramework> mesh;
    Amanzi::AmanziMesh::MeshFrameworkFactory meshfactory(comm);

    double x0( 0.0), y0( 0.0);
    double x1(10.0), y1(10.0);
    int nx(10), ny(10);

    // The MSTK framework, if available, will always generate
    if (Amanzi::AmanziMesh::framework_enabled(
            Amanzi::AmanziMesh::Framework::MSTK)) {
      pref.clear(); pref.push_back(Amanzi::AmanziMesh::Framework::MSTK);
      meshfactory.set_preference(pref);
      mesh = meshfactory.create(x0, y0,
              x1, y1,
              nx, ny);
      CHECK(!mesh.is_null());
      mesh.reset();
    }

    // The Simple framework is always available, but
    // cannot generate 2D meshes
    pref.clear(); pref.push_back(Amanzi::AmanziMesh::Framework::SIMPLE);
    meshfactory.set_preference(pref);

    CHECK_THROW(mesh = meshfactory.create(x0, y0,
                                    x1, y1,
                                    nx, ny),
                Amanzi::AmanziMesh::Message);
    mesh.reset();

    // The MOAB framework cannot generate
    if (Amanzi::AmanziMesh::framework_enabled(
            Amanzi::AmanziMesh::Framework::MOAB)) {
      pref.clear(); pref.push_back(Amanzi::AmanziMesh::Framework::MOAB);
      meshfactory.set_preference(pref);
      CHECK_THROW(mesh = meshfactory.create(x0, y0,
              x1, y1,
              nx, ny),
                  Amanzi::AmanziMesh::Message);
      mesh.reset();
    }
  }


  TEST (ParameterGenerate3)
  {
    auto comm = Amanzi::getDefaultComm();
    bool parallel(comm->NumProc() > 1);

    Amanzi::AmanziMesh::Preference pref;
    Teuchos::RCP<Amanzi::AmanziMesh::MeshFramework> mesh;
    Amanzi::AmanziMesh::MeshFrameworkFactory meshfactory(comm);

    // make a parameter list to try out

    Teuchos::ParameterList parameter_list;
    Teuchos::Array<int> ncells = Teuchos::tuple<int>(10,10,10);
    Teuchos::Array<double> low = Teuchos::tuple<double>(0.0, 0.0, 0.0);
    Teuchos::Array<double> high = Teuchos::tuple<double>(1.0, 1.0, 1.0);
    parameter_list.set< Teuchos::Array<int> >("number of cells", ncells);
    parameter_list.set< Teuchos::Array<double> >("domain low coordinate", low);
    parameter_list.set< Teuchos::Array<double> >("domain high coordinate", high);

    pref.clear(); pref.push_back(Amanzi::AmanziMesh::Framework::SIMPLE);
    meshfactory.set_preference(pref);

    if (parallel) {
      CHECK_THROW(mesh = meshfactory.create(parameter_list),
                  Amanzi::AmanziMesh::Message);
      mesh.reset();
    } else {
      mesh = meshfactory.create(parameter_list);
      CHECK(!mesh.is_null());
      mesh.reset();
    }

    if (Amanzi::AmanziMesh::framework_enabled(
            Amanzi::AmanziMesh::Framework::MSTK)) {
      pref.clear(); pref.push_back(Amanzi::AmanziMesh::Framework::MSTK);
      meshfactory.set_preference(pref);
      mesh = meshfactory.create(parameter_list);
      CHECK(!mesh.is_null());
      mesh.reset();
    }

    // Other frameworks can't generate, so they should throw
    if (Amanzi::AmanziMesh::framework_enabled(
            Amanzi::AmanziMesh::Framework::MOAB)) {
      pref.clear();
      pref.push_back(Amanzi::AmanziMesh::Framework::MOAB);
      meshfactory.set_preference(pref);
      CHECK_THROW(mesh = meshfactory.create(parameter_list),
                  Amanzi::AmanziMesh::Message);
    }
  }


  TEST (ParameterGenerate2)
  {
    auto comm = Amanzi::getDefaultComm();

    Amanzi::AmanziMesh::Preference pref;
    Teuchos::RCP<Amanzi::AmanziMesh::MeshFramework> mesh;
    Amanzi::AmanziMesh::MeshFrameworkFactory meshfactory(comm);

    // make a parameter list to try out

    Teuchos::ParameterList parameter_list;
    Teuchos::Array<int> ncells = Teuchos::tuple<int>(10,10);
    Teuchos::Array<double> low = Teuchos::tuple<double>(0.0, 0.0);
    Teuchos::Array<double> high = Teuchos::tuple<double>(1.0, 1.0);
    parameter_list.set< Teuchos::Array<int> >("number of cells", ncells);
    parameter_list.set< Teuchos::Array<double> >("domain low coordinate", low);
    parameter_list.set< Teuchos::Array<double> >("domain high coordinate", high);


    // MSTK, if available, can generate 2D meshes
    if (Amanzi::AmanziMesh::framework_enabled(
            Amanzi::AmanziMesh::Framework::MSTK)) {
      pref.clear(); pref.push_back(Amanzi::AmanziMesh::Framework::MSTK);
      meshfactory.set_preference(pref);
      mesh = meshfactory.create(parameter_list);
      CHECK(!mesh.is_null());
      mesh.reset();
    }

    // Simple mesh is always available but cannot generate 2D meshes
    pref.clear(); pref.push_back(Amanzi::AmanziMesh::Framework::SIMPLE);
    meshfactory.set_preference(pref);

    CHECK_THROW(mesh = meshfactory.create(parameter_list),
                Amanzi::AmanziMesh::Message);
    mesh.reset();

    // Other frameworks can't generate, so they should throw
    if (Amanzi::AmanziMesh::framework_enabled(
            Amanzi::AmanziMesh::Framework::MOAB)) {
      pref.clear();
      pref.push_back(Amanzi::AmanziMesh::Framework::MOAB);
      meshfactory.set_preference(pref);
      CHECK_THROW(mesh = meshfactory.create(parameter_list),
                  Amanzi::AmanziMesh::Message);
    }
  }


  // The Simple framework cannot read anything, even if it exists
  TEST (ReadSimple) {
    auto comm = Amanzi::getDefaultComm();

    Amanzi::AmanziMesh::Preference pref;
    Teuchos::RCP<Amanzi::AmanziMesh::MeshFramework> mesh;
    Amanzi::AmanziMesh::MeshFrameworkFactory meshfactory(comm);

    pref.clear(); pref.push_back(Amanzi::AmanziMesh::Framework::SIMPLE);
    meshfactory.set_preference(pref);
    CHECK_THROW(mesh = meshfactory.create(BOGUS_TEST_FILE),
                Amanzi::AmanziMesh::Message);
    CHECK_THROW(mesh = meshfactory.create(MOAB_TEST_FILE),
                Amanzi::AmanziMesh::Message);
    CHECK_THROW(mesh = meshfactory.create(EXODUS_TEST_FILE),
                Amanzi::AmanziMesh::Message);
    CHECK_THROW(mesh = meshfactory.create(NEMESIS_TEST_FILE),
                Amanzi::AmanziMesh::Message);
  }

  // Try to read a MOAB HDF5 file, which can only be read by the MOAB
  // framework
  TEST (ReadMOABHDF5)
  {
    auto comm = Amanzi::getDefaultComm();

    Teuchos::RCP<Amanzi::AmanziMesh::MeshFramework> mesh;
    Amanzi::AmanziMesh::MeshFrameworkFactory meshfactory(comm);

    if (Amanzi::AmanziMesh::framework_enabled(
            Amanzi::AmanziMesh::Framework::MOAB)) {
      mesh = meshfactory.create(MOAB_TEST_FILE);
      CHECK(!mesh.is_null());
    } else {
      CHECK_THROW(mesh = meshfactory.create(MOAB_TEST_FILE),
                  Amanzi::AmanziMesh::Message);
    }
  }

  TEST (ReadExodus)
  {
    auto comm = Amanzi::getDefaultComm();
    bool parallel(comm->NumProc() > 1);

    Teuchos::RCP<Amanzi::AmanziMesh::MeshFramework> mesh;
    Amanzi::AmanziMesh::MeshFrameworkFactory meshfactory(comm);

    bool available = false;
    if (Amanzi::AmanziMesh::framework_enabled(
            Amanzi::AmanziMesh::Framework::MSTK)) {
      available = true;
    }

    if (Amanzi::AmanziMesh::framework_enabled(
            Amanzi::AmanziMesh::Framework::MOAB)) {
      if (!parallel) available = true;
    }

    if (available) {
      mesh = meshfactory.create(EXODUS_TEST_FILE);
      CHECK(!mesh.is_null());
    } else {
      CHECK_THROW(mesh = meshfactory.create(EXODUS_TEST_FILE),
                  Amanzi::AmanziMesh::Message);
    }
  }

  TEST (ReadNemesis)
  {
    auto comm = Amanzi::getDefaultComm();
    bool parallel(comm->NumProc() > 1);

    Teuchos::RCP<Amanzi::AmanziMesh::MeshFramework> mesh;
    Amanzi::AmanziMesh::MeshFrameworkFactory meshfactory(comm);

    bool available = false;
    if (Amanzi::AmanziMesh::framework_enabled(
            Amanzi::AmanziMesh::Framework::MSTK)) {
      if (parallel) available = true;
    }
    if (available) {
      mesh = meshfactory.create(NEMESIS_TEST_FILE);
      CHECK(!mesh.is_null());
    } else {
      CHECK_THROW(mesh = meshfactory.create(NEMESIS_TEST_FILE),
                  Amanzi::AmanziMesh::Message);
    }
  }


  TEST (Extract3)
  {
    auto comm = Amanzi::getDefaultComm();
    bool parallel(comm->NumProc() > 1);

    double x0( 0.0), y0( 0.0), z0( 0.0);
    double x1(10.0), y1(10.0), z1(10.0);
    int nx(10), ny(10), nz(10);

  // create a sublist name Regions and put a reference to it in
  // reg_spec and other sublists as references. Turns out it is
  // important to define reg_spec and other lists below as references
  // - otherwise, a new copy is made of the sublist that is retrieved

    Teuchos::ParameterList parameterlist;

    Teuchos::ParameterList& reg_spec = parameterlist.sublist("regions");

    Teuchos::ParameterList& top_surface = reg_spec.sublist("Top Surface");
    Teuchos::ParameterList& top_surface_def = top_surface.sublist("region: plane");
    Teuchos::Array<double> loc1 = Teuchos::tuple(0.0,0.0,10.0);
    Teuchos::Array<double> dir1 = Teuchos::tuple(0.0,0.0,-1.0);
    top_surface_def.set< Teuchos::Array<double> >("point",loc1);
    top_surface_def.set< Teuchos::Array<double> >("normal",dir1);

    Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
        Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, reg_spec, *comm));

    std::string topsurfname("Top Surface");
    std::vector<std::string> setnames;
    setnames.push_back(topsurfname);

    Amanzi::AmanziMesh::Preference pref;
    Amanzi::AmanziMesh::MeshFrameworkFactory meshfactory(comm,gm);

    bool flatten = true;
    bool extrude = false;

    // Simple mesh CANNOT extract a mesh from another mesh

    pref.clear(); pref.push_back(Amanzi::AmanziMesh::Framework::SIMPLE);
    meshfactory.set_preference(pref);

    if (parallel) {
      Teuchos::RCP<Amanzi::AmanziMesh::MeshFramework> mesh, newmesh;
      CHECK_THROW(mesh = meshfactory.create(x0, y0, z0,
                                      x1, y1, z1,
                                      nx, ny, nz),
                  Amanzi::AmanziMesh::Message);
    } else {
      Teuchos::RCP<Amanzi::AmanziMesh::MeshFramework> mesh, newmesh;
      mesh = meshfactory.create(x0, y0, z0,
              x1, y1, z1,
              nx, ny, nz);

      CHECK(!mesh.is_null());
      Amanzi::AmanziMesh::Entity_ID_List ids("ids", 1);
      auto mesh_cache = Teuchos::rcp(new Amanzi::AmanziMesh::Mesh(mesh)); 
      CHECK_THROW(newmesh = meshfactory.create(mesh_cache,ids,
              Amanzi::AmanziMesh::Entity_kind::FACE,
              flatten),
                  Amanzi::AmanziMesh::Message);
    }


    // The MSTK framework, if available, will always generate
    if (Amanzi::AmanziMesh::framework_enabled(
            Amanzi::AmanziMesh::Framework::MSTK)) {
      Teuchos::RCP<Amanzi::AmanziMesh::MeshFramework> mesh, newmesh;
      pref.clear(); pref.push_back(Amanzi::AmanziMesh::Framework::MSTK);
      meshfactory.set_preference(pref);
      mesh = meshfactory.create(x0, y0, z0,
              x1, y1, z1,
              nx, ny, nz);
      CHECK(!mesh.is_null());

      Amanzi::AmanziMesh::Entity_ID_List ents("ents", 1);
      auto mesh_cache = Teuchos::rcp(new Amanzi::AmanziMesh::Mesh(mesh)); 
      newmesh = meshfactory.create(mesh_cache,ents,Amanzi::AmanziMesh::Entity_kind::FACE,flatten);
      CHECK(!newmesh.is_null());
    }

    // The MOAB framework cannot generate
    if (Amanzi::AmanziMesh::framework_enabled(
            Amanzi::AmanziMesh::Framework::MOAB)) {
      Teuchos::RCP<Amanzi::AmanziMesh::MeshFramework> mesh;
      pref.clear(); pref.push_back(Amanzi::AmanziMesh::Framework::MOAB);
      meshfactory.set_preference(pref);
      CHECK_THROW(mesh = meshfactory.create(x0, y0, z0,
                                      x1, y1, z1,
                                      nx, ny, nz),
                  Amanzi::AmanziMesh::Message);
      mesh.reset();
    }
  }


} // SUITE
