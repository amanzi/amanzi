#include <iostream>
#include <UnitTest++.h>

#include "AmanziTypes.hh"
#include "Teuchos_DefaultMpiComm.hpp"
#include "dbc.hh"
#include "../MeshFactory.hh"
#include "../FrameworkTraits.hh"

// Check to see if we have some files to read

#ifndef BOGUS_TEST_FILE
#define BOGUS_TEST_FILE "bogus.file"
#endif

#ifndef EXODUS_TEST_FILE
#error EXODUS_TEST_FILE must be defined
#endif

#ifndef NEMESIS_TEST_FILE
#error EXODUS_TEST_FILE must be defined
#endif

#ifndef MOAB_TEST_FILE
#error EXODUS_TEST_FILE must be defined
#endif

using namespace Amanzi;

// -------------------------------------------------------------
// check_preference
// -------------------------------------------------------------

static void
check_preference(AmanziMesh::MeshFactory& mesh_factory, 
                 const AmanziMesh::Framework& f)
{
  AmanziMesh::FrameworkPreference pref;
  pref.push_back(f);
  if (AmanziMesh::framework_available(f)) {
    mesh_factory.preference(pref);
    CHECK(mesh_factory.preference().front() == f);
  } else {
    CHECK_THROW(mesh_factory.preference(pref), AmanziMesh::Message);
  }
}


SUITE (MeshFramework)
{

  // This tests setting the Mesh Factory framework preferences.  If
  // only one framework is preferred, and it is not available, an
  // exception should be thrown, while setting preferences
  TEST (PreferenceThrow)
  {
    auto comm = Comm_ptr_type(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
    bool parallel(comm->getSize() > 1);
    
    AmanziMesh::MeshFactory mesh_factory(comm);
    AmanziMesh::FrameworkPreference pref(mesh_factory.preference());

    // The Simple framework should always be there
    check_preference(mesh_factory, AmanziMesh::Simple);
    check_preference(mesh_factory, AmanziMesh::MOAB);
    check_preference(mesh_factory, AmanziMesh::STKMESH);
    check_preference(mesh_factory, AmanziMesh::MSTK);
  }
    
  TEST (Generate)
  {
    auto comm = Comm_ptr_type(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
    bool parallel(comm->getSize() > 1);
    
    AmanziMesh::FrameworkPreference pref;
    Teuchos::RCP<AmanziMesh::Mesh> mesh;
    AmanziMesh::MeshFactory mesh_factory(comm);

    double x0( 0.0), y0( 0.0), z0( 0.0);
    double x1(10.0), y1(10.0), z1(10.0);
    int nx(10), ny(10), nz(10);

    // The Simple framework is always available, but will only
    // generate in serial

    pref.clear(); pref.push_back(AmanziMesh::Simple);
    mesh_factory.preference(pref);

    if (parallel) {
      CHECK_THROW(mesh = mesh_factory(x0, y0, z0,
                                      x1, y1, z1,
                                      nx, ny, nz),
                  AmanziMesh::Message);
      mesh.reset();
    } else {
      mesh = mesh_factory(x0, y0, z0,
                          x1, y1, z1,
                          nx, ny, nz);
      CHECK(!mesh.is_null());
      mesh.reset();
    }

    // The STK, if available, framework will always generate

    if (framework_available(AmanziMesh::STKMESH)) {
      pref.clear(); pref.push_back(AmanziMesh::STKMESH);
      mesh_factory.preference(pref);
      mesh = mesh_factory(x0, y0, z0,
                          x1, y1, z1,
                          nx, ny, nz);
      CHECK(!mesh.is_null());
      mesh.reset();
    }

    // The MSTK framework, if available, will always generate

    if (framework_available(AmanziMesh::MSTK)) {
      pref.clear(); pref.push_back(AmanziMesh::MSTK);
      mesh_factory.preference(pref);
      mesh = mesh_factory(x0, y0, z0,
                          x1, y1, z1,
                          nx, ny, nz);
      CHECK(!mesh.is_null());
      mesh.reset();
    }

    // The MOAB framework cannot generate


    if (framework_available(AmanziMesh::MOAB)) {
      pref.clear(); pref.push_back(AmanziMesh::MOAB);
      mesh_factory.preference(pref);
      CHECK_THROW(mesh = mesh_factory(x0, y0, z0,
                                      x1, y1, z1,
                                      nx, ny, nz),
                  AmanziMesh::Message);
      mesh.reset();
    }
  }

  TEST (Generate2D)
  {
    auto comm = Comm_ptr_type(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
    bool parallel(comm->getSize() > 1);
    
    AmanziMesh::FrameworkPreference pref;
    Teuchos::RCP<AmanziMesh::Mesh> mesh;
    AmanziMesh::MeshFactory mesh_factory(comm);

    double x0( 0.0), y0( 0.0);
    double x1(10.0), y1(10.0);
    int nx(10), ny(10);

    // The MSTK framework, if available, will always generate

    if (framework_available(AmanziMesh::MSTK)) {
      pref.clear(); pref.push_back(AmanziMesh::MSTK);
      mesh_factory.preference(pref);
      mesh = mesh_factory(x0, y0,
                          x1, y1,
                          nx, ny);
      CHECK(!mesh.is_null());
      mesh.reset();
    }

    // The Simple framework is always available, but 
    // cannot generate 2D meshes

    pref.clear(); pref.push_back(AmanziMesh::Simple);
    mesh_factory.preference(pref);

    CHECK_THROW(mesh = mesh_factory(x0, y0,
                                    x1, y1,
                                    nx, ny),
                AmanziMesh::Message);
    mesh.reset();

    // The STK, even if available cannot generate 2D meshes

    if (framework_available(AmanziMesh::STKMESH)) {
      pref.clear(); pref.push_back(AmanziMesh::STKMESH);
      mesh_factory.preference(pref);
      CHECK_THROW(mesh = mesh_factory(x0, y0,
                                      x1, y1,
                                      nx, ny),
                  AmanziMesh::Message);
      mesh.reset();
    }

    // The MOAB framework cannot generate


    if (framework_available(AmanziMesh::MOAB)) {
      pref.clear(); pref.push_back(AmanziMesh::MOAB);
      mesh_factory.preference(pref);
      CHECK_THROW(mesh = mesh_factory(x0, y0,
                                      x1, y1,
                                      nx, ny),
                  AmanziMesh::Message);
      mesh.reset();
    }
  }

  TEST (ParameterGenerate3)
  {
    auto comm = Comm_ptr_type(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
    bool parallel(comm->getSize() > 1);
    
    AmanziMesh::FrameworkPreference pref;
    Teuchos::RCP<AmanziMesh::Mesh> mesh;
    AmanziMesh::MeshFactory mesh_factory(comm);

    // make a parameter list to try out
    
    Teuchos::ParameterList parameter_list;
    Teuchos::Array<int> ncells = Teuchos::tuple<int>(10,10,10);
    Teuchos::Array<double> low = Teuchos::tuple<double>(0.0, 0.0, 0.0);
    Teuchos::Array<double> high = Teuchos::tuple<double>(1.0, 1.0, 1.0);
    parameter_list.set< Teuchos::Array<int> >("number of cells", ncells);
    parameter_list.set< Teuchos::Array<double> >("domain low coordinate", low);
    parameter_list.set< Teuchos::Array<double> >("domain high coordinate", high);
    
    pref.clear(); pref.push_back(AmanziMesh::Simple);
    mesh_factory.preference(pref);

    if (parallel) {
      CHECK_THROW(mesh = mesh_factory(parameter_list),
                  AmanziMesh::Message);
      mesh.reset();
    } else {
      mesh = mesh_factory(parameter_list);
      CHECK(!mesh.is_null());
      mesh.reset();
    }

    // The STK, if available, framework will always generate

    if (framework_available(AmanziMesh::STKMESH)) {
      pref.clear(); pref.push_back(AmanziMesh::STKMESH);
      mesh_factory.preference(pref);
      mesh = mesh_factory(parameter_list);
      CHECK(!mesh.is_null());
      mesh.reset();
    }

    if (framework_available(AmanziMesh::MSTK)) {
      pref.clear(); pref.push_back(AmanziMesh::MSTK);
      mesh_factory.preference(pref);
      mesh = mesh_factory(parameter_list);
      CHECK(!mesh.is_null());
      mesh.reset();
    }

    // Other frameworks can't generate, so they should throw
    pref.clear(); 
    if (framework_available(AmanziMesh::MOAB)) {
      pref.push_back(AmanziMesh::MOAB);
    }
    mesh_factory.preference(pref);
    CHECK_THROW(mesh = mesh_factory(parameter_list),
                AmanziMesh::Message);
  }


  TEST (ParameterGenerate2)
  {
    auto comm = Comm_ptr_type(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
    bool parallel(comm->getSize() > 1);
    
    AmanziMesh::FrameworkPreference pref;
    Teuchos::RCP<AmanziMesh::Mesh> mesh;
    AmanziMesh::MeshFactory mesh_factory(comm);

    // make a parameter list to try out
    
    Teuchos::ParameterList parameter_list;
    Teuchos::Array<int> ncells = Teuchos::tuple<int>(10,10);
    Teuchos::Array<double> low = Teuchos::tuple<double>(0.0, 0.0);
    Teuchos::Array<double> high = Teuchos::tuple<double>(1.0, 1.0);
    parameter_list.set< Teuchos::Array<int> >("number of cells", ncells);
    parameter_list.set< Teuchos::Array<double> >("domain low coordinate", low);
    parameter_list.set< Teuchos::Array<double> >("domain high coordinate", high);
    

    // MSTK, if available, can generate 2D meshes

    if (framework_available(AmanziMesh::MSTK)) {
      pref.clear(); pref.push_back(AmanziMesh::MSTK);
      mesh_factory.preference(pref);
      mesh = mesh_factory(parameter_list);
      CHECK(!mesh.is_null());
      mesh.reset();
    }

    // Simple mesh is always available but cannot generate 2D meshes

    pref.clear(); pref.push_back(AmanziMesh::Simple);
    mesh_factory.preference(pref);

    CHECK_THROW(mesh = mesh_factory(parameter_list),
                AmanziMesh::Message);
    mesh.reset();

    // The STK, even if available, cannot generate 2D meshes

    if (framework_available(AmanziMesh::STKMESH)) {
      pref.clear(); pref.push_back(AmanziMesh::STKMESH);
      mesh_factory.preference(pref);
      CHECK_THROW(mesh = mesh_factory(parameter_list),
                  AmanziMesh::Message);
      mesh.reset();
    }

    // Other frameworks can't generate, so they should throw
    pref.clear(); 
    if (framework_available(AmanziMesh::MOAB)) {
      pref.push_back(AmanziMesh::MOAB);
    }
    mesh_factory.preference(pref);
    CHECK_THROW(mesh = mesh_factory(parameter_list),
                AmanziMesh::Message);
  }


  // The Simple framework cannot read anything, even if it exists
  TEST (ReadSimple) {
    auto comm = Comm_ptr_type(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
    bool parallel(comm->getSize() > 1);
    
    AmanziMesh::FrameworkPreference pref;
    Teuchos::RCP<AmanziMesh::Mesh> mesh;
    AmanziMesh::MeshFactory mesh_factory(comm);

    pref.clear(); pref.push_back(AmanziMesh::Simple);
    mesh_factory.preference(pref);
    CHECK_THROW(mesh = mesh_factory(BOGUS_TEST_FILE),
                AmanziMesh::Message);
    CHECK_THROW(mesh = mesh_factory(MOAB_TEST_FILE),
                AmanziMesh::Message);
    CHECK_THROW(mesh = mesh_factory(EXODUS_TEST_FILE),
                AmanziMesh::Message);
    CHECK_THROW(mesh = mesh_factory(NEMESIS_TEST_FILE),
                AmanziMesh::Message);
  }

  // Try to read a MOAB HDF5 file, which can only be read by the MOAB
  // framework
  TEST (ReadMOABHDF5)
  {
    auto comm = Comm_ptr_type(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
    bool parallel(comm->getSize() > 1);

    Teuchos::RCP<AmanziMesh::Mesh> mesh;
    AmanziMesh::MeshFactory mesh_factory(comm);

    if (AmanziMesh::framework_available(AmanziMesh::MOAB)) {
      mesh = mesh_factory(MOAB_TEST_FILE);
      CHECK(!mesh.is_null());
    } else {
      CHECK_THROW(mesh = mesh_factory(MOAB_TEST_FILE),
                  AmanziMesh::Message);
    }

    // Try it with another framework just for grins

    if (AmanziMesh::framework_available(AmanziMesh::STKMESH)) {
      AmanziMesh::FrameworkPreference pref;
      pref.push_back(AmanziMesh::STKMESH);
      mesh_factory.preference(pref);
      CHECK_THROW(mesh = mesh_factory(MOAB_TEST_FILE),
                  AmanziMesh::Message);
    }
  }

  TEST (ReadExodus) 
  {
    auto comm = Comm_ptr_type(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
    bool parallel(comm->getSize() > 1);

    Teuchos::RCP<AmanziMesh::Mesh> mesh;
    AmanziMesh::MeshFactory mesh_factory(comm);

    bool available =
        (AmanziMesh::framework_available(AmanziMesh::STKMESH) && !parallel) ||
        AmanziMesh::framework_available(AmanziMesh::MSTK);

    if (available) {
      mesh = mesh_factory(EXODUS_TEST_FILE);
      CHECK(!mesh.is_null());
    } else {
      CHECK_THROW(mesh = mesh_factory(EXODUS_TEST_FILE),
                  AmanziMesh::Message);
    }
  }
  TEST (ReadNemesis) 
  {
    auto comm = Comm_ptr_type(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
    bool parallel(comm->getSize() > 1);
    
    Teuchos::RCP<AmanziMesh::Mesh> mesh;
    AmanziMesh::MeshFactory mesh_factory(comm);
    if ((AmanziMesh::framework_available(AmanziMesh::STKMESH) ||
         AmanziMesh::framework_available(AmanziMesh::MSTK)) && 
        parallel) {
      mesh = mesh_factory(NEMESIS_TEST_FILE);
      CHECK(!mesh.is_null());
    } else {
      CHECK_THROW(mesh = mesh_factory(NEMESIS_TEST_FILE),
                  AmanziMesh::Message);
    }
  }      


  TEST (Extract3)
  {
    auto comm = Comm_ptr_type(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
    bool parallel(comm->getSize() > 1);
    
    AmanziMesh::FrameworkPreference pref;
    Teuchos::RCP<AmanziMesh::Mesh> mesh, newmesh;
    AmanziMesh::MeshFactory mesh_factory(comm);

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

    Teuchos::RCP<AmanziGeometry::GeometricModel> gm =
        Teuchos::rcp(new AmanziGeometry::GeometricModel(3, reg_spec, comm));

    std::string topsurfname("Top Surface");
    std::vector<std::string> setnames;
    setnames.push_back(topsurfname);

    bool flatten = true;
    bool extrude = false;

    // Simple mesh CANNOT extract a mesh from another mesh

    pref.clear(); pref.push_back(AmanziMesh::Simple);
    mesh_factory.preference(pref);

    if (parallel) {
      CHECK_THROW(mesh = mesh_factory(x0, y0, z0,
                                      x1, y1, z1,
                                      nx, ny, nz,gm),
                  AmanziMesh::Message);
      mesh.reset();
    } else {
      mesh = mesh_factory(x0, y0, z0,
                          x1, y1, z1,
                          nx, ny, nz);

      CHECK_THROW(newmesh = mesh_factory(mesh.get(),setnames,
                                         AmanziMesh::FACE,
                                         flatten,extrude),
                  AmanziMesh::Message);
      mesh.reset();
      newmesh.reset();
    }

    // The STK, if available, framework will always generate

    if (framework_available(AmanziMesh::STKMESH)) {
      pref.clear(); pref.push_back(AmanziMesh::STKMESH);
      mesh_factory.preference(pref);
      mesh = mesh_factory(x0, y0, z0,
                          x1, y1, z1,
                          nx, ny, nz,gm);
      CHECK(!mesh.is_null());

      CHECK_THROW(newmesh = mesh_factory(mesh.get(),setnames,
                                         AmanziMesh::FACE,
                                         flatten,extrude),
                  AmanziMesh::Message);
      mesh.reset();
      newmesh.reset();
    }

    // The MSTK framework, if available, will always generate

    if (framework_available(AmanziMesh::MSTK)) {
      pref.clear(); pref.push_back(AmanziMesh::MSTK);
      mesh_factory.preference(pref);
      mesh = mesh_factory(x0, y0, z0,
                          x1, y1, z1,
                          nx, ny, nz,gm);
      CHECK(!mesh.is_null());

      newmesh = mesh_factory(mesh.get(),setnames,AmanziMesh::FACE,
                             flatten,extrude);
      CHECK(!newmesh.is_null());
      mesh.reset();
      newmesh.reset();
    }

    // The MOAB framework cannot generate


    if (framework_available(AmanziMesh::MOAB)) {
      pref.clear(); pref.push_back(AmanziMesh::MOAB);
      mesh_factory.preference(pref);
      CHECK_THROW(mesh = mesh_factory(x0, y0, z0,
                                      x1, y1, z1,
                                      nx, ny, nz,gm),
                  AmanziMesh::Message);
      mesh.reset();
    }
  }

    
}
