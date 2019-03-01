// -------------------------------------------------------------
// file: test_meshfactory.cc
// -------------------------------------------------------------

#include <iostream>
#include <UnitTest++.h>

#include <AmanziComm.hh>

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

// -------------------------------------------------------------
// check_preference
// -------------------------------------------------------------

static void
check_preference(Amanzi::AmanziMesh::MeshFactory& meshfactory, 
                 const Amanzi::AmanziMesh::Framework& f)
{
  Amanzi::AmanziMesh::FrameworkPreference pref;
  pref.push_back(f);
  if (Amanzi::AmanziMesh::framework_available(f)) {
    meshfactory.set_preference(pref);
    CHECK(meshfactory.preference().front() == f);
  } else {
    CHECK_THROW(meshfactory.set_preference(pref), Amanzi::AmanziMesh::Message);
  }
}


SUITE (MeshFramework)
{

  // This tests setting the Mesh Factory framework preferences.  If
  // only one framework is preferred, and it is not available, an
  // exception should be thrown, while setting preferences
  TEST (PreferenceThrow)
  {
    auto comm = Amanzi::getDefaultComm();
    bool parallel(comm->NumProc() > 1);
    
    Amanzi::AmanziMesh::MeshFactory meshfactory(comm);
    Amanzi::AmanziMesh::FrameworkPreference pref(meshfactory.preference());

    // The Simple framework should always be there
    check_preference(meshfactory, Amanzi::AmanziMesh::Simple);
    check_preference(meshfactory, Amanzi::AmanziMesh::MOAB);
    check_preference(meshfactory, Amanzi::AmanziMesh::STKMESH);
    check_preference(meshfactory, Amanzi::AmanziMesh::MSTK);
  }
    
  TEST (Generate3D)
  {
    auto comm = Amanzi::getDefaultComm();
    bool parallel(comm->NumProc() > 1);
    
    Amanzi::AmanziMesh::FrameworkPreference pref;
    Amanzi::AmanziMesh::MeshFactory meshfactory(comm);

    double x0( 0.0), y0( 0.0), z0( 0.0);
    double x1(10.0), y1(10.0), z1(10.0);
    int nx(10), ny(10), nz(10);

    // The Simple framework is always available, but will only
    // generate in serial

    pref.clear(); pref.push_back(Amanzi::AmanziMesh::Simple);
    meshfactory.set_preference(pref);

    if (parallel) {
      Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh;
      CHECK_THROW(mesh = meshfactory.create(x0, y0, z0,
                                      x1, y1, z1,
                                      nx, ny, nz),
                  Amanzi::AmanziMesh::Message);
    } else {
      Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh;
      mesh = meshfactory.create(x0, y0, z0,
                          x1, y1, z1,
                          nx, ny, nz);
    }

    // The STK, if available, framework will always generate

    if (framework_available(Amanzi::AmanziMesh::STKMESH)) {
      pref.clear(); pref.push_back(Amanzi::AmanziMesh::STKMESH);
      meshfactory.set_preference(pref);
      Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh =
          meshfactory.create(x0, y0, z0,
                             x1, y1, z1,
                             nx, ny, nz);
      CHECK(!mesh.is_null());
    }

    // The MSTK framework, if available, will always generate

    if (framework_available(Amanzi::AmanziMesh::MSTK)) {
      pref.clear(); pref.push_back(Amanzi::AmanziMesh::MSTK);
      meshfactory.set_preference(pref);
      Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh =
          meshfactory.create(x0, y0, z0,
                             x1, y1, z1,
                             nx, ny, nz);
      CHECK(!mesh.is_null());
    }

    // The MOAB framework cannot generate


    if (framework_available(Amanzi::AmanziMesh::MOAB)) {
      pref.clear(); pref.push_back(Amanzi::AmanziMesh::MOAB);
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
    bool parallel(comm->NumProc() > 1);
    
    Amanzi::AmanziMesh::FrameworkPreference pref;
    Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh;
    Amanzi::AmanziMesh::MeshFactory meshfactory(comm);

    double x0( 0.0), y0( 0.0);
    double x1(10.0), y1(10.0);
    int nx(10), ny(10);

    // The MSTK framework, if available, will always generate

    if (framework_available(Amanzi::AmanziMesh::MSTK)) {
      pref.clear(); pref.push_back(Amanzi::AmanziMesh::MSTK);
      meshfactory.set_preference(pref);
      mesh = meshfactory.create(x0, y0,
                          x1, y1,
                          nx, ny);
      CHECK(!mesh.is_null());
      mesh.reset();
    }

    // The Simple framework is always available, but 
    // cannot generate 2D meshes

    pref.clear(); pref.push_back(Amanzi::AmanziMesh::Simple);
    meshfactory.set_preference(pref);

    CHECK_THROW(mesh = meshfactory.create(x0, y0,
                                    x1, y1,
                                    nx, ny),
                Amanzi::AmanziMesh::Message);
    mesh.reset();

    // The STK, even if available cannot generate 2D meshes

    if (framework_available(Amanzi::AmanziMesh::STKMESH)) {
      pref.clear(); pref.push_back(Amanzi::AmanziMesh::STKMESH);
      meshfactory.set_preference(pref);
      CHECK_THROW(mesh = meshfactory.create(x0, y0,
                                      x1, y1,
                                      nx, ny),
                  Amanzi::AmanziMesh::Message);
      mesh.reset();
    }

    // The MOAB framework cannot generate


    if (framework_available(Amanzi::AmanziMesh::MOAB)) {
      pref.clear(); pref.push_back(Amanzi::AmanziMesh::MOAB);
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
    
    Amanzi::AmanziMesh::FrameworkPreference pref;
    Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh;
    Amanzi::AmanziMesh::MeshFactory meshfactory(comm);

    // make a parameter list to try out
    
    Teuchos::ParameterList parameter_list;
    Teuchos::Array<int> ncells = Teuchos::tuple<int>(10,10,10);
    Teuchos::Array<double> low = Teuchos::tuple<double>(0.0, 0.0, 0.0);
    Teuchos::Array<double> high = Teuchos::tuple<double>(1.0, 1.0, 1.0);
    parameter_list.set< Teuchos::Array<int> >("number of cells", ncells);
    parameter_list.set< Teuchos::Array<double> >("domain low coordinate", low);
    parameter_list.set< Teuchos::Array<double> >("domain high coordinate", high);
    
    pref.clear(); pref.push_back(Amanzi::AmanziMesh::Simple);
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

    // The STK, if available, framework will always generate

    if (framework_available(Amanzi::AmanziMesh::STKMESH)) {
      pref.clear(); pref.push_back(Amanzi::AmanziMesh::STKMESH);
      meshfactory.set_preference(pref);
      mesh = meshfactory.create(parameter_list);
      CHECK(!mesh.is_null());
      mesh.reset();
    }

    if (framework_available(Amanzi::AmanziMesh::MSTK)) {
      pref.clear(); pref.push_back(Amanzi::AmanziMesh::MSTK);
      meshfactory.set_preference(pref);
      mesh = meshfactory.create(parameter_list);
      CHECK(!mesh.is_null());
      mesh.reset();
    }

    // Other frameworks can't generate, so they should throw
    pref.clear(); 
    if (framework_available(Amanzi::AmanziMesh::MOAB)) {
      pref.push_back(Amanzi::AmanziMesh::MOAB);
    }
    meshfactory.set_preference(pref);
    CHECK_THROW(mesh = meshfactory.create(parameter_list),
                Amanzi::AmanziMesh::Message);
  }

  TEST (ParameterGenerate2)
  {
    auto comm = Amanzi::getDefaultComm();
    bool parallel(comm->NumProc() > 1);
    
    Amanzi::AmanziMesh::FrameworkPreference pref;
    Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh;
    Amanzi::AmanziMesh::MeshFactory meshfactory(comm);

    // make a parameter list to try out
    
    Teuchos::ParameterList parameter_list;
    Teuchos::Array<int> ncells = Teuchos::tuple<int>(10,10);
    Teuchos::Array<double> low = Teuchos::tuple<double>(0.0, 0.0);
    Teuchos::Array<double> high = Teuchos::tuple<double>(1.0, 1.0);
    parameter_list.set< Teuchos::Array<int> >("number of cells", ncells);
    parameter_list.set< Teuchos::Array<double> >("domain low coordinate", low);
    parameter_list.set< Teuchos::Array<double> >("domain high coordinate", high);
    

    // MSTK, if available, can generate 2D meshes

    if (framework_available(Amanzi::AmanziMesh::MSTK)) {
      pref.clear(); pref.push_back(Amanzi::AmanziMesh::MSTK);
      meshfactory.set_preference(pref);
      mesh = meshfactory.create(parameter_list);
      CHECK(!mesh.is_null());
      mesh.reset();
    }

    // Simple mesh is always available but cannot generate 2D meshes

    pref.clear(); pref.push_back(Amanzi::AmanziMesh::Simple);
    meshfactory.set_preference(pref);

    CHECK_THROW(mesh = meshfactory.create(parameter_list),
                Amanzi::AmanziMesh::Message);
    mesh.reset();

    // The STK, even if available, cannot generate 2D meshes

    if (framework_available(Amanzi::AmanziMesh::STKMESH)) {
      pref.clear(); pref.push_back(Amanzi::AmanziMesh::STKMESH);
      meshfactory.set_preference(pref);
      CHECK_THROW(mesh = meshfactory.create(parameter_list),
                  Amanzi::AmanziMesh::Message);
      mesh.reset();
    }

    // Other frameworks can't generate, so they should throw
    pref.clear(); 
    if (framework_available(Amanzi::AmanziMesh::MOAB)) {
      pref.push_back(Amanzi::AmanziMesh::MOAB);
    }
    meshfactory.set_preference(pref);
    CHECK_THROW(mesh = meshfactory.create(parameter_list),
                Amanzi::AmanziMesh::Message);
  }


  // The Simple framework cannot read anything, even if it exists
  TEST (ReadSimple) {
    auto comm = Amanzi::getDefaultComm();
    bool parallel(comm->NumProc() > 1);
    
    Amanzi::AmanziMesh::FrameworkPreference pref;
    Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh;
    Amanzi::AmanziMesh::MeshFactory meshfactory(comm);

    pref.clear(); pref.push_back(Amanzi::AmanziMesh::Simple);
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
    bool parallel(comm->NumProc() > 1);

    Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh;
    Amanzi::AmanziMesh::MeshFactory meshfactory(comm);

    if (Amanzi::AmanziMesh::framework_available(Amanzi::AmanziMesh::MOAB)) {
      mesh = meshfactory.create(MOAB_TEST_FILE);
      CHECK(!mesh.is_null());
    } else {
      CHECK_THROW(mesh = meshfactory.create(MOAB_TEST_FILE),
                  Amanzi::AmanziMesh::Message);
    }

    // Try it with another framework just for grins

    if (Amanzi::AmanziMesh::framework_available(Amanzi::AmanziMesh::STKMESH)) {
      Amanzi::AmanziMesh::FrameworkPreference pref;
      pref.push_back(Amanzi::AmanziMesh::STKMESH);
      meshfactory.set_preference(pref);
      CHECK_THROW(mesh = meshfactory.create(MOAB_TEST_FILE),
                  Amanzi::AmanziMesh::Message);
    }
  }

  TEST (ReadExodus) 
  {
    auto comm = Amanzi::getDefaultComm();
    bool parallel(comm->NumProc() > 1);

    Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh;
    Amanzi::AmanziMesh::MeshFactory meshfactory(comm);

    bool available =
        (Amanzi::AmanziMesh::framework_available(Amanzi::AmanziMesh::STKMESH) && !parallel) ||
        Amanzi::AmanziMesh::framework_available(Amanzi::AmanziMesh::MSTK);

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
    
    Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh;
    Amanzi::AmanziMesh::MeshFactory meshfactory(comm);
    if ((Amanzi::AmanziMesh::framework_available(Amanzi::AmanziMesh::STKMESH) ||
         Amanzi::AmanziMesh::framework_available(Amanzi::AmanziMesh::MSTK)) && 
        parallel) {
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
    
    Amanzi::AmanziMesh::FrameworkPreference pref;
    Amanzi::AmanziMesh::MeshFactory meshfactory(comm);

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

    bool flatten = true;
    bool extrude = false;

    // Simple mesh CANNOT extract a mesh from another mesh

    pref.clear(); pref.push_back(Amanzi::AmanziMesh::Simple);
    meshfactory.set_preference(pref);

    if (parallel) {
      Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh, newmesh;
      CHECK_THROW(mesh = meshfactory.create(x0, y0, z0,
                                      x1, y1, z1,
                                      nx, ny, nz,gm),
                  Amanzi::AmanziMesh::Message);
    } else {
      Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh, newmesh;
      mesh = meshfactory.create(x0, y0, z0,
              x1, y1, z1,
              nx, ny, nz, gm);

      CHECK(!mesh.is_null());
        CHECK_THROW(newmesh = meshfactory.create(mesh,setnames,
                                         Amanzi::AmanziMesh::FACE,
                                         flatten,extrude),
                    Amanzi::AmanziMesh::Message);
    }

    
    // The STK, if available, framework will always generate

    if (framework_available(Amanzi::AmanziMesh::STKMESH)) {
      Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh, newmesh;
      pref.clear(); pref.push_back(Amanzi::AmanziMesh::STKMESH);
      meshfactory.set_preference(pref);
      mesh = meshfactory.create(x0, y0, z0,
                          x1, y1, z1,
                          nx, ny, nz,gm);
      CHECK(!mesh.is_null());

      CHECK_THROW(newmesh = meshfactory.create(mesh,setnames,
                                         Amanzi::AmanziMesh::FACE,
                                         flatten,extrude),
                  Amanzi::AmanziMesh::Message);
    }

    // The MSTK framework, if available, will always generate
    if (framework_available(Amanzi::AmanziMesh::MSTK)) {
      Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh, newmesh;
      pref.clear(); pref.push_back(Amanzi::AmanziMesh::MSTK);
      meshfactory.set_preference(pref);
      mesh = meshfactory.create(x0, y0, z0,
                          x1, y1, z1,
                          nx, ny, nz,gm);
      CHECK(!mesh.is_null());

      newmesh = meshfactory.create(mesh,setnames,Amanzi::AmanziMesh::FACE,
                             flatten,extrude);
      CHECK(!newmesh.is_null());
    }

    // The MOAB framework cannot generate
    if (framework_available(Amanzi::AmanziMesh::MOAB)) {
      Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh;
      pref.clear(); pref.push_back(Amanzi::AmanziMesh::MOAB);
      meshfactory.set_preference(pref);
      CHECK_THROW(mesh = meshfactory.create(x0, y0, z0,
                                      x1, y1, z1,
                                      nx, ny, nz,gm),
                  Amanzi::AmanziMesh::Message);
      mesh.reset();
    }
  }

    
} // SUITE
