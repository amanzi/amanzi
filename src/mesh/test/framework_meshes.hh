/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

// Generates Mesh objects for use in testing.

// NOTE: DO NOT CHANGE THESE.  If you need a new mesh, add it instead.  The
// details of these meshes, including topology, geometry, etc, are used by
// various tests, and those tests may fail if you change these meshes.

#pragma once

#include "AmanziComm.hh"
#include "GeometricModel.hh"
#include "MeshFramework.hh"
#include "MeshFrameworkFactory.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"


using namespace Amanzi;

inline
Teuchos::RCP<AmanziMesh::MeshFramework> createFrameworkStructuredUnitQuad(
  const AmanziMesh::Preference& pref, int nx, int ny,
  Comm_ptr_type comm=Teuchos::null,
  const Teuchos::RCP<AmanziGeometry::GeometricModel>& gm=Teuchos::null,
  const Teuchos::RCP<Teuchos::ParameterList>& plist=Teuchos::null,
  double dx=1.0, double dy=1.0)
{
  if (comm == Teuchos::null) comm = getDefaultComm();
  AmanziMesh::MeshFrameworkFactory fac(comm, gm, plist);
  fac.set_preference(pref);
  return fac.create(0.0,0.0, dx,dy, nx,ny);
}

inline
Teuchos::RCP<AmanziMesh::Mesh>  createStructuredUnitQuad(
  const AmanziMesh::Preference& pref, int nx, int ny,
  Comm_ptr_type comm=Teuchos::null,
  const Teuchos::RCP<AmanziGeometry::GeometricModel>& gm=Teuchos::null,
  const Teuchos::RCP<Teuchos::ParameterList>& plist=Teuchos::null,
  double dx=1.0, double dy=1.0)
{
  auto mesh_fw = createFrameworkStructuredUnitQuad(pref, nx, ny, comm, gm, plist, dx, dy);
  auto mesh = Teuchos::rcp(new AmanziMesh::Mesh(mesh_fw, Teuchos::rcp(new AmanziMesh::MeshFrameworkAlgorithms()), Teuchos::null));
  return mesh;
}



inline
Teuchos::RCP<AmanziMesh::MeshFramework> createFrameworkStructuredUnitHex(
  const AmanziMesh::Preference& pref,
  int nx, int ny, int nz,
  Comm_ptr_type comm=Teuchos::null,
  const Teuchos::RCP<AmanziGeometry::GeometricModel>& gm=Teuchos::null,
  const Teuchos::RCP<Teuchos::ParameterList>& plist=Teuchos::null,
  double dx=1.0, double dy=1.0, double dz=1.0)
{
  if (comm == Teuchos::null) comm = getDefaultComm();
  AmanziMesh::MeshFrameworkFactory fac(comm, gm, plist);
  fac.set_preference(pref);
  return fac.create(0.0,0.0,0.0,dx,dy,dz,nx,ny,nz);
}

inline
Teuchos::RCP<AmanziMesh::Mesh> createStructuredUnitHex(
  const AmanziMesh::Preference& pref,
  int nx, int ny, int nz,
  Comm_ptr_type comm=Teuchos::null,
  const Teuchos::RCP<AmanziGeometry::GeometricModel>& gm=Teuchos::null,
  const Teuchos::RCP<Teuchos::ParameterList>& plist=Teuchos::null,
  double dx=1.0, double dy=1.0, double dz=1.0)
{
  auto mesh_fw = createFrameworkStructuredUnitHex(pref, nx,ny,nz, comm, gm, plist, dx,dy,dz);
  auto mesh = Teuchos::rcp(new AmanziMesh::Mesh(mesh_fw, Teuchos::rcp(new AmanziMesh::MeshFrameworkAlgorithms()), Teuchos::null));
  return mesh;
}


inline
Teuchos::RCP<AmanziMesh::MeshFramework> createFrameworkUnstructured(
  const AmanziMesh::Preference& pref,
  const std::string& filename,
  Comm_ptr_type comm=Teuchos::null,
  const Teuchos::RCP<AmanziGeometry::GeometricModel>& gm=Teuchos::null,
  const Teuchos::RCP<Teuchos::ParameterList>& plist=Teuchos::null)
{
  if (comm == Teuchos::null) comm = getDefaultComm();
  AmanziMesh::MeshFrameworkFactory fac(comm, gm, plist);
  fac.set_preference(pref);
  return fac.create(filename);
}

inline
Teuchos::RCP<AmanziMesh::Mesh> createUnstructured(
  const AmanziMesh::Preference& pref,
  const std::string& filename,
  Comm_ptr_type comm=Teuchos::null,
  const Teuchos::RCP<AmanziGeometry::GeometricModel>& gm=Teuchos::null,
  const Teuchos::RCP<Teuchos::ParameterList>& plist=Teuchos::null)
{
  AmanziMesh::MeshFactory fac(comm, gm, plist);
  fac.set_preference(pref);
  auto mesh_fw = createFrameworkUnstructured(pref, filename, comm, gm, plist);
  auto mesh = Teuchos::rcp(new AmanziMesh::Mesh(mesh_fw, Teuchos::rcp(new AmanziMesh::MeshFrameworkAlgorithms()), Teuchos::null));
  return mesh;
}



