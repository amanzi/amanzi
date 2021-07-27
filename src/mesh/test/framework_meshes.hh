/*
  Copyright 2010-201x held jointly by LANL, ORNL, LBNL, and PNNL.
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

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "GeometricModel.hh"
#include "AmanziComm.hh"

using namespace Amanzi;

inline
Teuchos::RCP<AmanziMesh::Mesh> createFrameworkStructuredUnitSquare(
  const AmanziMesh::Preference& pref, int nx, int ny,
  Comm_ptr_type comm=Teuchos::null,
  const Teuchos::RCP<const AmanziGeometry::GeometricModel>& gm=Teuchos::null,
  const Teuchos::RCP<Teuchos::ParameterList>& plist=Teuchos::null)
{
  if (comm == Teuchos::null) comm = getDefaultComm();
  AmanziMesh::MeshFactory fac(comm, gm, plist);
  fac.set_preference(pref);
  return fac.create(0.0,0.0,1.0,1.0,nx,ny);
}


inline
Teuchos::RCP<AmanziMesh::Mesh> createFrameworkStructuredUnitCube(
  const AmanziMesh::Preference& pref, int nx, int ny, int nz,
  Comm_ptr_type comm=Teuchos::null,
  const Teuchos::RCP<const AmanziGeometry::GeometricModel>& gm=Teuchos::null,
  const Teuchos::RCP<Teuchos::ParameterList>& plist=Teuchos::null)
{
  if (comm == Teuchos::null) comm = getDefaultComm();
  AmanziMesh::MeshFactory fac(comm, gm, plist);
  fac.set_preference(pref);
  return fac.create(0.0,0.0,0.0,1.0,1.0,1.0,nx,ny,nz);
}


inline
Teuchos::RCP<AmanziMesh::Mesh> createFrameworkUnstructured(
  const AmanziMesh::Preference& pref,
  const std::string& filename,
  Comm_ptr_type comm=Teuchos::null,
  const Teuchos::RCP<const AmanziGeometry::GeometricModel>& gm=Teuchos::null,
  const Teuchos::RCP<Teuchos::ParameterList>& plist=Teuchos::null)
{
  if (comm == Teuchos::null) comm = getDefaultComm();
  AmanziMesh::MeshFactory fac(comm, gm, plist);
  fac.set_preference(pref);
  return fac.create(filename);
}



