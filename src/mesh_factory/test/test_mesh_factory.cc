// -------------------------------------------------------------
// file: test_mesh_factory.cc
// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created March 18, 2011 by William A. Perkins
// Last Change: Fri Mar 18 15:32:30 2011 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------


static const char* SCCS_ID = "$Id$ Battelle PNL";

#include <iostream>
#include <UnitTest++.h>

#include <Epetra_MpiComm.h>

#include "dbc.hh"
#include "../MeshFactory.hh"


SUITE (MeshFramework)
{
  TEST (Generate)
  {
    Epetra_MpiComm comm(MPI_COMM_WORLD);
    bool parallel(comm.NumProc() > 1);
    
    Mesh::FrameworkPreference pref;
    Teuchos::RCP<Mesh_maps_base> mesh;
    Mesh::MeshFactory mesh_factory(comm);

    double x0( 0.0), y0( 0.0), z0( 0.0);
    double x1(10.0), y1(10.0), z1(10.0);
    int nx(10), ny(10), nz(10);

    pref.clear(); pref.push_back(Mesh::Simple);
    mesh = mesh_factory(x0, y0, z0,
                        x1, y1, z1,
                        nx, ny, nz);
    CHECK(!mesh.is_null());

    pref.clear(); pref.push_back(Mesh::STK);
    mesh = mesh_factory(x0, y0, z0,
                        x1, y1, z1,
                        nx, ny, nz);
    CHECK(!mesh.is_null());
   
  }
}
