/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
// -------------------------------------------------------------
/**
 * @file   test_column_mesh.cc
 * @author Rao V. Garimella
 * @date   Thu Mar 18, 2015
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------

#include <UnitTest++.h>

#include <mpi.h>
#include <iostream>

#include <Epetra_MpiComm.h>

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "FrameworkTraits.hh"
#include "Geometry.hh"
#include "../Column_mesh.hh"


TEST(COLUMN_MESH_3D)
{

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  const int nproc(comm.NumProc());
  const int me(comm.MyPID());

  // We are not including MOAB since Mesh_MOAB.cc does not have
  // routines for generating a mesh

  const Amanzi::AmanziMesh::Framework frameworks[] = {  
    Amanzi::AmanziMesh::MSTK, Amanzi::AmanziMesh::STKMESH, Amanzi::AmanziMesh::Simple
  };
  const char *framework_names[] = {
    "MSTK","STKMesh","Simple"
  };
  
  const int numframeworks = sizeof(frameworks)/sizeof(Amanzi::AmanziMesh::Framework);


  Amanzi::AmanziMesh::Framework the_framework;
  for (int i = 0; i < numframeworks; i++) {

    int nx = 4, ny = 4, nz = 4;
    int lx = 4, ly = 4, lz = 4;
    int dx = 1.0, dy = 1.0, dz = 1.0;

    // Set the framework

    the_framework = frameworks[i];

    if (!Amanzi::AmanziMesh::framework_available(the_framework)) continue;

    std::cerr << "Testing column mesh with " << framework_names[i] << std::endl;


    // Create the mesh

    Amanzi::AmanziMesh::MeshFactory factory(&comm);
    Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh;

    int ierr = 0;
    int aerr = 0;
    try {
      Amanzi::AmanziMesh::FrameworkPreference prefs(factory.preference());
      prefs.clear(); 
      prefs.push_back(the_framework);

      factory.preference(prefs);

      mesh = factory(0.0,0.0,0.0,lx,ly,lz,nx,ny,nz);

    } catch (const Amanzi::AmanziMesh::Message& e) {
      std::cerr << ": mesh error: " << e.what() << std::endl;
      ierr++;
    } catch (const std::exception& e) {
      std::cerr << ": error: " << e.what() << std::endl;
      ierr++;
    }

    comm.SumAll(&ierr, &aerr, 1);

    CHECK_EQUAL(aerr,0);

    // Perturb the nodes above the base layer just a bit
    
    int nnodes = mesh->num_entities(Amanzi::AmanziMesh::NODE,
                                    Amanzi::AmanziMesh::OWNED);
    
    for (int n = 0; n < nnodes; n++) {
      Amanzi::AmanziGeometry::Point xyz(3);
      mesh->node_get_coordinates(n,&xyz);
      xyz[0] += 0.01*xyz[2];
      xyz[1] += 0.01*xyz[2];
      xyz[2] += 0.005*xyz[0]*xyz[1];
      mesh->node_set_coordinates(n,xyz);
    }

    int status;
    int ncolumns = mesh->num_columns();
    CHECK_EQUAL(16,ncolumns);


    // Create a column mesh from one of the columns

    Amanzi::AmanziMesh::Column_mesh colmesh(*mesh,10);

    // Verify column mesh topology

    int ncells = colmesh.num_entities(Amanzi::AmanziMesh::CELL,
                                      Amanzi::AmanziMesh::OWNED);
    CHECK_EQUAL(4,ncells);

    int nfaces = colmesh.num_entities(Amanzi::AmanziMesh::FACE,
                                      Amanzi::AmanziMesh::OWNED);
    CHECK_EQUAL(5,nfaces);

    nnodes = colmesh.num_entities(Amanzi::AmanziMesh::NODE,
                                  Amanzi::AmanziMesh::OWNED);

    CHECK_EQUAL(20,nnodes);

    for (int j = 0; j < ncells; j++) {
      int expcellabove = (j < ncells-1) ? j+1 : -1;
      CHECK_EQUAL(expcellabove, colmesh.cell_get_cell_above(j));
      int expcellbelow = (j > 0) ? j-1 : -1;
      CHECK_EQUAL(expcellbelow, colmesh.cell_get_cell_below(j));

      Amanzi::AmanziMesh::Entity_ID_List cfaces;
      std::vector<int> cfdirs;
      colmesh.cell_get_faces_and_dirs(j,&cfaces,&cfdirs);

      CHECK_EQUAL(2,cfaces.size());
      CHECK_EQUAL(j,cfaces[0]);
      CHECK_EQUAL(1,cfdirs[0]);
      CHECK_EQUAL(j,cfaces[1]);
      CHECK_EQUAL(0,cfdirs[1]);
    }

    for (int j = 0; j < nfaces; j++) {
      Amanzi::AmanziMesh::Entity_ID_List fcells;
      colmesh.face_get_cells(j,Amanzi::AmanziMesh::OWNED,&fcells);
      
      if (j == 0) {
        CHECK_EQUAL(1,fcells.size());
        CHECK_EQUAL(0,fcells[0]);
      }
      else if (j == nfaces-1) {
        CHECK_EQUAL(1,fcells.size());
        CHECK_EQUAL(nfaces-2,fcells[0]);
      }
      else {
        CHECK_EQUAL(2,fcells.size());
        CHECK_EQUAL(j-1,fcells[0]);
        CHECK_EQUAL(j,fcells[1]);
      }
    }


    // Verify column mesh geometry

    // centroid of base face 

    Amanzi::AmanziGeometry::Point fcenbase(3);
    colmesh.face_centroid(0,&fcenbase);

    // Make sure centroids of other faces are stacked up 
    // exactly above that of the base face

    for (int j = 1; j < nfaces; j++) {
      Amanzi::AmanziGeometry::Point fcen(3);
      colmesh.face_centroid(j,&fcen);
      CHECK_EQUAL(fcenbase[0],fcen[0]);
      CHECK_EQUAL(fcenbase[1],fcen[1]);
    }

    // Make sure centroids of cells are stacked up 
    // exactly above that of the base face and that 
    // their z value is exactly between the z-values
    // of their corresponding bottom and top faces

    for (int i = 0; i < ncells; i++) {
      Amanzi::AmanziGeometry::Point ccen(3), fcen0(3), fcen1(3);

      colmesh.cell_centroid(i,&ccen);

      colmesh.face_centroid(i,&fcen0);
      colmesh.face_centroid(i+1,&fcen1);

      CHECK_EQUAL(fcenbase[0],ccen[0]);
      CHECK_EQUAL(fcenbase[1],ccen[1]);

      CHECK_EQUAL((fcen0[2]+fcen1[2])/2.0,ccen[2]);
    }
      
  } // for each framework i

}
