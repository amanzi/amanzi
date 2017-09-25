/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
// -------------------------------------------------------------
/**
 * @file   test_deform.cc
 * @author Rao V. Garimella
 * @date   Thu Apr 12, 2012
 * 
 * @brief  
 */
// -------------------------------------------------------------

#include <UnitTest++.h>

#include <mpi.h>
#include <iostream>

#include <Epetra_MpiComm.h>

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "FrameworkTraits.hh"
#include "Geometry.hh"


TEST(MESH_COLUMNS)
{

  Epetra_MpiComm comm_(MPI_COMM_WORLD);
  const int nproc(comm_.NumProc());
  const int me(comm_.MyPID());

  // We are not including MOAB since Mesh_MOAB.cc does not have
  // routines for generating a mesh

  const Amanzi::AmanziMesh::Framework frameworks[] = {  
    Amanzi::AmanziMesh::MSTK
  };
  const char *framework_names[] = {
    "MSTK"
  };
  
  const int numframeworks = sizeof(frameworks)/sizeof(Amanzi::AmanziMesh::Framework);


  Amanzi::AmanziMesh::Framework the_framework;
  for (int i = 0; i < numframeworks; i++) {

    // Set the framework

    the_framework = frameworks[i];

    if (!Amanzi::AmanziMesh::framework_available(the_framework)) continue;

    std::cerr << "Testing columns with " << framework_names[i] << std::endl;


    // Create the mesh

    Amanzi::AmanziMesh::MeshFactory factory(&comm_);
    Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh;

    int ierr = 0;
    int aerr = 0;
    try {
      Amanzi::AmanziMesh::FrameworkPreference prefs(factory.preference());
      prefs.clear(); 
      prefs.push_back(the_framework);

      factory.preference(prefs);
      factory.partitioner(Amanzi::AmanziMesh::Partitioner_type::ZOLTAN_RCB);

      mesh = factory(0.0,0.0,0.0,1.0,1.0,1.0,4,4,4);

    } catch (const Amanzi::AmanziMesh::Message& e) {
      std::cerr << ": mesh error: " << e.what() << std::endl;
      ierr++;
    } catch (const std::exception& e) {
      std::cerr << ": error: " << e.what() << std::endl;
      ierr++;
    }

    comm_.SumAll(&ierr, &aerr, 1);

    CHECK_EQUAL(aerr,0);

    // Explicitly call build columns method
    mesh->build_columns();

    int status;
    int ncells = mesh->num_entities(Amanzi::AmanziMesh::CELL,
                                    Amanzi::AmanziMesh::OWNED);

    double dz = 0.25;  // difference in height between cell centroids
                       // or between node points

    // Identify cell above and cell below 

    for (int c = 0; c < ncells; c++) {
      Amanzi::AmanziGeometry::Point ccen = mesh->cell_centroid(c);
      
      int expcellabove = -1, expcellbelow = -1;
      bool found_above = false, found_below = false;

      if (fabs(ccen[2] - dz/2.0) < 1.e-10) found_below = true;  // bottom layer
      if (fabs(ccen[2] - (1.0-dz/2.0)) < 1.e-10) found_above = true;  // top layer

      std::vector<int> adjcells;
      mesh->cell_get_node_adj_cells(c, Amanzi::AmanziMesh::OWNED, &adjcells);
      int nadjcells = adjcells.size();
      for (int k = 0; k < nadjcells && (!found_above || !found_below); k++) {
        int c2 = adjcells[k];
        if (c == c2) continue;
        
        Amanzi::AmanziGeometry::Point ccen2 = mesh->cell_centroid(c2);

        if (fabs(ccen2[0]-ccen[0]) > 1.0e-10 ||
            fabs(ccen2[1]-ccen[1]) > 1.0e-10) continue;

        if (!found_above && (fabs(ccen2[2]-dz - ccen[2]) < 1.0e-10)) {
          expcellabove = c2;
          found_above = true;
        }
        if (!found_below && (fabs(ccen[2]-dz - ccen2[2]) < 1.0e-10)) {
          expcellbelow = c2;
          found_below = true;
        }
      }
      CHECK(found_below);
      CHECK(found_above);

      CHECK_EQUAL(expcellabove,mesh->cell_get_cell_above(c));
      CHECK_EQUAL(expcellbelow,mesh->cell_get_cell_below(c));
    }

    
    int nnodes = mesh->num_entities(Amanzi::AmanziMesh::CELL,
                                    Amanzi::AmanziMesh::OWNED);

    // Identify node_above and node_below with brute force, n^2 search
    // since there are only 125 nodes

    for (int n = 0; n < nnodes; n++) {
      Amanzi::AmanziGeometry::Point coord;
      mesh->node_get_coordinates(n, &coord);

      int expnodeabove = -1;
      bool found_above = false;

      if (coord[2] == 1.0) found_above = true;  // top surface node

      // Get connected cells of node, and check if one of their nodes
      // qualifies as the node above or node below
      
      std::vector<int> nodecells;
      mesh->node_get_cells(n, Amanzi::AmanziMesh::OWNED, &nodecells);
      int nnodecells = nodecells.size();
      
      for (int k = 0; k < nnodecells && !found_above; k++) {
        int c = nodecells[k];

        std::vector<int> cnodes;
        mesh->cell_get_nodes(c, &cnodes);
        int ncnodes = cnodes.size();
        for (int l = 0; l < ncnodes && !found_above; l++) {
          int n2 = cnodes[l];
          if (n == n2) continue;
          
          Amanzi::AmanziGeometry::Point coord2;
          mesh->node_get_coordinates(n2, &coord2);
          
          if (fabs(coord[0] - coord2[0]) > 1e-10 ||
              fabs(coord[1] - coord2[1]) > 1e-10) continue;

          if (!found_above && (fabs(coord2[2]-dz - coord[2]) < 1e-10)) {
            expnodeabove = n2;
            found_above = true;
          }
        }
      }
      CHECK(found_above);
          
      CHECK_EQUAL(expnodeabove,mesh->node_get_node_above(n));
    }
      
  } // for each framework i

}
