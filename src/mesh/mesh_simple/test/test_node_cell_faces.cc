#include <iostream>
#include "stdlib.h"
#include "math.h"
#include "UnitTest++.h"
#include "../Mesh_simple.hh"
#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include "Epetra_SerialComm.h"

TEST(NODE_CELL_FACES) {
  
  using namespace std;

#ifdef HAVE_MPI
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm *comm = new Epetra_SerialComm();
#endif

  const unsigned int exp_nnode = 27;

  Amanzi::AmanziMesh::Mesh_simple Mm(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2, 2, 2, comm); 


  for (int i = 0; i < exp_nnode; i++)
    {

      Amanzi::AmanziMesh::Entity_ID node = i;
      
      Amanzi::AmanziMesh::Entity_ID_List cells;

      Mm.node_get_cells(node, Amanzi::AmanziMesh::OWNED, &cells);

      unsigned int ncells = cells.size();

      for (int j = 0; j < ncells; j++)
	{
	  Amanzi::AmanziMesh::Entity_ID cell = cells[j];

	  Amanzi::AmanziMesh::Entity_ID_List faces;

	  Mm.node_get_cell_faces(node, cell, Amanzi::AmanziMesh::OWNED, &faces);

	  // This is a hex mesh. In any given cell, number of faces
	  // connected to a node should be 3

	  CHECK_EQUAL(3,faces.size());

	  for (int k = 0; k < 3; k++) 
	    {
	      
	      Amanzi::AmanziMesh::Entity_ID face = faces[k];
		
	      Amanzi::AmanziMesh::Entity_ID_List fnodes;

	      Mm.face_get_nodes(face, &fnodes);
	      
	      unsigned int nfnodes = fnodes.size();
	      
	      unsigned int found = 0;

	      for (int n = 0; n < nfnodes; n++)
		{
		  if (fnodes[n] == node)
		    {
		      found = 1;
		      break;
		    }
		}
	      
	      CHECK_EQUAL(1,found);
	    }
	  
	}
    }

}

