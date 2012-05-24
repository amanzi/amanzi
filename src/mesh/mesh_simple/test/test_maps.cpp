#include <iostream>
#include "stdlib.h"
#include "math.h"
#include "UnitTest++.h"
#include "../Mesh_simple.hh"
#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include "Epetra_SerialComm.h"
#include "GenerationSpec.hh"

SUITE (MeshSimple) {
TEST(MAPS) {
  
  using namespace std;

#ifdef HAVE_MPI
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm *comm = new Epetra_SerialComm();
#endif


  Amanzi::AmanziMesh::Mesh_simple Mm(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1, 1, 1, comm); 

  double xc[] = { 2.0, 2.0, 2.0 };
  Mm.node_set_coordinates(7,xc);

  cout << "number of cells = " << Mm.num_entities(Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::OWNED) << endl;
  cout << "number of faces = " << Mm.num_entities(Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::OWNED) << endl;
  cout << "number of nodes = " << Mm.num_entities(Amanzi::AmanziMesh::NODE,Amanzi::AmanziMesh::OWNED) << endl;

  vector<double> x(24);
  vector<Amanzi::AmanziMesh::Entity_ID> nodes(8);
  vector<Amanzi::AmanziMesh::Entity_ID> faces(6);
  vector<int> face_dirs(6);
  
  for (Amanzi::AmanziMesh::Entity_ID i=0; i<Mm.num_entities(Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::OWNED); i++)
    {

      cout << "cell_to_nodes" << endl;
      Mm.cell_to_nodes(i, nodes.begin(), nodes.end());
      
      for (int j=0; j<8; j++)
	cout << " " << nodes[j];
      cout << endl << endl;
      
      cout << "node_to_coordinates" << endl;
      for (int j=0; j<8; j++) {
	Mm.node_to_coordinates(nodes[j],x.begin(),x.begin()+3);
	for (int k=0; k<3; k++)
	  cout << " " << x[k];
	cout << endl << endl;
      }

      cout << "cell_get_faces" << endl;
      Mm.cell_get_faces_and_dirs(i, &faces, &face_dirs, true);
      double xx[4][3];
      for (int j=0; j<6; j++) {
	Mm.face_to_coordinates(faces[j],x.begin(), x.begin()+12);
	
	for (int k=0; k<4; k++) {
	  cout << x[3*k] << " " << x[3*k+1] << " " << x[3*k+2] << endl;
	
	  for (int ii=0; ii<3; ii++)
	    xx[k][ii] = x[3*k+ii];
	 


	}

	  double v1[3];
	  double v2[3];
	  
	  v1[0] = xx[0][0] - xx[2][0];
	  v1[1] = xx[0][1] - xx[2][1];
	  v1[2] = xx[0][2] - xx[2][2];

	  v2[0] = xx[1][0] - xx[3][0];
	  v2[1] = xx[1][1] - xx[3][1];
	  v2[2] = xx[1][2] - xx[3][2];

	  double cp = (v1[1]*v2[2] - v1[2]*v2[1]) - (v1[0]*v2[2]-v1[2]*v2[0]) + (v1[0]*v2[1]-v1[1]*v2[0]);  
	  cout << cp/abs(cp) << endl;	
	cout << endl;
      }


      cout << "cell_get_face_dirs" << endl;
      for (int j=0; j<6; j++) {

	cout << face_dirs[j] << " ";
      }
      cout << endl << endl;


      cout << "cell_to_coordinates" << endl;
      Mm.cell_to_coordinates(i, x.begin(), x.begin()+24);
      for (int j=0; j<8; j++)
	cout << x[3*j] << " " << x[3*j+1] << " " << x[3*j+2] << endl;
      cout << endl;

    }
  
}

}
