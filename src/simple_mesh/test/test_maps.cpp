#include <iostream>
#include "stdlib.h"
#include "math.h"
#include "UnitTest++.h"
#include "../Mesh_maps_simple.hh"
#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include "Epetra_SerialComm.h"

TEST(MAPS) {
  
  using namespace std;

#ifdef HAVE_MPI
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm *comm = new Epetra_SerialComm();
#endif


  Mesh_maps_simple Mm(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1, 2, 1, comm); 

  
  cout << "number of cells = " << Mm.count_entities(Mesh_data::CELL,OWNED) << endl;
  cout << "number of faces = " << Mm.count_entities(Mesh_data::FACE,OWNED) << endl;
  cout << "number of nodes = " << Mm.count_entities(Mesh_data::NODE,OWNED) << endl;

  vector<double> x(24);
  vector<int> nodes(8);
  vector<int> faces(6);
  
  for (int i=0; i<Mm.count_entities(Mesh_data::CELL,OWNED); i++)
    {
      Mm.cell_to_nodes(i, nodes.begin(), nodes.end());
      
      for (int j=0; j<8; j++)
	cout << " " << nodes[j];
      cout << endl << endl;
      
      for (int j=0; j<8; j++) {
	Mm.node_to_coordinates(nodes[j],x.begin(),x.begin()+3);
	for (int k=0; k<3; k++)
	  cout << " " << x[k];
	cout << endl << endl;
      }

      Mm.cell_to_faces(i, faces.begin(), faces.end());
      for (int j=0; j<6; j++) {
	Mm.face_to_coordinates(faces[j],x.begin(), x.begin()+12);
	
	for (int k=0; k<4; k++) 
	  cout << x[3*k] << " " << x[3*k+1] << " " << x[3*k+2] << endl;
	
	cout << endl;
      }

      Mm.cell_to_coordinates(i, x.begin(), x.begin()+24);
      for (int j=0; j<8; j++)
	cout << x[3*j] << " " << x[3*j+1] << " " << x[3*j+2] << endl;
      cout << endl;

    }
  
  vector<int> ss;
  
  for (int is=0; is<Mm.num_sets(Mesh_data::FACE); is++)  {
    ss.resize(Mm.get_set_size(is,Mesh_data::FACE,OWNED));
    Mm.get_set(is, Mesh_data::FACE,OWNED,ss.begin(), ss.end());

    for (int k=0; k<ss.size(); k++) {
      cout << is << " : " << ss[k] << endl;
    }
    cout<< endl;
  }



}

