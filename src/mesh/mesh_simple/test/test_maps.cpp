#include <iostream>
#include "stdlib.h"
#include "math.h"
#include "UnitTest++.h"
#include "../Mesh_simple.hh"
#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include "Epetra_SerialComm.h"

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
  Mm.set_coordinate(7,xc,xc+3);

  cout << "number of cells = " << Mm.num_entities(Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::OWNED) << endl;
  cout << "number of faces = " << Mm.num_entities(Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::OWNED) << endl;
  cout << "number of nodes = " << Mm.num_entities(Amanzi::AmanziMesh::NODE,Amanzi::AmanziMesh::OWNED) << endl;

  vector<double> x(24);
  vector<unsigned int> nodes(8);
  vector<unsigned int> faces(6);
  vector<int> face_dirs(6);
  
  for (unsigned int i=0; i<Mm.num_entities(Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::OWNED); i++)
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

      cout << "cell_to_faces" << endl;
      Mm.cell_to_faces(i, faces.begin(), faces.end());
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


      cout << "cell_to_face_dirs" << endl;
      Mm.cell_to_face_dirs(i, face_dirs.begin(), face_dirs.end());
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
  
  vector<unsigned int> ss;
  
  for (int is=0; is<Mm.num_sets(Amanzi::AmanziMesh::FACE); is++)  {
    ss.resize(Mm.get_set_size(is,Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::OWNED));
    Mm.get_set(is, Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::OWNED,ss.begin(), ss.end());

    for (int k=0; k<ss.size(); k++) {
      cout << is << " : " << ss[k] << endl;
    }
    cout<< endl;
  }



}

TEST (ParameterList) {

#ifdef HAVE_MPI
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm *comm = new Epetra_SerialComm();
#endif

    // make a parameter list to try out
    
    Teuchos::ParameterList parameter_list;
    parameter_list.set<int>("Numer of Cells in X", 10);
    parameter_list.set<int>("Numer of Cells in Y", 10);
    parameter_list.set<int>("Numer of Cells in Z", 10);
    
    parameter_list.set<double>("X_Min", 0);
    parameter_list.set<double>("X_Max", 1);
    
    parameter_list.set<double>("Y_Min", 0);
    parameter_list.set<double>("Y_Max", 1);
    
    parameter_list.set<double>("Z_Min", 0);
    parameter_list.set<double>("Z_Max", 1);

    Teuchos::ParameterList sublist1;
    sublist1.set<double>("Z0", 0.1);
    sublist1.set<double>("Z1", 0.3);
    parameter_list.set("Mesh block 1", sublist1);

    Teuchos::ParameterList sublist2;
    sublist2.set<double>("Z0", 0.7);
    sublist2.set<double>("Z1", 1.0);
    parameter_list.set("Mesh block 2", sublist2);

    parameter_list.set<int>("Number of mesh blocks", 2);

    Teuchos::RCP<Amanzi::AmanziMesh::Mesh> 
      mesh(new Amanzi::AmanziMesh::Mesh_simple(parameter_list, comm));
    CHECK(!mesh.is_null());
    CHECK_EQUAL(3, mesh->num_sets(Amanzi::AmanziMesh::CELL));
    mesh.reset();
}
}
