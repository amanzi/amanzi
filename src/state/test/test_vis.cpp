#include "UnitTest++.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Epetra_MpiComm.h"
#include "Epetra_Vector.h"
#include "State.hpp"
#include "Mesh_STK.hh"
#include "Vis.hpp"


SUITE(VISUALIZATION) {

  TEST(VIZ_DUMP_REQUIRED) {
    
    Teuchos::ParameterList plist;

    plist.set<string>("File Name Base","visdump");

    Teuchos::Array<int> sps(3);
    sps[0] = 0;
    sps[1] = 4;
    sps[2] = 10;
    plist.set<Teuchos::Array<int> >("Start_Period_Stop",sps);
    
    Epetra_MpiComm comm(MPI_COMM_WORLD);
    Amanzi::Vis V(plist, &comm);

    
    // test the cycle stuff, the expected result is in cycles_ and 
    // we store the computed result in cycles
    
    int cycles_[31] = { 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0,
                        0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0  };
    int cycles [31];
    for (int ic = 0; ic<=30; ic++)
      {
	cycles[ic] = V.dump_requested(ic);
      }
    CHECK_ARRAY_EQUAL(cycles_, cycles, 31);
  }



  TEST(DUMP_MESH_AND_DATA) 
  {
    // here we just check that the code does not crash when 
    // the mesh and data files are written


    Teuchos::ParameterList plist;
    plist.set<string>("File Name Base","visdump");
    plist.set<bool>("Enable Gnuplot",true);
    Teuchos::ParameterList& i1_ = plist.sublist("Cycle Data");
    
    i1_.set<int>("Start",0);
    i1_.set<int>("End",10);
    i1_.set<int>("Interval",1);


    Epetra_MpiComm comm(MPI_COMM_WORLD);
    Amanzi::Vis V(plist, &comm);   

    // make a simple mesh
    Teuchos::RCP<Amanzi::AmanziMesh::Mesh_STK> Mesh
      = Teuchos::rcp(new Amanzi::AmanziMesh::Mesh_STK(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 8, 1, 1, &comm));

    V.create_files(*Mesh);

    State S(1, Mesh);
    S.set_cycle(3);

    // create some auxillary data
    Epetra_MultiVector aux(*S.get_total_component_concentration());
    
    std::vector<string> compnames(1);
    std::vector<string> auxnames(1);

    compnames[0] = "comp test";
    auxnames[0] = "aux test";

    V.set_compnames(compnames);
    V.set_auxnames(auxnames);

    V.dump_state(S, &aux);

  }



}
