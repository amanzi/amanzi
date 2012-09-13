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
    plist.set<int>("File Name Digits",5);

    Teuchos::ParameterList& clist = plist.sublist("cycle start period stop").sublist("some name");
    Teuchos::Array<int> csps(3);
    csps[0] = 0;
    csps[1] = 4;
    csps[2] = 10;
    clist.set<Teuchos::Array<int> >("start period stop", csps);
    
    Teuchos::ParameterList& tlist = plist.sublist("time start period stop").sublist("some name");
    Teuchos::Array<double> tsps(3);
    tsps[0] = 0.0;
    tsps[1] = 4.0;
    tsps[2] = 10.0;
    tlist.set<Teuchos::Array<double> >("start period stop", tsps);    

    Teuchos::Array<double> times(2);
    times[0] = 1.0;
    times[1] = 3.0;
    plist.set<Teuchos::Array<double> >("times",times);

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

    // test the time sps stuff
    CHECK_EQUAL(true, V.dump_requested(0.0));
    CHECK_EQUAL(true, V.dump_requested(1.0));
    CHECK_EQUAL(true, V.dump_requested(3.0));
    CHECK_EQUAL(true, V.dump_requested(4.0));
    CHECK_EQUAL(true, V.dump_requested(8.0));
    
    CHECK_EQUAL(false, V.dump_requested(0.5));
    CHECK_EQUAL(false, V.dump_requested(1.1));
    CHECK_EQUAL(false, V.dump_requested(3.2));
    CHECK_EQUAL(false, V.dump_requested(3.99));
    CHECK_EQUAL(false, V.dump_requested(10.0));    

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
    
    State S(1, 2, Mesh);
    std::vector<std::string> names(2);
    names.at(0) = "Aoeui";
    names.at(1) = "Snthd";
    S.set_mineral_names(names);
    S.set_cycle(3);

    // create some auxillary data
    Teuchos::RCP<Epetra_MultiVector> aux = S.get_total_component_concentration();
    
    std::vector<string> compnames(1);
    std::vector<string> auxnames(1);

    compnames[0] = "comp test";
    auxnames[0] = "aux test";

    S.set_compnames(compnames);
    
    S.write_vis(V, aux, auxnames, false);

  }



}
