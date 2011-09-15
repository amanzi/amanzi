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

    plist.set<string>("file base name","visdump");

    Teuchos::ParameterList& i1 = plist.sublist("Interval 1");
    Teuchos::ParameterList& i1_ = i1.sublist("time range");
    
    i1_.set<double>("start time",0.0);
    i1_.set<double>("end time",10.0);
    i1_.set<double>("time frequency",4.0);

    Teuchos::ParameterList& i2 = plist.sublist("Interval 2");
    Teuchos::ParameterList& i2_ = i2.sublist("time range");
    
    i2_.set<double>("start time",10.0);
    i2_.set<double>("end time",20.0);
    i2_.set<double>("time frequency",3.0);

    
    Teuchos::ParameterList& i3 = plist.sublist("Interval 3");
    Teuchos::ParameterList& i3_ = i3.sublist("cycle range");
    
    i3_.set<int>("start cycle",0);
    i3_.set<int>("end cycle",10);
    i3_.set<int>("cycle frequency",3);

    Teuchos::ParameterList& i4 = plist.sublist("Interval 4");
    Teuchos::ParameterList& i4_ = i4.sublist("cycle range");
    
    i4_.set<int>("start cycle",15);
    i4_.set<int>("end cycle",30);
    i4_.set<int>("cycle frequency",10);
   
    
    Epetra_MpiComm comm(MPI_COMM_WORLD);
    Amanzi::Vis V(plist, &comm);

    
    // test the cycle stuff, the expected result is in cycles_ and 
    // we store the computed result in cycles
    
    int cycles_[31] = { 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0,
                        0, 0, 0, 0,
                        1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0  };
    int cycles [31];
    for (int ic = 0; ic<=30; ic++)
      {
	cycles[ic] = V.dump_requested(0.0,0.0,ic);
      }
    CHECK_ARRAY_EQUAL(cycles_, cycles, 31);

    
    // now we test the times stuff

    CHECK_EQUAL(1, V.dump_requested(-0.1,0.0,-1));
    CHECK_EQUAL(1, V.dump_requested(3.9,4.0,-1));
    CHECK_EQUAL(1, V.dump_requested(7.9,8.1,-1));
    CHECK_EQUAL(1, V.dump_requested(9.9,10.1,-1));
    CHECK_EQUAL(1, V.dump_requested(12.9,13.0,-1));
    CHECK_EQUAL(1, V.dump_requested(15.9,16.9,-1));
    CHECK_EQUAL(1, V.dump_requested(18.9,19.9999,-1));
    
    
    CHECK_EQUAL(0, V.dump_requested(18.9,18.9999,-1));
    CHECK_EQUAL(0, V.dump_requested(0.0,21.0,-1)); 
  }



  TEST(DUMP_MESH_AND_DATA) 
  {
    // here we just check that the code does not crash when 
    // the mesh and data files are written


    Teuchos::ParameterList plist;
    plist.set<string>("file base name","visdump");
    plist.set<bool>("enable gnuplot",true);
    Teuchos::ParameterList& i1 = plist.sublist("Interval 1");
    Teuchos::ParameterList& i1_ = i1.sublist("cycle range");
    
    i1_.set<int>("start cycle",0);
    i1_.set<int>("end cycle",10);
    i1_.set<int>("cycle frequency",1);


    Epetra_MpiComm comm(MPI_COMM_WORLD);
    Amanzi::Vis V(plist, &comm);   

    // make a simple mesh
    Teuchos::RCP<Amanzi::AmanziMesh::Mesh_STK> Mesh
      = Teuchos::rcp(new Amanzi::AmanziMesh::Mesh_STK(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 8, 1, 1, &comm));

    V.create_files(*Mesh);

    State S(1, Mesh);

    // create some auxillary data
    Epetra_MultiVector aux(*S.get_total_component_concentration());
    
    std::vector<string> compnames(1);
    std::vector<string> auxnames(1);

    compnames[0] = "comp test";
    auxnames[0] = "aux test";

    V.set_compnames(compnames);
    V.set_auxnames(auxnames);

    V.dump_state(0.0, 1.0, 0, S, &aux);

  }



}
