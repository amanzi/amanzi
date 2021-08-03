/*
  Shallow water PK
 
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.
*/

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

// Amanzi
#include "LeastSquare.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "ShallowWater_PK.hh"
#include "OutputXDMF.hh"

// General
#define _USE_MATH_DEFINES
#include "math.h"


/* **************************************************************** */
TEST(SHALLOW_WATER_MANUFACTURED_SOLUTION) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::ShallowWater;
    
  Comm_ptr_type comm = Amanzi::getDefaultComm();
    
  int MyPID = comm -> MyPID();
    
  if (MyPID == 0) {
    std::cout<<"Test: 2D Shallow water"<<std::endl;
  }
    
  // Read parameter list
  std::string xmlFilename = "test/shallow_water_manufactured_solution.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFilename);
    
  // Create a mesh framework
  ParameterList regions_list = plist -> get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, regions_list, *comm));
    
  // Creat a mesh
  bool request_faces = true, request_edges = false;
  MeshFactory meshfactory (comm, gm);
  meshfactory.set_preference(Preference ({Framework::MSTK}));
    
  std::vector<double> dx, Linferror, L1error, L2error;
    
  // Rectangular mesh
  RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 80, 80, request_faces, request_edges);
    
  // Polygonal meshes
//  RCP<Mesh> mesh = meshfactory.create ("test/median15x16.exo");
//  RCP<Mesh> mesh = meshfactory.create ("test/random40.exo");
//  RCP<Mesh> mesh = meshfactory.create ("test/triangular16.exo");
  
  // Create a state
        
  Teuchos::ParameterList state_list = plist -> sublist ("state");
  RCP<State> S = rcp(new State(state_list));
  S->RegisterMesh("surface", mesh);
  S->set_time(0.0);
        
  Teuchos::RCP<TreeVector> soln = Teuchos::rcp (new TreeVector());
          
  Teuchos::ParameterList sw_list = plist->sublist("PKs").sublist("shallow water");

  // Create a shallow water PK
  ShallowWater_PK SWPK(sw_list, plist, S, soln);
  SWPK.Setup(S.ptr());
  S->Setup();
  S->InitializeFields();
  S->InitializeEvaluators();
  SWPK.Initialize(S.ptr());

  S->CheckAllFieldsInitialized();
        
  const Epetra_MultiVector &B = *S->GetFieldData("surface-bathymetry")->ViewComponent("cell");
  const Epetra_MultiVector &Bn = *S->GetFieldData("surface-bathymetry")->ViewComponent("node");
  const Epetra_MultiVector &hh = *S->GetFieldData("surface-ponded_depth")->ViewComponent("cell");
  const Epetra_MultiVector &ht = *S->GetFieldData("surface-total_depth")->ViewComponent("cell");
  const Epetra_MultiVector &vel = *S->GetFieldData("surface-velocity")->ViewComponent("cell");
  const Epetra_MultiVector &q = *S->GetFieldData("surface-discharge")->ViewComponent("cell");
  const Epetra_MultiVector &p = *S->GetFieldData("surface-ponded_pressure")->ViewComponent("cell");
        
  // Create a pid vector
  Epetra_MultiVector pid(B);
        
  for (int c = 0; c < pid.MyLength(); ++c) {
    pid[0][c] = MyPID;
  }
        
  // Create screen io
  auto vo = Teuchos::rcp (new Amanzi::VerboseObject("ShallowWater", *plist));
  WriteStateStatistics (*S, *vo);
        
  // Advance in time
  double t_old(0.0), t_new(0.0), dt;
        
  // Initialize io
  Teuchos::ParameterList iolist;
  std::string fname;
  fname = "SW_sol";
  iolist.get<std::string>("file name base", fname);
  OutputXDMF io(iolist, mesh, true, false);
        
  std::string passwd("state");
        
  int iter = 0;
        
  while (iter < 100 ) {
   
    double t_out = t_new;
            
    Epetra_MultiVector ht_ex(ht);
    Epetra_MultiVector vel_ex(vel);
            
    if (iter % 5 == 0) {
               
      io.InitializeCycle(t_out, iter, "");
                
      io.WriteVector(*hh(0), "depth", AmanziMesh::CELL);
      io.WriteVector(*ht(0), "total_depth", AmanziMesh::CELL);
      io.WriteVector(*vel(0), "vx", AmanziMesh::CELL);
      io.WriteVector(*vel(1), "vy", AmanziMesh::CELL);
      io.WriteVector(*q(0), "qx", AmanziMesh::CELL);
      io.WriteVector(*q(1), "qy", AmanziMesh::CELL);
      io.WriteVector(*B(0), "B", AmanziMesh::CELL);
      io.WriteVector(*p(0), "hyd pressure", AmanziMesh::CELL);
      io.WriteVector(*Bn(0), "B_n", AmanziMesh::NODE);
      io.WriteVector(*pid(0), "pid", AmanziMesh::CELL);
                
      io.FinalizeCycle();
      std::cout<<"time: "<<t_new<<std::endl;
    }
            
    dt = SWPK.get_dt();
            
    if (iter < 10) {
      dt = 0.01 * dt;
    }
            
    t_new = t_old + dt;
            
    SWPK.AdvanceStep(t_old, t_new);
    SWPK.CommitStep(t_old, t_new, S);
            
    t_old = t_new;
    iter += 1;
  } // time loop
        
  if (MyPID == 0) {
    std::cout<<"Time-stepping finished. "<<std::endl;
  }
        
  double t_out  = t_new;
        
  io.InitializeCycle(t_out, iter, "");
        
  io.WriteVector(*hh(0), "depth", AmanziMesh::CELL);
  io.WriteVector(*ht(0), "total_depth", AmanziMesh::CELL);
  io.WriteVector(*vel(0), "vx", AmanziMesh::CELL);
  io.WriteVector(*vel(1), "vy", AmanziMesh::CELL);
  io.WriteVector(*q(0), "qx", AmanziMesh::CELL);
  io.WriteVector(*q(1), "qy", AmanziMesh::CELL);
  io.WriteVector(*B(0), "B", AmanziMesh::CELL);
  io.WriteVector(*Bn(0), "B_n", AmanziMesh::NODE);
  io.WriteVector(*p(0), "hyd pressure", AmanziMesh::CELL);
  io.WriteVector(*pid(0), "pid", AmanziMesh::CELL);
        
  io.FinalizeCycle();
  
}

