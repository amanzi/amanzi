/*
 Shallow water PK
 
 Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
 Amanzi is released under the three-clause BSD License.
 The terms of use and "as is" disclaimer for this license are
 provided in the top-level COPYRIGHT file.
 
 Author: Svetlana Tokareva (tokareva@lanl.gov)
 */

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

// Amanzi
#include "GMVMesh.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"

#include "ShallowWater_PK.hh"

#include "OutputXDMF.hh"


/* **************************************************************** */
TEST(SHALLOW_WATER_1D) {
    using namespace Teuchos;
    using namespace Amanzi;
    using namespace Amanzi::AmanziMesh;
    using namespace Amanzi::AmanziGeometry;
    using namespace Amanzi::ShallowWater;
    
    Comm_ptr_type comm = Amanzi::getDefaultComm();
    int MyPID = comm->MyPID();
    if (MyPID == 0) std::cout << "Test: 1D shallow water" << std::endl;
    
    // read parameter list
    std::string xmlFileName = "test/shallow_water_1D.xml";
    Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

    /* create a mesh framework */
    Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
    Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2));
    if (MyPID == 0) std::cout << "Geometric model created." << std::endl;

    // create a mesh
    bool request_faces = true, request_edges = true;
    MeshFactory meshfactory(comm,gm);
    meshfactory.set_preference(Preference({Framework::MSTK, Framework::STK}));
    if (MyPID == 0) std::cout << "Mesh factory created." << std::endl;

    RCP<const Mesh> mesh;
    mesh = meshfactory.create(0.0, 0.0, 10.0, 1.0, 100, 10, request_faces,
			   request_edges);
//    mesh = meshfactory.create("median63x64.exo",request_faces,request_edges);
    if (MyPID == 0) std::cout << "Mesh created." << std::endl;

    // create a state
    RCP<State> S = rcp(new State());
    //S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));
    S->RegisterMesh("surface",rcp_const_cast<Mesh>(mesh));
    S->set_time(0.0);
    if (MyPID == 0) std::cout << "State created." << std::endl;
    
    Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
    
    Teuchos::ParameterList pk_tree = plist->sublist("PK tree").sublist("Shallow water");
    
    // create a shallow water PK
    ShallowWater_PK SWPK(pk_tree,plist,S,soln);
    SWPK.Setup(S.ptr());
    S->Setup();
    //SWPK.CreateDefaultState(mesh, 1);
    S->InitializeFields();
    S->InitializeEvaluators();
    SWPK.Initialize(S.ptr());
    if (MyPID == 0) std::cout << "Shallow water PK created." << std::endl;
    
    // create screen io
    auto vo = Teuchos::rcp(new Amanzi::VerboseObject("ShallowWater", *plist));
    S->WriteStatistics(vo);
    
    // advance in time
    double t_old(0.0), t_new(0.0), dt;
//    Teuchos::RCP<Epetra_MultiVector>
//    tcc = S->GetFieldData("total_component_concentration", passwd)->ViewComponent("cell", false);
    
    std::string passwd("state");
    
    Teuchos::RCP<Epetra_MultiVector> h_vec = S->GetFieldData("surface-ponded_depth",passwd)->ViewComponent("cell");
    Teuchos::RCP<Epetra_MultiVector> u_vec = S->GetFieldData("surface-velocity-x",passwd)->ViewComponent("cell");
    Teuchos::RCP<Epetra_MultiVector> v_vec = S->GetFieldData("surface-velocity-y",passwd)->ViewComponent("cell");
    
    int iter = 0;
    bool flag = true;
    
    while (t_new < 0.5) {

        // initialize io
        Teuchos::ParameterList iolist;
        std::string fname;
        fname = "SW_sol_"+std::to_string(iter);
        iolist.get<std::string>("file name base", fname);
        OutputXDMF io(iolist, mesh, true, false);

        // cycle 1, time t
        double t_out = t_new;

        const Epetra_MultiVector& hh = *S->GetFieldData("surface-ponded_depth",passwd)->ViewComponent("cell");
        const Epetra_MultiVector& ht = *S->GetFieldData("surface-total_depth",passwd)->ViewComponent("cell");
        const Epetra_MultiVector& vx = *S->GetFieldData("surface-velocity-x",passwd)->ViewComponent("cell");
        const Epetra_MultiVector& vy = *S->GetFieldData("surface-velocity-y",passwd)->ViewComponent("cell");
        const Epetra_MultiVector& pid = *S->GetFieldData("surface-PID",passwd)->ViewComponent("cell");

        io.InitializeCycle(t_out, 1);
        io.WriteVector(*hh(0), "depth", AmanziMesh::CELL);
        io.WriteVector(*ht(0), "total_depth", AmanziMesh::CELL);
        io.WriteVector(*vx(0), "vx", AmanziMesh::CELL);
        io.WriteVector(*vy(0), "vy", AmanziMesh::CELL);
        io.WriteVector(*pid(0), "pid", AmanziMesh::CELL);
        io.FinalizeCycle();

        dt = SWPK.get_dt();

        t_new = t_old + dt;

//        Teuchos::RCP<Epetra_MultiVector> tmp(v_vec), F(v_vec);

        std::cout << "h_vec.MyLength() = " << (*h_vec).MyLength() << std::endl;


        SWPK.AdvanceStep(t_old, t_new);
        SWPK.CommitStep(t_old, t_new, S);

        t_old = t_new;
        iter++;

    }

    if (MyPID == 0) std::cout << "Time-stepping finished." << std::endl;
    
    std::cout << "MyPID = " << MyPID << ", iter = " << iter << std::endl;


//    if (MyPID == 0) {
//        GMV::open_data_file(*mesh, (std::string)"transport.gmv");
//        GMV::start_data();
//        GMV::write_cell_data(*h_vec, 0, "depth");
//        GMV::close_data_file();
//    }
    
    // initialize io
    Teuchos::ParameterList iolist;
    std::string fname;
    fname = "SW_sol_"+std::to_string(iter);
    iolist.get<std::string>("file name base", fname);
    OutputXDMF io(iolist, mesh, true, false);

    // cycle 1, time t
    double t_out = t_new;

    const Epetra_MultiVector& hh = *S->GetFieldData("surface-ponded_depth",passwd)->ViewComponent("cell");
    const Epetra_MultiVector& ht = *S->GetFieldData("surface-total_depth",passwd)->ViewComponent("cell");
    const Epetra_MultiVector& vx = *S->GetFieldData("surface-velocity-x",passwd)->ViewComponent("cell");
    const Epetra_MultiVector& vy = *S->GetFieldData("surface-velocity-y",passwd)->ViewComponent("cell");
    const Epetra_MultiVector& pid = *S->GetFieldData("surface-PID",passwd)->ViewComponent("cell");

    io.InitializeCycle(t_out, 1);
    io.WriteVector(*hh(0), "depth", AmanziMesh::CELL);
    io.WriteVector(*ht(0), "total_depth", AmanziMesh::CELL);
    io.WriteVector(*vx(0), "vx", AmanziMesh::CELL);
    io.WriteVector(*vy(0), "vy", AmanziMesh::CELL);
    io.WriteVector(*pid(0), "pid", AmanziMesh::CELL);
    io.FinalizeCycle();
    
}
