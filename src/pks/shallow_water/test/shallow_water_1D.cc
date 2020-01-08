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

    /* create a mesh framework */
    Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
    Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2));
    if (MyPID == 0) std::cout << "Geometric model created." << std::endl;

    // create a mesh
    MeshFactory meshfactory(comm,gm);
    meshfactory.set_preference(Preference({Framework::MSTK, Framework::STK}));
    if (MyPID == 0) std::cout << "Mesh factory created." << std::endl;

    RCP<const Mesh> mesh;
    mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 10, 10);
    if (MyPID == 0) std::cout << "Mesh created." << std::endl;
    
    // create a state
    RCP<State> S = rcp(new State());
    S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));
    S->set_time(0.0);
    S->set_intermediate_time(0.0);
    if (MyPID == 0) std::cout << "State created." << std::endl;
    
    // create a shallow water PK
    ShallowWater_PK SWPK;
    SWPK.Setup(S.ptr());
    //SWPK.CreateDefaultState(mesh, 1);
    S->InitializeFields();
    S->InitializeEvaluators();
    if (MyPID == 0) std::cout << "Shallow water PK created." << std::endl;
    
    // advance in time
    double t_old(0.0), t_new(0.0), dt;
//    Teuchos::RCP<Epetra_MultiVector>
//    tcc = S->GetFieldData("total_component_concentration", passwd)->ViewComponent("cell", false);
    
    int iter = 0;
    bool flag = true;
    while (t_new < 0.25) {
        dt = SWPK.get_dt();
        t_new = t_old + dt;
        
        SWPK.AdvanceStep(t_old, t_new);
        SWPK.CommitStep(t_old, t_new, S);
        
        t_old = t_new;
        iter++;
    }
    if (MyPID == 0) std::cout << "Time-stepping finished." << std::endl;
    
}
