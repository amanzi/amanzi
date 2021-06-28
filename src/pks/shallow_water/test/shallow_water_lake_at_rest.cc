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
#include "GMVMesh.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"

#include "ShallowWater_PK.hh"

#include "OutputXDMF.hh"

// General

#define _USE_MATH_DEFINES

#include "math.h"

double H_inf = 0.5; // Lake at rest total height

// ---- ---- ---- ---- ---- ---- ---- ----
// Exact Solution
// ---- ---- ---- ---- ---- ---- ---- ----

void lake_at_rest_exact (double t, double x, double y, double &ht, double &u, double &v)
{
    ht = H_inf;
    u = 0;
    v = 0;
}

void lake_at_rest_exact_field (Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh,
                               Epetra_MultiVector &ht_ex, Epetra_MultiVector &vel_ex, double t)
{
    double x, y, ht, u, v;
    
    int ncells_owned = mesh -> num_entities (Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::OWNED);
    
    for (int c = 0; c < ncells_owned; ++c)
    {
        const Amanzi::AmanziGeometry::Point &xc = mesh -> cell_centroid(c);
        
        x = xc[0]; y = xc[1];
        
        lake_at_rest_exact (t, x, y, ht, u, v);
        ht_ex[0][c] = ht;
        vel_ex[0][c] = u;
        vel_ex[1][c] = v;
    }
}


// ---- ---- ---- ---- ---- ---- ---- ----
// Inital Conditions
// ---- ---- ---- ---- ---- ---- ---- ----

void lake_at_rest_setIC (Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh, Teuchos::RCP<Amanzi::State> &S)
{
    double pi = M_PI;
    
    int ncells_owned = mesh -> num_entities (Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::OWNED);
    std::string passwd = "state";
    
    Epetra_MultiVector &B_c = *S -> GetFieldData ("surface-bathymetry", passwd) -> ViewComponent ("cell");
    Epetra_MultiVector &h_c = *S -> GetFieldData ("surface-ponded_depth", passwd) -> ViewComponent ("cell");
    Epetra_MultiVector &ht_c = *S -> GetFieldData ("surface-total_depth", passwd) -> ViewComponent ("cell");
    Epetra_MultiVector &vel_c = *S -> GetFieldData ("surface-velocity", passwd) -> ViewComponent ("cell");
    Epetra_MultiVector &q_c = *S -> GetFieldData ("surface-discharge", "surface-discharge") -> ViewComponent ("cell");
    
    for (int c = 0; c < ncells_owned; ++c)
    {
        const Amanzi::AmanziGeometry::Point &xc = mesh -> cell_centroid(c);
        
        B_c[0][c] = std::max(0.0, 0.25 - 5 * ((xc[0] - 0.5) * (xc[0] - 0.5) + (xc[1] - 0.5) * (xc[1] - 0.5)));
//        B_c[0][c] = 0.0;
//        B_c[0][c] = 0.8 * std::exp( -5*(xc[0] - 0.9)*(xc[0] - 0.9) - 50*(xc[1] - 0.5)*(xc[1] - 0.5) );
//        if (  (xc[0] - 0.3)*(xc[0] - 0.3) + (xc[1] - 0.3)*(xc[1] - 0.3) < 0.1 * 0.1   ) // Perturb the solution
//        if ( std::abs(xc[0] - 0.1) < 0.05 )
//        {
//            ht_c[0][c] = H_inf + 0.01;
//        }
//        else
//        {
//            ht_c[0][c] = H_inf;
//        }
        ht_c[0][c] = H_inf;
        h_c[0][c] = ht_c[0][c] - B_c[0][c];
        vel_c[0][c] = 0.0;
        vel_c[1][c] = 0.0;
        q_c[0][c] = 0.0;
        q_c[1][c] = 0.0;
    }
    
}

// ---- ---- ---- ---- ---- ---- ---- ----
// Error
// ---- ---- ---- ---- ---- ---- ---- ----

void error (Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh,
            Epetra_MultiVector &ht_ex, Epetra_MultiVector &vel_ex,
            const Epetra_MultiVector &ht, const Epetra_MultiVector &vel,
            double &err_max, double &err_L1, double &hmax)
{
    int ncells_owned = mesh -> num_entities (Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::OWNED);
    
    err_max = 0.0;
    err_L1 = 0.0;
    hmax = 0.0;
    
    for (int c = 0; c < ncells_owned; ++c)
    {
        double tmp = std::abs (ht_ex[0][c] - ht[0][c]);
        err_max = std::max (err_max, tmp);
        err_L1 += tmp * ( mesh -> cell_volume(c) );
        hmax = std::sqrt(mesh -> cell_volume(c));
    }
    
    double err_max_tmp, err_L1_tmp;
    
    mesh -> get_comm() -> MaxAll(&err_max, &err_max_tmp, 1);
    mesh -> get_comm() -> SumAll(&err_L1, &err_L1_tmp, 1);
    
    err_max = err_max_tmp; err_L1 = err_L1_tmp;
    
    std::cout<<"err_max: "<<err_max<<std::endl;
    std::cout<<"err_L1: "<<err_L1<<std::endl;
}

// ---- ---- ---- ---- ---- ---- ---- ----
// Run code
// ---- ---- ---- ---- ---- ---- ---- ----
TEST(SHALLOW_WATER_LAKE_AT_REST)
{
    using namespace Teuchos;
    using namespace Amanzi;
    using namespace Amanzi::AmanziMesh;
    using namespace Amanzi::AmanziGeometry;
    using namespace Amanzi::ShallowWater;
    
    Comm_ptr_type comm = Amanzi::getDefaultComm();
    
    int MyPID = comm -> MyPID();
    
    if (MyPID == 0)
    {
        std::cout<<"Test: 2D Shallow water: Lake at rest"<<std::endl;
    }
    
    // Read parameter list
    
    std::string xmlFilename = "test/shallow_water_lake_at_rest.xml";
    Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile (xmlFilename);
    
    // Create a mesh framework
    
    ParameterList regions_list = plist -> get<Teuchos::ParameterList>("regions");
    auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, regions_list, *comm));
    
    // Creat a mesh ----
    
    bool request_faces = true, request_edges = false;
    MeshFactory meshfactory (comm, gm);
    meshfactory.set_preference (Preference ({Framework::MSTK}));
    
    
    std::vector<double> dx, Linferror, L1error, L2error;
    
    for (int NN = 10; NN <= 40; NN *= 2)
    {
        RCP<Mesh> mesh = meshfactory.create (0.0, 0.0, 1.0, 1.0, NN, NN, request_faces, request_edges);
        
        // Create a state
        
        Teuchos::ParameterList state_list = plist -> sublist ("state");
        RCP<State> S = rcp (new State(state_list));
        S -> RegisterMesh ("surface", mesh);
        S -> set_time(0.0);
        
        Teuchos::RCP<TreeVector> soln = Teuchos::rcp (new TreeVector());
        
        Teuchos::ParameterList pk_tree = plist -> sublist ("PK tree").sublist("shallow water");
        
        // Create a shallow water PK
        
        ShallowWater_PK SWPK (pk_tree, plist, S, soln);
        SWPK.Setup(S.ptr());
        S -> Setup();
        S -> InitializeFields();
        S -> InitializeEvaluators();
        SWPK.Initialize(S.ptr());
        lake_at_rest_setIC (mesh, S);
        
        const Epetra_MultiVector &B = *S -> GetFieldData ("surface-bathymetry") -> ViewComponent ("cell");
        const Epetra_MultiVector &hh = *S -> GetFieldData ("surface-ponded_depth") -> ViewComponent ("cell");
        const Epetra_MultiVector &ht = *S -> GetFieldData ("surface-total_depth") -> ViewComponent ("cell");
        const Epetra_MultiVector &vel = *S -> GetFieldData ("surface-velocity") -> ViewComponent ("cell");
        const Epetra_MultiVector &q = *S -> GetFieldData ("surface-discharge") -> ViewComponent ("cell");
        
        // Create a pid vector
        
        Epetra_MultiVector pid(B);
        
        for (int c = 0; c < pid.MyLength(); ++c)
        {
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
        
        while (t_new < 0.48)
        {
            double t_out = t_new;
            
            Epetra_MultiVector ht_ex(ht);
            Epetra_MultiVector vel_ex(vel);
            
            lake_at_rest_exact_field (mesh, ht_ex, vel_ex, t_out);
            
            if (iter % 5 == 0)
            {
                io.InitializeCycle(t_out, iter, "");
                
                io.WriteVector(*hh(0), "depth", AmanziMesh::CELL);
                io.WriteVector(*ht(0), "total_depth", AmanziMesh::CELL);
                io.WriteVector(*vel(0), "vx", AmanziMesh::CELL);
                io.WriteVector(*vel(1), "vy", AmanziMesh::CELL);
                io.WriteVector(*q(0), "qx", AmanziMesh::CELL);
                io.WriteVector(*q(1), "qy", AmanziMesh::CELL);
                io.WriteVector(*B(0), "B", AmanziMesh::CELL);
                io.WriteVector(*pid(0), "pid", AmanziMesh::CELL);
                
                io.WriteVector(*ht_ex(0), "hh_ex", AmanziMesh::CELL);
                io.WriteVector(*vel_ex(0), "vx_ex", AmanziMesh::CELL);
                io.WriteVector(*vel_ex(1), "vy_ex", AmanziMesh::CELL);
                
                io.FinalizeCycle();
            }
            
            dt = SWPK.get_dt();
            
            if (iter < 10)
            {
                dt = 0.01 * dt;
            }
            
            t_new = t_old + dt;
            
            SWPK.AdvanceStep(t_old, t_new);
            SWPK.CommitStep(t_old, t_new, S);
            
            t_old = t_new;
            
            iter += 1;
            
        } // time loop
        
        if (MyPID == 0)
        {
            std::cout<<"Time-stepping finished. "<<std::endl;
        }
        
        double t_out  = t_new;
        
        Epetra_MultiVector ht_ex(ht);
        Epetra_MultiVector vel_ex(vel);
        
        lake_at_rest_exact_field (mesh, ht_ex, vel_ex, t_out);
        
        double err_max, err_L1, hmax;
        
        error (mesh, ht_ex, vel_ex, ht, vel, err_max, err_L1, hmax);
        
        std::cout<<"hmax: "<<hmax<<std::endl;
        
        dx.push_back (hmax);
        Linferror.push_back (err_max);
        L1error.push_back (err_L1);
        
        io.InitializeCycle(t_out, iter, "");
        
        io.WriteVector(*hh(0), "depth", AmanziMesh::CELL);
        io.WriteVector(*ht(0), "total_depth", AmanziMesh::CELL);
        io.WriteVector(*vel(0), "vx", AmanziMesh::CELL);
        io.WriteVector(*vel(1), "vy", AmanziMesh::CELL);
        io.WriteVector(*q(0), "qx", AmanziMesh::CELL);
        io.WriteVector(*q(1), "qy", AmanziMesh::CELL);
        io.WriteVector(*B(0), "B", AmanziMesh::CELL);
        io.WriteVector(*pid(0), "pid", AmanziMesh::CELL);
        
        io.WriteVector(*ht_ex(0), "hh_ex", AmanziMesh::CELL);
        io.WriteVector(*vel_ex(0), "vx_ex", AmanziMesh::CELL);
        io.WriteVector(*vel_ex(1), "vy_ex", AmanziMesh::CELL);
        
        io.FinalizeCycle();
        
    } // NN loop
    
    
    std::cout<<"* ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- *"<<std::endl;
    
    double L1_order = Amanzi::Utils::bestLSfit (dx, L1error);
    
    std::cout<<"computed order (L_1): "<<L1_order<<std::endl;
    
    double Linf_order = Amanzi::Utils::bestLSfit (dx, Linferror);
    
    std::cout<<"computed order (L_inf): "<<Linf_order<<std::endl;
    
    // CHECK_CLOSE (1.5, L1_order, 0.2);
    // CHECK_CLOSE (1.5, Linf_order, 0.2);
}


























