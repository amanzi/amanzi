//#include "mpi.h"
#include "UnitTest++.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Epetra_SerialComm.h"
#include "Epetra_MpiComm.h"

#include "AztecOO.h"

#include "Mesh_maps_simple.hh"
//#include "DiffusionMatrix.hpp"
//#include "DiffusionPrecon.hpp"
#include "DarcyProblem.hpp"
//#include "DarcyMatvec.hpp"

#include "gmv_mesh.hh"

struct problem_setup
{
  Epetra_Comm *comm;
  Teuchos::RCP<Mesh_maps_simple> mesh;
  Teuchos::ParameterList bc_params;
  DarcyProblem *problem;
  AztecOO *solver;
  Epetra_Vector *solution;
  // parameters for analytic pressure function
  double p0, pgrad[3];
  
  enum Side { LEFT, RIGHT, FRONT, BACK, BOTTOM, TOP };
  
  problem_setup()
  {
#ifdef HAVE_MPI
    comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
    comm = new Epetra_SerialComm;
#endif
    
    // Create the mesh.
    mesh = Teuchos::rcp(new Mesh_maps_simple(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 4, 4, 4, comm));
    
    // Define the default BC parameter list: no flow on all sides.
    // Can overwrite the BC on selected sides before creating the problem.
    bc_params.set("number of BCs", 6);
    set_bc(LEFT,   "No Flow");
    set_bc(RIGHT,  "No Flow");
    set_bc(FRONT,  "No Flow");
    set_bc(BACK,   "No Flow");
    set_bc(BOTTOM, "No Flow");
    set_bc(TOP,    "No Flow");
  }
  
  ~problem_setup()
  {
    delete solution;
    delete solver;
    delete problem;
  }  
  
  void set_bc(Side side, std::string type, double value = 0.0) {
    Teuchos::ParameterList *bc;
    switch (side) {
    case LEFT:
      bc = &bc_params.sublist("BC00");
      bc->set("Side set ID", 3);
      break;
    case RIGHT:
      bc = &bc_params.sublist("BC01");
      bc->set("Side set ID", 1);
      break;
    case FRONT:
      bc = &bc_params.sublist("BC02");
      bc->set("Side set ID", 0);
      break;
    case BACK:
      bc = &bc_params.sublist("BC03");
      bc->set("Side set ID", 2);
      break;
    case BOTTOM:
      bc = &bc_params.sublist("BC04");
      bc->set("Side set ID", 4);
      break;
    case TOP:
      bc = &bc_params.sublist("BC05");
      bc->set("Side set ID", 5);
      break;
    }
    bc->set("Type", type);
    bc->set("BC value", value);
  }
  
  void create_problem()
  {
    // Create the Darcy problem parameter list.
    Teuchos::ParameterList pl;
    //Teuchos::ParameterList &precon_pl = pl.sublist("Diffusion Preconditioner");
    //Teuchos::ParameterList &ml_pl = precon_pl.sublist("ML Parameters");
    //ml_pl.set("default values", "SA");

    // Create the flow BCs from the BC parameter list.
    Teuchos::RCP<FlowBC> bc(new FlowBC(bc_params, mesh));
  
    // Create the problem.
    problem = new DarcyProblem(mesh, pl, bc);

    // Set Darcy model defaults; these can be overwritten before solving the problem.
    problem->SetFluidDensity(1.0);
    problem->SetFluidViscosity(1.0);
    problem->SetPermeability(1.0);
    double g[3] = { 0.0 };
    problem->SetGravity(g);
  }
  
  void solve_problem()
  {
    problem->Assemble();

    solver = new AztecOO;
    solver->SetAztecOption(AZ_solver, AZ_cg);
    solver->SetUserOperator(&(problem->Matvec()));
    solver->SetPrecOperator(&(problem->Precon()));

    // Define the RHS.
    Epetra_Vector b(problem->RHS()); // make a copy, Aztec00 will muck with it
    solver->SetRHS(&b);

    // Register the initial solution guess; will be overwritten with the solution.
    solution = new Epetra_Vector(problem->Map());
    solver->SetLHS(solution);

    solver->Iterate(50, 1.0e-10);
    std::cout << "Solver performed " << solver->NumIters() << " iterations." << std::endl
              << "Norm of true residual = " << solver->TrueResidual() << std::endl;
  }
  
  void set_pressure_constants(double p0_, double pgrad_[])
  {
    p0 = p0_;
    for (int i = 0; i < 3; ++i) pgrad[i] = pgrad_[i];
  }
  
  double pressure(double *x)
  {
    return p0 + pgrad[0]*x[0] + pgrad[1]*x[1] + pgrad[2]*x[2];
  }
  
  void cell_pressure_error(double &l2error)
  {
    Epetra_Vector &p_solve = *(problem->CreateCellView(*solution));
    Epetra_Vector  p_error(problem->CellMap());
        
    double xdata[24];
    Epetra_SerialDenseMatrix x(View, xdata, 3, 3, 8);
    
    for (int j = 0; j < p_error.MyLength(); ++j) {
      mesh->cell_to_coordinates((unsigned int) j, xdata, xdata+24);
      // Compute the cell centroid xc.
      double xc[3];
      for (int k = 0; k < 3; ++k) {
        double s = 0.0;
        for (int i = 0; i < 8; ++i) s += x(k,i);
        xc[k] = s / 8.0;
      }
      // Evaluate the pressure error.
      p_error[j] = pressure(xc) - p_solve[j];
    }
  p_error.Norm2(&l2error);
  delete &p_solve;
  }
  
  void face_pressure_error(double &l2error)
  {
    Epetra_Vector &p_solve = *(problem->CreateFaceView(*solution));
    Epetra_Vector  p_error(problem->FaceMap());
        
    double xdata[12];
    Epetra_SerialDenseMatrix x(View, xdata, 3, 3, 4);
    
    for (int j = 0; j < p_error.MyLength(); ++j) {
      mesh->face_to_coordinates((unsigned int) j, xdata, xdata+12);
      // Compute the face centroid xc.
      double xc[3];
      for (int k = 0; k < 3; ++k) {
        double s = 0.0;
        for (int i = 0; i < 4; ++i) s += x(k,i);
        xc[k] = s / 4.0;
      }
      // Evaluate the pressure error.
      p_error[j] = pressure(xc) - p_solve[j];
    }
  p_error.Norm2(&l2error);
  delete &p_solve;
  }
  
};

SUITE(Simple_1D_Flow) {
  
  TEST_FIXTURE(problem_setup, x_p_p)
  {

    // Set non-default BC before create_problem().
    set_bc(LEFT,  "Pressure Constant", 1.0);
    set_bc(RIGHT, "Pressure Constant", 0.0);

    create_problem();

    // Set non-default model parameters before solve_problem().

    solve_problem();

    // Define the analytic pressure solution:
    // p0 = pressure at (0,0,0).
    // pgrad = constant pressure gradient.
    // Domain is [0,1]^3.
    // pgrad = rho * g - (mu/k) * q, where g is the gravity vector and
    //   q is the expected constant Darcy velocity
    double p0 = 1.0;
    double pgrad[3] = {-1.0, 0.0, 0.0};
    set_pressure_constants(p0, pgrad);

    double error;
    cell_pressure_error(error);
    CHECK(error < 1.0e-8);
    
    face_pressure_error(error);
    CHECK(error < 1.0e-8);
  }

  TEST_FIXTURE(problem_setup, xg_p_p)
  {

    // Set non-default BC before create_problem().
    set_bc(LEFT,  "Static Head", 1.0);
    set_bc(RIGHT, "Static Head", 0.0);

    create_problem();

    // Set non-default model parameters before solve_problem().
    double g[3] = {0.0, 0.0, -1.0};
    problem->SetGravity(g);

    solve_problem();

    // Define the analytic pressure solution:
    // p0 = pressure at (0,0,0).
    // pgrad = constant pressure gradient.
    // Domain is [0,1]^3.
    // pgrad = rho * g - (mu/k) * q, where g is the gravity vector and
    //   q is the expected constant Darcy velocity
    double p0 = 1.0;
    double pgrad[3] = {-1.0, 0.0, -1.0};
    set_pressure_constants(p0, pgrad);

    double error;
    cell_pressure_error(error);
    CHECK(error < 1.0e-8);
    
    face_pressure_error(error);
    CHECK(error < 1.0e-8);
  }

  TEST_FIXTURE(problem_setup, x_q_p)
  {

    // Set non-default BC before create_problem().
    set_bc(LEFT,  "Darcy Constant", -1.0);
    set_bc(RIGHT, "Pressure Constant", 0.0);

    create_problem();

    // Set non-default model parameters before solve_problem().

    solve_problem();

    // Define the analytic pressure solution:
    // p0 = pressure at (0,0,0).
    // pgrad = constant pressure gradient.
    // Domain is [0,1]^3.
    // pgrad = rho * g - (mu/k) * q, where g is the gravity vector and
    //   q is the expected constant Darcy velocity
    double p0 = 1.0;
    double pgrad[3] = {-1.0, 0.0, 0.0};
    set_pressure_constants(p0, pgrad);

    double error;
    cell_pressure_error(error);
    CHECK(error < 1.0e-8);
    
    face_pressure_error(error);
    CHECK(error < 1.0e-8);
  }


  TEST_FIXTURE(problem_setup, xg_q_p)
  {

    // Set non-default BC before create_problem().
    set_bc(LEFT,  "Darcy Constant", -1.0);
    set_bc(RIGHT, "Static Head", 0.0);

    create_problem();

    // Set non-default model parameters before solve_problem().
    double g[3] = {0.0, 0.0, -1.0};
    problem->SetGravity(g);

    solve_problem();

    // Define the analytic pressure solution:
    // p0 = pressure at (0,0,0).
    // pgrad = constant pressure gradient.
    // Domain is [0,1]^3.
    // pgrad = rho * g - (mu/k) * q, where g is the gravity vector and
    //   q is the expected constant Darcy velocity
    double p0 = 1.0;
    double pgrad[3] = {-1.0, 0.0, -1.0};
    set_pressure_constants(p0, pgrad);

    double error;
    cell_pressure_error(error);
    CHECK(error < 1.0e-8);
    
    face_pressure_error(error);
    CHECK(error < 1.0e-8);
  }

  
  TEST_FIXTURE(problem_setup, y_p_p)
  {

    // Set non-default BC before create_problem().
    set_bc(BACK,  "Pressure Constant", 1.0);
    set_bc(FRONT, "Pressure Constant", 0.0);

    create_problem();

    // Set non-default model parameters before solve_problem().

    solve_problem();

    // Define the analytic pressure solution:
    // p0 = pressure at (0,0,0).
    // pgrad = constant pressure gradient.
    // Domain is [0,1]^3.
    // pgrad = rho * g - (mu/k) * q, where g is the gravity vector and
    //   q is the expected constant Darcy velocity
    double p0 = 0.0;
    double pgrad[3] = {0.0, 1.0, 0.0};
    set_pressure_constants(p0, pgrad);

    double error;
    cell_pressure_error(error);
    CHECK(error < 1.0e-8);
    
    face_pressure_error(error);
    std::cout << error << std::endl;
    CHECK(error < 1.0e-8);
  }

  TEST_FIXTURE(problem_setup, yg_p_p)
  {

    // Set non-default BC before create_problem().
    set_bc(BACK,  "Static Head", 1.0);
    set_bc(FRONT, "Static Head", 0.0);

    create_problem();

    // Set non-default model parameters before solve_problem().
    double g[3] = {0.0, 0.0, -1.0};
    problem->SetGravity(g);

    solve_problem();

    // Define the analytic pressure solution:
    // p0 = pressure at (0,0,0).
    // pgrad = constant pressure gradient.
    // Domain is [0,1]^3.
    // pgrad = rho * g - (mu/k) * q, where g is the gravity vector and
    //   q is the expected constant Darcy velocity
    double p0 = 0.0;
    double pgrad[3] = {0.0, 1.0, -1.0};
    set_pressure_constants(p0, pgrad);

    double error;
    cell_pressure_error(error);
    CHECK(error < 1.0e-8);
    
    face_pressure_error(error);
    CHECK(error < 1.0e-8);
  }

  TEST_FIXTURE(problem_setup, y_q_p)
  {

    // Set non-default BC before create_problem().
    set_bc(BACK,  "Darcy Constant", -1.0);
    set_bc(FRONT, "Pressure Constant", 0.0);

    create_problem();

    // Set non-default model parameters before solve_problem().

    solve_problem();

    // Define the analytic pressure solution:
    // p0 = pressure at (0,0,0).
    // pgrad = constant pressure gradient.
    // Domain is [0,1]^3.
    // pgrad = rho * g - (mu/k) * q, where g is the gravity vector and
    //   q is the expected constant Darcy velocity
    double p0 = 0.0;
    double pgrad[3] = {0.0, 1.0, 0.0};
    set_pressure_constants(p0, pgrad);

    double error;
    cell_pressure_error(error);
    CHECK(error < 1.0e-8);
    
    face_pressure_error(error);
    CHECK(error < 1.0e-8);
  }


  TEST_FIXTURE(problem_setup, yg_q_p)
  {

    // Set non-default BC before create_problem().
    set_bc(BACK,  "Darcy Constant", -1.0);
    set_bc(FRONT, "Static Head", 0.0);

    create_problem();

    // Set non-default model parameters before solve_problem().
    double g[3] = {0.0, 0.0, -1.0};
    problem->SetGravity(g);

    solve_problem();

    // Define the analytic pressure solution:
    // p0 = pressure at (0,0,0).
    // pgrad = constant pressure gradient.
    // Domain is [0,1]^3.
    // pgrad = rho * g - (mu/k) * q, where g is the gravity vector and
    //   q is the expected constant Darcy velocity
    double p0 = 0.0;
    double pgrad[3] = {0.0, 1.0, -1.0};
    set_pressure_constants(p0, pgrad);

    double error;
    cell_pressure_error(error);
    CHECK(error < 1.0e-8);
    
    face_pressure_error(error);
    CHECK(error < 1.0e-8);
  }
  
  
  TEST_FIXTURE(problem_setup, z_p_p)
  {

    // Set non-default BC before create_problem().
    set_bc(TOP,    "Pressure Constant", 0.0);
    set_bc(BOTTOM, "Pressure Constant", 1.0);

    create_problem();

    // Set non-default model parameters before solve_problem().

    solve_problem();

    // Define the analytic pressure solution:
    // p0 = pressure at (0,0,0).
    // pgrad = constant pressure gradient.
    // Domain is [0,1]^3.
    // pgrad = rho * g - (mu/k) * q, where g is the gravity vector and
    //   q is the expected constant Darcy velocity
    double p0 = 1.0;
    double pgrad[3] = {0.0, 0.0, -1.0};
    set_pressure_constants(p0, pgrad);

    double error;
    cell_pressure_error(error);
    CHECK(error < 1.0e-8);
    
    face_pressure_error(error);
    CHECK(error < 1.0e-8);
  }
  
  
  TEST_FIXTURE(problem_setup, zg_p_p)
  {

    // Set non-default BC before create_problem().
    set_bc(TOP,    "Pressure Constant", 0.0);
    set_bc(BOTTOM, "Pressure Constant", 2.0);

    create_problem();

    // Set non-default model parameters before solve_problem().
    double g[3] = {0.0, 0.0, -1.0};
    problem->SetGravity(g);

    solve_problem();

    // Define the analytic pressure solution:
    // p0 = pressure at (0,0,0).
    // pgrad = constant pressure gradient.
    // Domain is [0,1]^3.
    // pgrad = rho * g - (mu/k) * q, where g is the gravity vector and
    //   q is the expected constant Darcy velocity
    double p0 = 2.0;
    double pgrad[3] = {0.0, 0.0, -2.0};
    set_pressure_constants(p0, pgrad);

    double error;
    cell_pressure_error(error);
    CHECK(error < 1.0e-8);
    
    face_pressure_error(error);
    CHECK(error < 1.0e-8);
  }
  
  
  TEST_FIXTURE(problem_setup, z_q_p)
  {

    // Set non-default BC before create_problem().
    set_bc(TOP,    "Darcy Constant", 1.0);
    set_bc(BOTTOM, "Pressure Constant", 1.0);

    create_problem();

    // Set non-default model parameters before solve_problem().

    solve_problem();

    // Define the analytic pressure solution:
    // p0 = pressure at (0,0,0).
    // pgrad = constant pressure gradient.
    // Domain is [0,1]^3.
    // pgrad = rho * g - (mu/k) * q, where g is the gravity vector and
    //   q is the expected constant Darcy velocity
    double p0 = 1.0;
    double pgrad[3] = {0.0, 0.0, -1.0};
    set_pressure_constants(p0, pgrad);

    double error;
    cell_pressure_error(error);
    CHECK(error < 1.0e-8);
    
    face_pressure_error(error);
    CHECK(error < 1.0e-8);
  }
  
  
  TEST_FIXTURE(problem_setup, zg_q_p)
  {

    // Set non-default BC before create_problem().
    set_bc(TOP,    "Darcy Constant", 1.0);
    set_bc(BOTTOM, "Pressure Constant", 2.0);

    create_problem();

    // Set non-default model parameters before solve_problem().
    double g[3] = {0.0, 0.0, -1.0};
    problem->SetGravity(g);

    solve_problem();

    // Define the analytic pressure solution:
    // p0 = pressure at (0,0,0).
    // pgrad = constant pressure gradient.
    // Domain is [0,1]^3.
    // pgrad = rho * g - (mu/k) * q, where g is the gravity vector and
    //   q is the expected constant Darcy velocity
    double p0 = 2.0;
    double pgrad[3] = {0.0, 0.0, -2.0};
    set_pressure_constants(p0, pgrad);

    double error;
    cell_pressure_error(error);
    CHECK(error < 1.0e-8);
    
    face_pressure_error(error);
    CHECK(error < 1.0e-8);
  }


}
