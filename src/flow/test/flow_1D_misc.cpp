#include "UnitTest++.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Epetra_SerialComm.h"
#include "Epetra_MpiComm.h"

#include "AztecOO.h"

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "DarcyProblem.hpp"
#include "cell_geometry.hh"

#include <iostream>

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;

struct problem_setup
{
  Epetra_MpiComm *comm;
  //Teuchos::RCP<Mesh_simple> mesh;
  Teuchos::RCP<Mesh> mesh;
  GeometricModel *gm;
  Teuchos::ParameterList params;
  Amanzi::DarcyProblem *problem;
  AztecOO *solver;
  Epetra_Vector *solution;
  // parameters for analytic pressure function
  double p0, pgrad[3];
  // other parameters
  double rho, mu, g[3];

  problem_setup()
  {
    comm = new Epetra_MpiComm(MPI_COMM_WORLD);

    // Generate the mesh file name.
    std::string filename;
    if (comm->NumProc() == 1)
      filename = "test/4x4x4.exo";
    else 
      filename = "test/4x4x4.par";
    
    // Create the geometric model.
    Teuchos::ParameterList regions;
    regions.sublist("FRONT").sublist("Region: Labeled Set").
        set("File",filename).set("Format","Exodus II").set("Label","1").set("Entity","Face");
    regions.sublist("RIGHT").sublist("Region: Labeled Set").
        set("File",filename).set("Format","Exodus II").set("Label","2").set("Entity","Face");
    regions.sublist("BACK").sublist("Region: Labeled Set").
        set("File",filename).set("Format","Exodus II").set("Label","3").set("Entity","Face");
    regions.sublist("LEFT").sublist("Region: Labeled Set").
        set("File",filename).set("Format","Exodus II").set("Label","4").set("Entity","Face");
    regions.sublist("BOTTOM").sublist("Region: Labeled Set").
        set("File",filename).set("Format","Exodus II").set("Label","5").set("Entity","Face");
    regions.sublist("TOP").sublist("Region: Labeled Set").
        set("File",filename).set("Format","Exodus II").set("Label","6").set("Entity","Face");
    gm = new GeometricModel(3,regions);

    // Create the mesh.
    MeshFactory mesh_fact(*comm);
    mesh = mesh_fact(filename, gm);

    // Set default model parameters.
    rho = 1.0;
    mu = 1.0;
    for (int i=0; i<3; ++i) g[i] = 0.0;
    
    // Default boundary conditions are no-flux
    params.sublist("boundary conditions");
  }

  ~problem_setup()
  {
    delete solution;
    delete solver;
    delete problem;
  }

  void set_bc(const char *side, const char *type, double value) {
    Teuchos::Array<std::string> reg(1,side);
    std::string func_list_name;
    if (type == "pressure") {
      func_list_name = "boundary pressure";
    } else if (type == "static head") {
      func_list_name = "water table elevation";
    } else if (type == "mass flux") {
      func_list_name = "outward mass flux";
    }
    params.sublist("boundary conditions").sublist(type).sublist(side).set("regions",reg).
        sublist(func_list_name).sublist("function-constant").set("value",value);
  }
  
  void setFluidDensity(double value) { rho = value; }
  void setFluidViscosity(double value) { mu = value; }
  void setGravity(double *value) { for (int i=0; i<3; ++i) g[i] = value[i]; }

  void create_problem()
  {
    // Create the Darcy problem parameter list.
    //Teuchos::ParameterList &precon_pl = pl.sublist("Diffusion Preconditioner");
    //Teuchos::ParameterList &ml_pl = precon_pl.sublist("ML Parameters");
    //ml_pl.set("default values", "SA");
    
    params.set("fluid density", rho);
    params.set("fluid viscosity", mu);
    params.set("gravity", -g[2]);

    // Create the problem.
    problem = new Amanzi::DarcyProblem(mesh, params);

    // Other model parameters; we won't be messing with these.
    problem->SetPermeability(1.0);
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
    //std::cout << "Solver performed " << solver->NumIters() << " iterations." << std::endl
    //          << "Norm of true residual = " << solver->TrueResidual() << std::endl;
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

    double xdata[24], xc[3];
    Epetra_SerialDenseMatrix x(View, xdata, 3, 3, 8);

    for (int j = 0; j < p_error.MyLength(); ++j) {
      mesh->cell_to_coordinates((unsigned int) j, xdata, xdata+24);
      cell_geometry::hex_centroid(x, xc);
      p_error[j] = pressure(xc) - p_solve[j];
    }

    p_error.Norm2(&l2error);
    delete &p_solve;
  }

  void face_pressure_error(double &l2error)
  {
    Epetra_Vector &p_solve = *(problem->CreateFaceView(*solution));
    Epetra_Vector  p_error(problem->FaceMap());

    double xdata[12], xc[3];
    Epetra_SerialDenseMatrix x(View, xdata, 3, 3, 4);

    for (int j = 0; j < p_error.MyLength(); ++j) {
      mesh->face_to_coordinates((unsigned int) j, xdata, xdata+12);
      cell_geometry::quad_face_centroid(x, xc);
      p_error[j] = pressure(xc) - p_solve[j];
    }
  p_error.Norm2(&l2error);
  delete &p_solve;
  }

  void darcy_flux_error(double q[3], double &error1, double &error2)
  {
    Epetra_Vector qflux(problem->FaceMap());
    problem->DeriveDarcyFlux(*solution, qflux, error2);

    double a[3];
    double x[12];
    for (unsigned int j = 0; j < qflux.MyLength(); ++j) {
      mesh->face_to_coordinates(j, x, x+12);
      cell_geometry::quad_face_normal(x, a);
      qflux[j] = qflux[j] - cell_geometry::dot_product(q, a, 3);
    }
    qflux.Norm2(&error1);
  }

  void darcy_velocity_error(double q[3], double &error)
  {
    Epetra_MultiVector qcell(problem->CellMap(),3);
    problem->DeriveDarcyVelocity(*solution, qcell);
    for (int j=0; j<qcell.MyLength(); ++j) {
      for (int k=0; k<3; ++k) qcell[k][j] -= q[k];
    }

    double multi_error[3];
    qcell.Norm2(multi_error);

    error = 0.0;
    for (int k=0; k<3; k++) error += std::pow(multi_error[0], 2.0);
    error = std::sqrt(error);
  }

};

SUITE(Simple_1D_Flow) {

  TEST_FIXTURE(problem_setup, x_p_p)
  {
    std::cout <<"Flow 1D: test 1" << std::endl;
    // Set non-default BC before create_problem().
    set_bc("LEFT",  "pressure", 1.0);
    set_bc("RIGHT", "pressure", 0.0);

    // Set non-default model parameters before create_problem().

    create_problem();
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
    //std::cout << "cell pressure error=" << error << std::endl;

    face_pressure_error(error);
    CHECK(error < 1.0e-8);
    //std::cout << "face pressure error=" << error << std::endl;
  }

  TEST_FIXTURE(problem_setup, xg_p_p)
  {
    std::cout <<"Flow 1D: test 2" << std::endl;

    // Set non-default BC before create_problem().
    set_bc("LEFT",  "static head", 1.0);
    set_bc("RIGHT", "static head", 0.0);

    // Set non-default model parameters before create_problem().
    double g[3] = {0.0, 0.0, -1.0};
    setGravity(g);

    create_problem();
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
    std::cout <<"Flow 1D: test 3" << std::endl;

    // Set non-default BC before create_problem().
    set_bc("LEFT",  "mass flux", -1.0);
    set_bc("RIGHT", "pressure", 0.0);

    // Set non-default model parameters before create_problem().
    
    create_problem();
    solve_problem();

    // Define the analytic pressure solution:
    // p0 = pressure at (0,0,0).
    // pgrad = constant pressure gradient.
    // Domain is [0,1]^3.
    // pgrad = rho * g - (mu/(rho*k)) * f, where g is the gravity
    //   vector and f is the expected constant mass flux
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
    std::cout <<"Flow 1D: test 4" << std::endl;

    // Set non-default BC before create_problem().
    set_bc("LEFT",  "mass flux", -1.0);
    set_bc("RIGHT", "static head", 0.0);

    // Set non-default model parameters before create_problem().
    double g[3] = {0.0, 0.0, -1.0};
    setGravity(g);

    create_problem();
    solve_problem();

    // Define the analytic pressure solution:
    // p0 = pressure at (0,0,0).
    // pgrad = constant pressure gradient.
    // Domain is [0,1]^3.
    // pgrad = rho * g - (mu/(rho*k)) * f, where g is the gravity
    //   vector and f is the expected constant mass flux
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
    std::cout <<"Flow 1D: test 5" << std::endl;

    // Set non-default BC before create_problem().
    set_bc("BACK",  "pressure", 1.0);
    set_bc("FRONT", "pressure", 0.0);

    // Set non-default model parameters before create_problem().

    create_problem();
    solve_problem();

    // Define the analytic pressure solution:
    // p0 = pressure at (0,0,0).
    // pgrad = constant pressure gradient.
    // Domain is [0,1]^3.
    // pgrad = rho * g - (mu/(rho*k)) * f, where g is the gravity
    //   vector and f is the expected constant mass flux
    double p0 = 0.0;
    double pgrad[3] = {0.0, 1.0, 0.0};
    set_pressure_constants(p0, pgrad);

    double error;
    cell_pressure_error(error);
    CHECK(error < 1.0e-8);

    face_pressure_error(error);
    CHECK(error < 1.0e-8);
  }

  TEST_FIXTURE(problem_setup, yg_p_p)
  {
    std::cout <<"Flow 1D: test 6" << std::endl;

    // Set non-default BC before create_problem().
    set_bc("BACK",  "static head", 1.0);
    set_bc("FRONT", "static head", 0.0);

    // Set non-default model parameters before create_problem().
    double g[3] = {0.0, 0.0, -1.0};
    setGravity(g);

    create_problem();
    solve_problem();

    // Define the analytic pressure solution:
    // p0 = pressure at (0,0,0).
    // pgrad = constant pressure gradient.
    // Domain is [0,1]^3.
    // pgrad = rho * g - (mu/(rho*k)) * f, where g is the gravity
    //   vector and f is the expected constant mass flux
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
    std::cout <<"Flow 1D: test 7" << std::endl;

    // Set non-default BC before create_problem().
    set_bc("BACK",  "mass flux", -1.0);
    set_bc("FRONT", "pressure", 0.0);

    // Set non-default model parameters before create_problem().

    create_problem();
    solve_problem();

    // Define the analytic pressure solution:
    // p0 = pressure at (0,0,0).
    // pgrad = constant pressure gradient.
    // Domain is [0,1]^3.
    // pgrad = rho * g - (mu/(rho*k)) * f, where g is the gravity
    //   vector and f is the expected constant mass flux
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
    std::cout <<"Flow 1D: test 8" << std::endl;

    // Set non-default BC before create_problem().
    set_bc("BACK",  "mass flux", -1.0);
    set_bc("FRONT", "static head", 0.0);

    // Set non-default model parameters before create_problem().
    double g[3] = {0.0, 0.0, -1.0};
    setGravity(g);

    create_problem();
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
    std::cout <<"Flow 1D: test 9" << std::endl;

    // Set non-default BC before create_problem().
    set_bc("TOP",    "pressure", 0.0);
    set_bc("BOTTOM", "pressure", 1.0);

    // Set non-default model parameters before create_problem().

    create_problem();
    solve_problem();

    // Define the analytic pressure solution:
    // p0 = pressure at (0,0,0).
    // pgrad = constant pressure gradient.
    // Domain is [0,1]^3.
    // pgrad = rho * g - (mu/(rho*k)) * f, where g is the gravity
    //   vector and f is the expected constant mass flux
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
    std::cout <<"Flow 1D: test 10" << std::endl;

    // Set non-default BC before create_problem().
    set_bc("TOP",    "pressure", 0.0);
    set_bc("BOTTOM", "pressure", 2.0);

    // Set non-default model parameters before create_problem().
    double g[3] = {0.0, 0.0, -1.0};
    setGravity(g);

    create_problem();
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
    std::cout <<"Flow 1D: test 11" << std::endl;

    // Set non-default BC before create_problem().
    set_bc("TOP",    "mass flux", 1.0);
    set_bc("BOTTOM", "pressure", 1.0);

    // Set non-default model parameters before create_problem().

    create_problem();
    solve_problem();

    // Define the analytic pressure solution:
    // p0 = pressure at (0,0,0).
    // pgrad = constant pressure gradient.
    // Domain is [0,1]^3.
    // pgrad = rho * g - (mu/(rho*k)) * f, where g is the gravity
    //   vector and f is the expected constant mass flux
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
    std::cout <<"Flow 1D: test 12" << std::endl;

    // Set non-default BC before create_problem().
    set_bc("TOP",    "mass flux", 1.0);
    set_bc("BOTTOM", "pressure", 2.0);

    // Set non-default model parameters before create_problem().
    double g[3] = {0.0, 0.0, -1.0};
    setGravity(g);

    create_problem();
    solve_problem();

    // Define the analytic pressure solution:
    // p0 = pressure at (0,0,0).
    // pgrad = constant pressure gradient.
    // Domain is [0,1]^3.
    // pgrad = rho * g - (mu/(rho*k)) * f, where g is the gravity
    //   vector and f is the expected constant mass flux
    double p0 = 2.0;
    double pgrad[3] = {0.0, 0.0, -2.0};
    set_pressure_constants(p0, pgrad);

    double error;
    cell_pressure_error(error);
    CHECK(error < 1.0e-8);
    //std::cout << "cell pressure error=" << error << std::endl;

    face_pressure_error(error);
    CHECK(error < 1.0e-8);
    //std::cout << "face pressure error=" << error << std::endl;
  }

}

SUITE(Darcy_Flux) {

  TEST_FIXTURE(problem_setup, Darcy_Flux_X)
  {
    std::cout <<"Darcy flux: test 1" << std::endl;

    // Set non-default BC before create_problem().
    set_bc("LEFT",  "pressure", 1.0);
    set_bc("RIGHT", "pressure", 0.0);

    // Set non-default model parameters before create_problem().

    create_problem();
    solve_problem();

    // Darcy velocity
    double q[3] = { 1.0, 0.0, 0.0 };
    double error1, error2;
    darcy_flux_error(q, error1, error2);
    CHECK(error1 < 1.0e-9); // flux error norm
    CHECK(error2 < 1.0e-9); // flux discrepancy norm
  }


  TEST_FIXTURE(problem_setup, Darcy_Flux_Y)
  {
    std::cout <<"Darcy flux: test 2" << std::endl;

    // Set non-default BC before create_problem().
    set_bc("FRONT", "pressure", 1.0);
    set_bc("BACK",  "pressure", 0.0);

    // Set non-default model parameters before create_problem().

    create_problem();
    solve_problem();

    // Darcy velocity
    double q[3] = { 0.0, 1.0, 0.0 };
    double error1, error2;
    darcy_flux_error(q, error1, error2);
    CHECK(error1 < 1.0e-9); // flux error norm
    CHECK(error2 < 1.0e-9); // flux discrepancy norm
  }


  TEST_FIXTURE(problem_setup, Darcy_Flux_Z)
  {
    std::cout <<"Darcy flux: test 3" << std::endl;

    // Set non-default BC before create_problem().
    set_bc("BOTTOM", "pressure", 1.0);
    set_bc("TOP",    "pressure", 0.0);

    // Set non-default model parameters before solve_problem().

    create_problem();
    solve_problem();

    // Darcy velocity
    double q[3] = { 0.0, 0.0, 1.0 };
    double error1, error2;
    darcy_flux_error(q, error1, error2);
    CHECK(error1 < 1.0e-9); // flux error norm
    CHECK(error2 < 1.0e-9); // flux discrepancy norm
  }
}


SUITE(Darcy_Velocity) {

  TEST_FIXTURE(problem_setup, Darcy_Velocity_X)
  {
    std::cout <<"Darcy velocity: test 1" << std::endl;

    // Set non-default BC before create_problem().
    set_bc("LEFT",  "pressure", 1.0);
    set_bc("RIGHT", "pressure", 0.0);

    // Set non-default model parameters before create_problem().
    setFluidViscosity(2.0);

    create_problem();
    problem->SetPermeability(2.0);
    solve_problem();

    // Darcy velocity
    double q[3] = {1.0, 0.0, 0.0};
    double error;
    darcy_velocity_error(q, error);
    CHECK(error < 1.0e-8);
  }


  TEST_FIXTURE(problem_setup, Darcy_Velocity_Y)
  {
    std::cout <<"Darcy velocity: test 2" << std::endl;

    // Set non-default BC before create_problem().
    set_bc("FRONT", "pressure", 1.0);
    set_bc("BACK",  "pressure", 0.0);

    // Set non-default model parameters before create_problem().
    setFluidViscosity(2.0);

    create_problem();
    problem->SetPermeability(2.0);
    solve_problem();

    // Darcy velocity
    double q[3] = { 0.0, 1.0, 0.0 };
    double error;
    darcy_velocity_error(q, error);
    CHECK(error < 1.0e-8);
  }


  TEST_FIXTURE(problem_setup, Darcy_Velocity_Z)
  {
    std::cout <<"Darcy velocity: test 3" << std::endl;

    // Set non-default BC before create_problem().
    set_bc("BOTTOM", "pressure", 1.0);
    set_bc("TOP",    "pressure", 0.0);

    // Set non-default model parameters before solve_problem().
    setFluidViscosity(2.0);

    create_problem();
    problem->SetPermeability(2.0);
    solve_problem();

    // Darcy velocity
    double q[3] = { 0.0, 0.0, 1.0 };
    double error;
    darcy_velocity_error(q, error);
    CHECK(error < 1.0e-8);
  }


  TEST_FIXTURE(problem_setup, Darcy_Velocity_X_Gravity)
  {
    std::cout <<"Darcy velocity: test 4" << std::endl;

    // Set non-default BC before create_problem().
    set_bc("LEFT",  "static head", 1.0);
    set_bc("RIGHT", "static head", 0.0);

    // Set non-default model parameters before create_problem().
    setFluidViscosity(2.0);
    setFluidDensity(0.5);
    double g[3] = {0.0, 0.0, -2.0};
    setGravity(g);

    create_problem();
    problem->SetPermeability(2.0);
    solve_problem();

    // Darcy velocity
    double q[3] = { 1.0, 0.0, 0.0 };
    double error;
    darcy_velocity_error(q, error);
    CHECK(error < 1.0e-8);
  }


  TEST_FIXTURE(problem_setup, Darcy_Velocity_Y_Gravity)
  {
    std::cout <<"Darcy velocity: test 5" << std::endl;

    // Set non-default BC before create_problem().
    set_bc("FRONT", "static head", 1.0);
    set_bc("BACK",  "static head", 0.0);

    // Set non-default model parameters before create_problem().
    setFluidViscosity(2.0);
    setFluidDensity(0.5);
    double g[3] = {0.0, 0.0, -2.0};
    setGravity(g);

    create_problem();
    problem->SetPermeability(2.0);
    solve_problem();

    // Darcy velocity
    double q[3] = {0.0, 1.0, 0.0};
    double error;
    darcy_velocity_error(q, error);
    CHECK(error < 1.0e-8);
  }


  TEST_FIXTURE(problem_setup, Darcy_Velocity_Z_Gravity)
  {
    std::cout <<"Darcy velocity: test 6" << std::endl;

    // Set non-default BC before create_problem().
    set_bc("BOTTOM", "pressure", 2.0);
    set_bc("TOP",    "pressure", 0.0);

    // Set non-default model parameters before create_problem().
    setFluidViscosity(2.0);
    setFluidDensity(0.5);
    double g[3] = {0.0, 0.0, -2.0};
    setGravity(g);

    create_problem();
    problem->SetPermeability(2.0);
    solve_problem();

    // Darcy velocity
    double q[3] = {0.0, 0.0, 1.0};
    double error;
    darcy_velocity_error(q, error);
    CHECK(error < 1.1e-8);
  }
}

