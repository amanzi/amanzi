/*
This is the flow component of the Amanzi code. 
License: BSD
Authors: Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

#include <iostream>

#include "UnitTest++.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Epetra_SerialComm.h"
#include "Epetra_MpiComm.h"

#include "Mesh.hh"
#include "Mesh_MSTK.hh"
//#include "Mesh_simple.hh"
#include "Darcy_PK.hpp"


using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;
using namespace Amanzi::AmanziFlow;

class DarcyProblem {
 public:
  Epetra_MpiComm* comm;
  Teuchos::RCP<AmanziMesh::Mesh> mesh;
  State *S;
  Teuchos::ParameterList dp_list;
  AmanziFlow::Darcy_PK* DPK;
  int MyPID;

  DarcyProblem() 
  {
    comm = new Epetra_MpiComm(MPI_COMM_WORLD);
    MyPID = comm->MyPID();

    Teuchos::ParameterList parameter_list;
    string xmlFileName = "test/flow_darcy_misc.xml";
    updateParametersFromXmlFile(xmlFileName, &parameter_list);

    // create an SIMPLE mesh framework 
    Teuchos::ParameterList region_list = parameter_list.get<Teuchos::ParameterList>("Regions");
    GeometricModelPtr gm = new GeometricModel(3, region_list, comm);
    //mesh = Teuchos::rcp(new Mesh_simple(0.0,0.0,-0.0, 1.0,1.0,1.0, 4, 4, 4, comm, gm)); 
    mesh = Teuchos::rcp(new Mesh_MSTK("test/hexes.exo", comm, gm)); 

    Teuchos::ParameterList flow_list = parameter_list.get<Teuchos::ParameterList>("Flow");
    dp_list = flow_list.get<Teuchos::ParameterList>("Darcy Problem");

    // create Darcy process kernel
    Teuchos::ParameterList state_list = parameter_list.get<Teuchos::ParameterList>("State");
    S = new State(state_list, mesh);
    S->set_time(0.0);

    Teuchos::RCP<Flow_State> FS = Teuchos::rcp(new Flow_State(*S));
    DPK = new Darcy_PK(dp_list, FS);
  }

  ~DarcyProblem() { delete S; delete DPK; delete comm; }

  void createBClist(
      const char* type, const char* bc_x, Teuchos::Array<std::string>& regions, double value) 
  {
    std::string func_list_name;
    if (type == "pressure") {
      func_list_name = "boundary pressure";
    } else if (type == "static head") {
      func_list_name = "water table elevation";
    } else if (type == "mass flux") {
      func_list_name = "outward mass flux";
    }
    Teuchos::ParameterList& bc_list = dp_list.get<Teuchos::ParameterList>("boundary conditions");
    Teuchos::ParameterList& type_list = bc_list.get<Teuchos::ParameterList>(type);

    Teuchos::ParameterList& bc_sublist = type_list.sublist(bc_x);
    bc_sublist.set("regions", regions);

    Teuchos::ParameterList& bc_sublist_named = bc_sublist.sublist(func_list_name);
    Teuchos::ParameterList& function_list = bc_sublist_named.sublist("function-constant");
    function_list.set("value", value);
  }
  
  double calculatePressureCellError(double p0, AmanziGeometry::Point& pressure_gradient)
  {
    Epetra_Vector& pressure = DPK->flow_state_next()->ref_pressure();

    double error_L2 = 0.0;
    int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
    for (int c=0; c<ncells; c++) {
      const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
      double pressure_exact = p0 + pressure_gradient * xc;
//if (MyPID==0) cout << c << " " << pressure[c] << " exact=" <<  pressure_exact << endl;
      error_L2 += std::pow(pressure[c] - pressure_exact, 2.0);
    }
    return sqrt(error_L2);
  }

  double calculatePressureFaceError(double p0, AmanziGeometry::Point& pressure_gradient)
  {
    Epetra_Vector& solution_faces = DPK->ref_solution_faces();

    double error_L2 = 0.0;
    int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
    for (int f=0; f<nfaces; f++) {
      const AmanziGeometry::Point& xf = mesh->face_centroid(f);
      double pressure_exact = p0 + pressure_gradient * xf;
      error_L2 += std::pow(solution_faces[f] - pressure_exact, 2.0);
    }
    return sqrt(error_L2);
  }

  double calculateDarcyMassFluxError(AmanziGeometry::Point& velocity_exact)
  {
    Epetra_Vector& darcy_mass_flux = DPK->flow_state_next()->ref_darcy_mass_flux();

    double error_L2 = 0.0;
    int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
    for (int f=0; f<nfaces; f++) {
      const AmanziGeometry::Point& normal = mesh->face_normal(f);      
//cout << f << " " << darcy_mass_flux[f] << " exact=" << velocity_exact * normal << endl;
      error_L2 += std::pow(darcy_mass_flux[f] - velocity_exact * normal, 2.0);
    }
    return sqrt(error_L2);
  }

  double calculateDarcyVelocityError(AmanziGeometry::Point& velocity_exact)
  {
    Epetra_Vector& darcy_mass_flux = DPK->flow_state_next()->ref_darcy_mass_flux();
    int dim = mesh->space_dimension();

    Epetra_MultiVector velocity(mesh->cell_map(false), dim);
    DPK->deriveDarcyVelocity(velocity);

    double error_L2 = 0.0;
    int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
    for (int c=0; c<ncells; c++) {
//cout << c << " " << velocity[0][c] << " " << velocity[2][c] << " exact=" << velocity_exact << endl;
      AmanziGeometry::Point velocity_num(velocity[0][c], velocity[1][c], velocity[2][c]);
      error_L2 += L22(velocity_num - velocity_exact);
    }
    return sqrt(error_L2);
  }
};


SUITE(Darcy_PK) {
  TEST_FIXTURE(DarcyProblem, DirichletDirichlet) {
    if (MyPID == 0) std::cout <<"Darcy PK: Dirichlet-Dirichlet" << std::endl;

    double rho = DPK->rho();  // set up analytic solution
    double mu = DPK->mu();

    double p0 = 1.0;
    AmanziGeometry::Point pressure_gradient(0.0, 0.0, -1.0);
    AmanziGeometry::Point velocity(3);
    velocity = -rho * (pressure_gradient - rho * DPK->gravity()) / mu;

    Teuchos::Array<std::string> regions(1);  // modify boundary conditions
    regions[0] = string("Top side");
    createBClist("pressure", "BC 1", regions, 0.0);

    regions[0] = string("Bottom side");
    createBClist("pressure", "BC 2", regions, 1.0);
    DPK->resetParameterList(dp_list);

    DPK->Init();  // setup the problem
    DPK->advance_to_steady_state();

    double error = calculatePressureCellError(p0, pressure_gradient);  // error checks
    CHECK(error < 1.0e-8);
    error = calculatePressureFaceError(p0, pressure_gradient);
    CHECK(error < 1.0e-8);
    error = calculateDarcyMassFluxError(velocity);
    CHECK(error < 1.0e-8);
  }


  TEST_FIXTURE(DarcyProblem, DirichletNeumann)
  {
    if (MyPID == 0) std::cout <<"Darcy PK: Dirichlet-Neumann" << std::endl;

    double rho = DPK->rho();  // set up analytic solution
    double mu = DPK->mu();

    double p0 = 1.0;
    AmanziGeometry::Point pressure_gradient(0.0, 0.0, -1.0);
    AmanziGeometry::Point velocity(3);
    velocity = -rho * (pressure_gradient - rho * DPK->gravity()) / mu;
    double u0 = velocity * AmanziGeometry::Point(0.0, 0.0, 1.0);

    Teuchos::Array<std::string> regions(1);  // modify boundary conditions
    regions[0] = string("Top side");
    createBClist("mass flux", "BC 1", regions, u0);

    regions[0] = string("Bottom side");
    createBClist("pressure", "BC 2", regions, 1.0);
    DPK->resetParameterList(dp_list);

    DPK->Init();  // setup the problem
    DPK->advance_to_steady_state();

    double error = calculatePressureCellError(p0, pressure_gradient);
    CHECK(error < 1.0e-8);
    error = calculatePressureFaceError(p0, pressure_gradient);
    CHECK(error < 1.0e-8);
    error = calculateDarcyMassFluxError(velocity);
    CHECK(error < 1.0e-8);
  }

  TEST_FIXTURE(DarcyProblem, StaticHeadDirichlet) {
    if (MyPID == 0) std::cout <<"Darcy PK: StaticHead-Dirichlet" << std::endl;

    double rho = DPK->rho();  // set up analytic solution
    double mu = DPK->mu();

    double p0 = 2.0;
    AmanziGeometry::Point pressure_gradient(0.0, 0.0, -1.0);
    AmanziGeometry::Point velocity(3);
    velocity = -rho * (pressure_gradient - rho * DPK->gravity()) / mu;

    Teuchos::Array<std::string> regions(1);  // modify boundary conditions
    regions[0] = string("Top side");
    createBClist("pressure", "BC 1", regions, 1.0);

    regions[0] = string("Bottom side");
    createBClist("static head", "BC 2", regions, 0.25);
    DPK->resetParameterList(dp_list);

    DPK->Init();  // setup the problem
    DPK->advance_to_steady_state();

    double error = calculatePressureCellError(p0, pressure_gradient);  // error checks
    CHECK(error < 1.0e-8);
    error = calculatePressureFaceError(p0, pressure_gradient);
    CHECK(error < 1.0e-8);
    error = calculateDarcyMassFluxError(velocity);
    CHECK(error < 1.0e-8);
  }
}

SUITE(Darcy_Velocity) {
  TEST_FIXTURE(DarcyProblem, DarcyVelocity)
  {
    if (MyPID == 0) std::cout <<"Darcy PK: velocity" << std::endl;

    double rho = DPK->rho();  // set up analytic solution
    double mu = DPK->mu();

    double p0 = 2.0;
    AmanziGeometry::Point pressure_gradient(0.0, 0.0, -1.0);
    AmanziGeometry::Point velocity(3);
    velocity = -rho * (pressure_gradient - rho * DPK->gravity()) / mu;

    Teuchos::Array<std::string> regions(1);  // modify boundary conditions
    regions[0] = string("Top side");
    createBClist("pressure", "BC 1", regions, 1.0);

    regions[0] = string("Bottom side");
    createBClist("static head", "BC 2", regions, 0.25);
    DPK->resetParameterList(dp_list);

    DPK->Init();  // setup the problem
    DPK->advance_to_steady_state();

    double error = calculateDarcyVelocityError(velocity);
    CHECK(error < 1.0e-8);
  }
}

