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
#include "Mesh_simple.hh"
#include "Richards_PK.hpp"


using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;
using namespace Amanzi::AmanziFlow;

class RichardsProblem {
 public:
  Epetra_MpiComm* comm;
  Teuchos::RCP<AmanziMesh::Mesh> mesh;
  Teuchos::ParameterList rp_list;
  AmanziFlow::Richards_PK* RPK; 

  RichardsProblem() 
  {
    comm = new Epetra_MpiComm(MPI_COMM_WORLD);

    Teuchos::ParameterList parameter_list;
    string xmlFileName = "test/flow_richards_misc.xml";
    updateParametersFromXmlFile(xmlFileName, &parameter_list);

    // create an SIMPLE mesh framework 
    Teuchos::ParameterList region_list = parameter_list.get<Teuchos::ParameterList>("Regions");
    GeometricModelPtr gm = new GeometricModel(3, region_list);
    mesh = Teuchos::rcp(new Mesh_simple(0.0,0.0,-10.0, 1.0,1.0,0.0, 2, 2, 80, comm, gm)); 

    Teuchos::ParameterList flow_list = parameter_list.get<Teuchos::ParameterList>("Flow");
    rp_list = flow_list.get<Teuchos::ParameterList>("Richards Problem");

    // create Richards process kernel
    Teuchos::ParameterList state_list = parameter_list.get<Teuchos::ParameterList>("State");
    State S(state_list, mesh);

    Teuchos::RCP<Flow_State> FS = Teuchos::rcp(new Flow_State(S));
    RPK = new Richards_PK(rp_list, FS);
    RPK->set_standalone_mode(true);
  }

  ~RichardsProblem() { delete RPK; delete comm; }

  void create_bc_list(
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
    Teuchos::ParameterList& bc_list = rp_list.get<Teuchos::ParameterList>("boundary conditions");
    Teuchos::ParameterList& type_list = bc_list.get<Teuchos::ParameterList>(type);

    Teuchos::ParameterList& bc_sublist = type_list.sublist(bc_x);
    bc_sublist.set("regions", regions);

    Teuchos::ParameterList& bc_sublist_named = bc_sublist.sublist(func_list_name);
    Teuchos::ParameterList& function_list = bc_sublist_named.sublist("function-constant");
    function_list.set("value", value);
  }
  
  double cell_pressure_error(double p0, AmanziGeometry::Point& pressure_gradient)
  {
    Epetra_Vector& solution_cells = RPK->get_solution_cells();

    double error_L2 = 0.0;
    int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
    for (int c=0; c<ncells; c++) {
      const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
      double pressure_exact = p0 + pressure_gradient * xc;
//cout << c << " " << solution_cells[c] << " exact=" <<  pressure_exact << endl;
      error_L2 += std::pow(solution_cells[c] - pressure_exact, 2.0);
    }
    return sqrt(error_L2);
  }

  double face_pressure_error(double p0, AmanziGeometry::Point& pressure_gradient)
  {
    Epetra_Vector& solution_faces = RPK->get_solution_faces();

    double error_L2 = 0.0;
    int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
    for (int f=0; f<nfaces; f++) {
      const AmanziGeometry::Point& xf = mesh->face_centroid(f);
      double pressure_exact = p0 + pressure_gradient * xf;
      error_L2 += std::pow(solution_faces[f] - pressure_exact, 2.0);
    }
    return sqrt(error_L2);
  }

  double darcy_flux_error(AmanziGeometry::Point& velocity_exact)
  {
    Epetra_Vector& darcy_flux = *(RPK->get_FS().get_darcy_flux());

    double error_L2 = 0.0;
    int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
    for (int f=0; f<nfaces; f++) {
      const AmanziGeometry::Point& normal = mesh->face_normal(f);      
//cout << f << " " << darcy_flux[f] << " exact=" << velocity_exact * normal << endl;
      error_L2 += std::pow(darcy_flux[f] - velocity_exact * normal, 2.0);
    }
    return sqrt(error_L2);
  }
};


SUITE(Simple_1D_Flow) {
  TEST_FIXTURE(RichardsProblem, DirichletDirichlet) {
    std::cout <<"Richards 1D: Dirichlet-Dirichlet" << std::endl;

    double rho = RPK->get_rho();  // set up analytic solution
    double mu = RPK->get_mu();

    double p0 = 1.0;
    AmanziGeometry::Point pressure_gradient(0.0, 0.0, -1.0);
    AmanziGeometry::Point velocity(3);
    velocity = -rho * (pressure_gradient - rho * RPK->get_gravity()) / mu;

    Teuchos::Array<std::string> regions(1);  // modify boundary conditions
    regions[0] = string("Top side");
    create_bc_list("pressure", "BC 1", regions, 0.0);

    regions[0] = string("Bottom side");
    create_bc_list("pressure", "BC 2", regions, 0.0);
    RPK->resetParameterList(rp_list);

    RPK->Init();  // setup the problem
    RPK->advance_to_steady_state();
for (int c=0; c<320; c+=4) cout << c << " " << RPK->get_solution_cells()[c] << endl;

    double error = cell_pressure_error(p0, pressure_gradient); // error checks
    CHECK(error < 1.0e-8);
    error = face_pressure_error(p0, pressure_gradient);
    CHECK(error < 1.0e-8);
    error = darcy_flux_error(velocity);
    CHECK(error < 1.0e-8);
  }
}

