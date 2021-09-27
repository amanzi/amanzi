/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "Analytic00.hh"
#include "Analytic02.hh"
#include "Analytic03b.hh"

void test(const std::string& prec_solver,
          const std::string& bc_type,
          const std::string& mesh_type, int dim, int nx,
          const std::string& disc_type,
          bool symmetric,
          AmanziMesh::Entity_kind scalar_coef, 
          double tol = 1.0e-12,
          int order = 1,
          const std::string& ana = "00",
          int niters = 1,
          Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::null)
{
  DiffusionFixture fix(plist);
  fix.Init(dim, nx, mesh_type);
  if (fix.get_comm()->getRank() == 0) {
    std::cout 
      << "================================================================================" << std::endl
      << "Diffusion Test (np=" << fix.comm->getSize() << "): " << disc_type 
      << ",  " << prec_solver << ",  mesh=" << mesh_type << std::endl
      << "--------------------------------------------------------------------------------"
      << std::endl;
  }

  if (ana == "00") 
    fix.ana = Teuchos::rcp(new Analytic00(order, 1.0, 1.0, 0.0));
  else if (ana == "02") 
    fix.ana = Teuchos::rcp(new Analytic02(dim));
  else if (ana == "03") 
    fix.ana = Teuchos::rcp(new Analytic03b());

  if (disc_type == "mixed")
    fix.Discretize<Amanzi::Operators::PDE_DiffusionMFD>(disc_type, scalar_coef);
  else if (disc_type == "fv")
    fix.Discretize<Amanzi::Operators::PDE_DiffusionFV>(disc_type, scalar_coef);

  if (bc_type == "Dirichlet") {
    fix.SetBCsDirichlet();
    std::cout<<"Set Dirichlet"<<std::endl;
  } else if (bc_type == "DirichletNeumann") {
    fix.SetBCsDirichletNeumann();
  } else if (bc_type == "DirichletNeumannRobin") {
    fix.SetBCsDirichletNeumannRobin();
  } else {
    AMANZI_ASSERT(false);
  }     
  fix.Setup(prec_solver, symmetric);


  auto start = std::chrono::high_resolution_clock::now();
  if (prec_solver == "mat-vec") {
    fix.MatVec(niters);
  } else {
    // use initial guess as the last iteration
    bool initial_guess = (disc_type == "nlfv") ? true : false;
    for (int i = 0; i < niters - 1; ++i) fix.Go(0.0, initial_guess);
    fix.Go(tol, initial_guess);
  }
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
  std::cout<<"Solver time: "<<double(duration.count() * 1e-6)<<" sec"<<std::endl;

  if (fix.get_comm()->getRank() == 0) {
    std::cout << "================================================================================" 
              << std::endl << std::endl;
  }
}

