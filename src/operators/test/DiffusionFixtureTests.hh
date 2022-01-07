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

  if (fix.getComm()->MyPID() == 0) {
    std::cout 
      << "================================================================================" << std::endl
      << "Diffusion Test (np=" << fix.comm->NumProc() << "): " << disc_type 
      << ",  " << prec_solver << ",  mesh=" << mesh_type << std::endl
      << "  bc type=" << bc_type << ", problem=Analytic" << ana << std::endl
      << "--------------------------------------------------------------------------------"
      << std::endl;
  }

  if (ana == "00") 
    fix.ana = Teuchos::rcp(new Analytic00(fix.mesh, order));
  else if (ana == "02") 
    fix.ana = Teuchos::rcp(new Analytic02(fix.mesh));
  else if (ana == "03") 
    fix.ana = Teuchos::rcp(new Analytic03(fix.mesh));

  fix.Discretize(disc_type, scalar_coef);

  if (bc_type == "Dirichlet") {
    fix.SetBCsDirichlet();
  } else if (bc_type == "DirichletNeumann") {
    fix.SetBCsDirichletNeumann();
  } else if (bc_type == "DirichletNeumannRobin") {
    fix.SetBCsDirichletNeumannRobin();
  } else {
    AMANZI_ASSERT(false);
  }      
  fix.Setup(prec_solver, symmetric);

  if (prec_solver == "mat-vec") {
    fix.MatVec(niters);
  } else {
    for (int i = 0; i < niters - 1; ++i) fix.Go(0.0);
    fix.Go(tol);
  }

  if (fix.getComm()->MyPID() == 0) {
    std::cout << "================================================================================" 
              << std::endl << std::endl;
  }
}


void testWGravity(double gravity,
                  const std::string& pc_type,
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

  if (fix.getComm()->MyPID() == 0) {
    std::cout
      << "================================================================================" << std::endl
      << "DiffusionWithGravity Test (np=" << fix.comm->NumProc() << "): "
      << disc_type << ", PC: " << pc_type << ", mesh=" << mesh_type << std::endl
      << "  bc type=" << bc_type << ", problem=Analytic" << ana << std::endl
      << "--------------------------------------------------------------------------------"
      << std::endl;
  }

  if (ana == "00") 
    fix.ana = Teuchos::rcp(new Analytic00(fix.mesh, order, gravity));
  else if (ana == "02") 
    fix.ana = Teuchos::rcp(new Analytic02(fix.mesh, gravity));
  else if (ana == "03") 
    fix.ana = Teuchos::rcp(new Analytic03(fix.mesh));

  fix.DiscretizeWithGravity(disc_type, gravity, scalar_coef);

  if (bc_type == "Dirichlet") {
    fix.SetBCsDirichlet();
  } else if (bc_type == "DirichletNeumann") {
    fix.SetBCsDirichletNeumann();
  } else {
    AMANZI_ASSERT(false);
  }      
  fix.Setup(pc_type, symmetric);
  for (int i = 0; i < niters - 1; ++i) fix.Go(0.0);
  fix.Go(tol);

  if (fix.getComm()->MyPID() == 0) {
    std::cout << "================================================================================" 
              << std::endl << std::endl;
  }
}

