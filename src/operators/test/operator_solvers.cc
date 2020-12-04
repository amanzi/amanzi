/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "UnitTest++.h"

// Amanzi
#include "AmanziTypes.hh"
#include "MeshFactory.hh"
#include "Tensor.hh"

// Operators
#include "Analytic00.hh"
#include "Analytic02.hh"
#include "Analytic03b.hh"

#include "OperatorDefs.hh"
#include "PDE_DiffusionFV.hh"
#include "PDE_DiffusionFVwithGravity.hh"
#include "PDE_DiffusionMFD.hh"
#include "PDE_DiffusionMFDwithGravity.hh"

#include "DiffusionFixture.hh"

#include "nvToolsExt.h"

struct TestHarness {};

template<class PDE_Diffusion_type>
void test(const Teuchos::RCP<AnalyticBase>& ana,
          const std::string& pc_type,
          const std::string& bc_type,
          const std::string& mesh_type,
          const std::string& disc_type,
          bool symmetric,
          AmanziMesh::Entity_kind scalar_coef, 
          double tol,
          int niters=1) {
  nvtxRangePushA("Generating mesh");
  auto plist = Teuchos::getParametersFromXmlFile("test/operator_solvers.xml");
  auto comm = Amanzi::getDefaultComm(); 
  auto gm = Teuchos::rcp(new AmanziGeometry::GeometricModel(ana->dimension(),
            plist->sublist("regions"), *comm));
  AmanziMesh::MeshFactory meshfactory(comm, gm);
  Teuchos::RCP<const AmanziMesh::Mesh> mesh = meshfactory.create(-1.0, -1.0, 1.0, 1.0, 1000, 100);  
  nvtxRangePop();

  nvtxRangePushA("Create fix");
  DiffusionFixture fix(ana, mesh_type, mesh);
  nvtxRangePop();

  std::cout << std::endl << std::endl << std::endl
            << "================================================================================" << std::endl
            << "Diffusion Test (np=" << fix.comm->getSize() << "): " << disc_type << ", "
            << ana->name() << ", " << pc_type << ", " << mesh_type << std::endl
            << "--------------------------------------------------------------------------------"
            << std::endl;
  nvtxRangePushA("Discretize fix");
  fix.discretize<PDE_Diffusion_type, AmanziMesh::FACE>(disc_type);
  nvtxRangePop();
  if (scalar_coef != AmanziMesh::Entity_kind::UNKNOWN) fix.scalarCoefficient(scalar_coef);
  if (bc_type == "Dirichlet") {
    fix.setBCsDirichlet();
  } else if (bc_type == "DirichletNeumannBox") {
    fix.setBCsDirichletNeumannBox();
  } else {
    AMANZI_ASSERT(false);
  }      
  nvtxRangePushA("Fix setup");
  fix.setup(pc_type, symmetric);
  nvtxRangePop();
  nvtxRangePushA("Fix go");
  for (int i=0; i!=niters-1; ++i) fix.go(0.0);
  fix.go(tol);
  nvtxRangePop();
  std::cout << "=============================================================================" << std::endl;
}

SUITE(SOLVERS) {

  // ifpack2: ILUT 
  TEST(Analytic00_Linear1_FV_DirichletNeumannBox_Generate2D_ILUT) {
    nvtxRangePushA("ifpack2: ILUT");
    auto ana = Teuchos::rcp(new Analytic00(1, 1.0, 1.0, 0.0));
    test<Operators::PDE_DiffusionFV>(
        ana, "ifpack2: ILUT", "Dirichlet", "Generate2D",
        "fv", true, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12);
    nvtxRangePop();
  }
  // ifpack2: RILUK (KSPILUK)
  TEST(Analytic00_Linear1_FV_DirichletNeumannBox_Generate2D_RILUK) {
    nvtxRangePushA("ifpack2: RILUK");
    auto ana = Teuchos::rcp(new Analytic00(1, 1.0, 1.0, 0.0));
    test<Operators::PDE_DiffusionFV>(
        ana, "ifpack2: RILUK", "Dirichlet", "Generate2D",
        "fv", true, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12);
    nvtxRangePop();
  }

}
