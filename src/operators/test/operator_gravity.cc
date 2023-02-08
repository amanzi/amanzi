/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

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
#include "MeshFactory.hh"
#include "GMVMesh.hh"
#include "Tensor.hh"
#include "WhetStoneDefs.hh"

// Operators
#include "OperatorDefs.hh"
#include "PDE_Diffusion.hh"
#include "PDE_DiffusionFactory.hh"
#include "UpwindFlux.hh"

#include "HeatConduction.hh"

/* *****************************************************************
* Comparison of gravity models with constant and vector density.
***************************************************************** */
void
RunTestGravity(std::string op_list_name)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "\nTest: check gravity induced rhs" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_diffusion.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  Teuchos::ParameterList plist = xmlreader.getParameters();
  Teuchos::ParameterList op_list = plist.sublist("PK operator").sublist(op_list_name);

  // create a mesh framework
  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  Teuchos::RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 3, 3);

  // create diffusion coefficient
  int ncells = mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces_wghost = mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::ALL);

  const WhetStone::Tensor Kc(2, 1);
  Teuchos::RCP<std::vector<WhetStone::Tensor>> K =
    Teuchos::rcp(new std::vector<WhetStone::Tensor>());

  for (int c = 0; c < ncells; c++) {
    Kc(0, 0) = 1.0 + fabs((mesh->getCellCentroid(c))[0]);
    K->push_back(Kc);
  }

  AmanziGeometry::Point g(0.0, -1.0);

  // create homogeneous boundary data
  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(mesh, AmanziMesh::Entity_kind::FACE, WhetStone::DOF_Type::SCALAR));
  std::vector<int>& bc_model = bc->bc_model();
  std::vector<double>& bc_value = bc->bc_value();

  // create fluid densities
  CompositeVectorSpace cvs;
  cvs.SetMesh(mesh)
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1)
    ->AddComponent("face", AmanziMesh::Entity_kind::FACE, 1);

  double rho(2.0);
  Teuchos::RCP<CompositeVector> rho_cv = Teuchos::rcp(new CompositeVector(cvs));
  rho_cv->PutScalar(2.0);

  // we need flux and dummy solution to populate nonlinear coefficient
  Teuchos::RCP<CompositeVector> flux = Teuchos::rcp(new CompositeVector(cvs));
  Epetra_MultiVector& flx = *flux->ViewComponent("face", true);

  Point velocity(-1.0, 0.0);
  for (int f = 0; f < nfaces_wghost; f++) {
    const Point& normal = mesh->getFaceNormal(f);
    flx[0][f] = velocity * normal;
  }
  CompositeVector u(cvs);

  // create nonlinear coefficient.
  Teuchos::RCP<HeatConduction> knc = Teuchos::rcp(new HeatConduction(mesh));

  // create upwind model
  Teuchos::ParameterList& ulist = plist.sublist("PK operator").sublist("upwind");
  UpwindFlux upwind(mesh);
  upwind.Init(ulist);

  knc->UpdateValues(*flux, bc_model, bc_value); // argument is not used
  upwind.Compute(*flux, u, bc_model, *knc->values());

  // create first diffusion operator using constant density
  Operators::PDE_DiffusionFactory opfactory;
  Teuchos::RCP<PDE_Diffusion> op1 = opfactory.Create(op_list, mesh, bc, rho, g);

  op1->Setup(K, knc->values(), knc->derivatives());
  op1->UpdateMatrices(flux.ptr(), Teuchos::null);

  // create and populate the second operator using vector density
  Teuchos::RCP<PDE_Diffusion> op2 = opfactory.Create(op_list, mesh, bc, rho_cv, g);

  op2->Setup(K, knc->values(), knc->derivatives());
  op2->UpdateMatrices(flux.ptr(), Teuchos::null);

  // check norm of the right-hand sides
  double a1, a2;
  CompositeVector& rhs1 = *op1->global_operator()->rhs();
  CompositeVector& rhs2 = *op2->global_operator()->rhs();

  rhs1.Norm2(&a1);
  rhs2.Norm2(&a2);

  if (MyPID == 0) {
    std::cout << "||rhs1||=" << a1 << std::endl;
    std::cout << "||rhs2||=" << a2 << std::endl;
  }
  CHECK_CLOSE(a1, a2, 1e-12);
}


/* *****************************************************************
* Two tests for MFD and FV methods.
* **************************************************************** */
TEST(OPERATOR_DIFFUSION_GRAVITY_MFD)
{
  RunTestGravity("diffusion operator gravity mfd");
}

TEST(OPERATOR_DIFFUSION_GRAVITY_FV)
{
  RunTestGravity("diffusion operator gravity fv");
}
