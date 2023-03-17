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
#include "Epetra_MultiVector.h"
#include "EpetraExt_RowMatrixOut.h"
#include "UnitTest++.h"

// Amanzi
#include "MeshFactory.hh"
#include "Mesh_MSTK.hh"
#include "Tensor.hh"

// Operators
#include "Analytic06.hh"

#include "OperatorDefs.hh"
#include "PDE_DiffusionFactory.hh"


/* *****************************************************************
* Exactness test for mixed diffusion solver.
***************************************************************** */
std::pair<double, double>
RunForwardProblem(const std::string& discretization, int nx, int ny)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();

  // create a mesh
  auto mesh_mstk = Teuchos::rcp(new Mesh_MSTK(0., 0., 1., 1., nx, ny, comm));
  auto mesh = Teuchos::rcp(new Mesh(mesh_mstk, Teuchos::rcp(new Amanzi::AmanziMesh::MeshFrameworkAlgorithms()), Teuchos::null)); 

  // modify diffusion coefficient
  Teuchos::RCP<std::vector<WhetStone::Tensor>> K =
    Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  int ncells = mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  int nfaces = mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
  int nfaces_wghost = mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::ALL);

  Analytic06 ana(mesh);

  for (int c = 0; c < ncells; c++) {
    const Point& xc = mesh->getCellCentroid(c);
    const WhetStone::Tensor& Kc = ana.TensorDiffusivity(xc, 0.0);
    K->push_back(Kc);
  }

  // create boundary data
  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(mesh, AmanziMesh::Entity_kind::FACE, WhetStone::DOF_Type::SCALAR));
  std::vector<int>& bc_model = bc->bc_model();
  std::vector<double>& bc_value = bc->bc_value();

  bool flag;
  for (int f = 0; f < nfaces_wghost; f++) {
    const Point& xf = mesh->getFaceCentroid(f);
    double area = mesh->getFaceArea(f);
    Point normal = ana.face_normal_exterior(f, &flag);

    if (fabs(xf[0]) < 1e-6) {
      bc_model[f] = Operators::OPERATOR_BC_NEUMANN;
      bc_value[f] = ana.velocity_exact(xf, 0.0) * normal / area;
    } else if (fabs(xf[1]) < 1e-6) {
      bc_model[f] = Operators::OPERATOR_BC_NEUMANN;
      bc_value[f] = ana.velocity_exact(xf, 0.0) * normal / area;
    } else if (fabs(xf[0] - 1.0) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6) {
      bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
      bc_value[f] = ana.pressure_exact(xf, 0.0);
    }
  }

  // create diffusion operator
  ParameterList op_list;
  op_list.set("discretization primary", discretization);
  op_list.set("nonlinear coefficient", "none");

  Operators::PDE_DiffusionFactory fac;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_c = mesh;
  Teuchos::RCP<PDE_Diffusion> op = fac.Create(op_list, mesh_c);
  op->SetBCs(bc, bc);
  op->SetTensorCoefficient(K);
  op->SetScalarCoefficient(Teuchos::null, Teuchos::null);
  op->UpdateMatrices(Teuchos::null, Teuchos::null);

  // create the rhs
  CompositeVector source(*op->global_operator()->rhs());
  Epetra_MultiVector& source_c = *source.ViewComponent("cell", false);
  for (int c = 0; c != ncells; ++c) {
    const AmanziGeometry::Point& xc = mesh->getCellCentroid(c);
    source_c[0][c] = ana.source_exact(xc, 0);
  }
  op->global_operator()->UpdateRHS(source, false);

  // apply boundary conditions
  op->ApplyBCs(true, true, true);

  // set up the vectors
  CompositeVector& rhs = *op->global_operator()->rhs();
  CompositeVector u(rhs), error(rhs), flux(rhs);

  // fill the solution
  u.PutScalar(0.0);
  Epetra_MultiVector& u_c = *u.ViewComponent("cell", false);
  for (int c = 0; c != ncells; ++c) {
    const AmanziGeometry::Point& xc = mesh->getCellCentroid(c);
    u_c[0][c] = ana.pressure_exact(xc, 0.0);
  }

  if (u.HasComponent("face")) {
    Epetra_MultiVector& u_f = *u.ViewComponent("face", false);

    for (int f = 0; f != nfaces; ++f) {
      const AmanziGeometry::Point& xf = mesh->getFaceCentroid(f);
      u_f[0][f] = ana.pressure_exact(xf, 0.0);
    }
  }

  if (u.HasComponent("boundary_face")) {
    int nboundary_faces =
      mesh->getNumEntities(AmanziMesh::Entity_kind::BOUNDARY_FACE, AmanziMesh::Parallel_kind::OWNED);
    Epetra_MultiVector& u_f = *u.ViewComponent("boundary_face", false);

    for (int bf = 0; bf != nboundary_faces; ++bf) {
      int f = mesh->getMap(AmanziMesh::Entity_kind::FACE,false).LID(mesh->getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE,false).GID(bf));
      const AmanziGeometry::Point& xf = mesh->getFaceCentroid(f);
      u_f[0][f] = ana.pressure_exact(xf, 0.0);
    }
  }

  // forward apply
  double pnorm(0.0);
  u.ViewComponent("cell", false)->Norm2(&pnorm);

  op->global_operator()->ComputeResidual(u, error, true);
  double l2(0.0);
  error.ViewComponent("cell", false)->Norm2(&l2);
  l2 /= pnorm;
  double linf(0.0);
  error.NormInf(&linf);

  if (comm->MyPID() == 0) {
    printf("[%4d, %6.12e, %6.12e],\n", (int)round(log2(nx)), log2(l2), log2(linf));
  }
  return std::make_pair(log2(l2), log2(linf));
}


std::pair<double, double>
RunInverseProblem(const std::string& discretization, int nx, int ny, bool write_matrix)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();

  // create a mesh
  auto mesh_mstk = Teuchos::rcp(new Mesh_MSTK(0.0, 0.0, 1.0, 1.0, nx, ny, comm));
  auto mesh = Teuchos::rcp(new Mesh(mesh_mstk, Teuchos::rcp(new Amanzi::AmanziMesh::MeshFrameworkAlgorithms()), Teuchos::null)); 

  // modify diffusion coefficient
  Teuchos::RCP<std::vector<WhetStone::Tensor>> K =
    Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  int ncells = mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  int nfaces = mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
  int nfaces_wghost = mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::ALL);

  Analytic06 ana(mesh);

  for (int c = 0; c < ncells; c++) {
    const Point& xc = mesh->getCellCentroid(c);
    const WhetStone::Tensor& Kc = ana.TensorDiffusivity(xc, 0.0);
    K->push_back(Kc);
  }

  // create boundary data
  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(mesh, AmanziMesh::Entity_kind::FACE, WhetStone::DOF_Type::SCALAR));
  std::vector<int>& bc_model = bc->bc_model();
  std::vector<double>& bc_value = bc->bc_value();

  bool flag;
  for (int f = 0; f < nfaces_wghost; f++) {
    const Point& xf = mesh->getFaceCentroid(f);
    double area = mesh->getFaceArea(f);
    Point normal = ana.face_normal_exterior(f, &flag);

    if (fabs(xf[0]) < 1e-6) {
      bc_model[f] = Operators::OPERATOR_BC_NEUMANN;
      bc_value[f] = ana.velocity_exact(xf, 0.0) * normal / area;
    } else if (fabs(xf[1]) < 1e-6) {
      bc_model[f] = Operators::OPERATOR_BC_NEUMANN;
      bc_value[f] = ana.velocity_exact(xf, 0.0) * normal / area;
    } else if (fabs(xf[0] - 1.0) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6) {
      bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
      bc_value[f] = ana.pressure_exact(xf, 0.0);
    }
  }

  // create diffusion operator
  ParameterList op_list;
  op_list.set("discretization primary", discretization);
  op_list.set("nonlinear coefficient", "none");

  Operators::PDE_DiffusionFactory fac;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_c = mesh;
  Teuchos::RCP<PDE_Diffusion> op = fac.Create(op_list, mesh_c);
  op->SetBCs(bc, bc);
  op->SetTensorCoefficient(K);
  op->SetScalarCoefficient(Teuchos::null, Teuchos::null);
  op->UpdateMatrices(Teuchos::null, Teuchos::null);

  // create the rhs
  CompositeVector source(*op->global_operator()->rhs());
  Epetra_MultiVector& source_c = *source.ViewComponent("cell", false);
  for (int c = 0; c != ncells; ++c) {
    const AmanziGeometry::Point& xc = mesh->getCellCentroid(c);
    source_c[0][c] = ana.source_exact(xc, 0.0);
  }
  op->global_operator()->UpdateRHS(source, false);

  // apply boundary conditions
  op->ApplyBCs(true, true, true);

  // set up the vectors
  CompositeVector& rhs = *op->global_operator()->rhs();
  CompositeVector u(rhs), error(rhs), flux(rhs);

  // fill the solution
  u.PutScalar(0.0);
  Epetra_MultiVector& u_c = *u.ViewComponent("cell", false);
  for (int c = 0; c != ncells; ++c) {
    const AmanziGeometry::Point& xc = mesh->getCellCentroid(c);
    u_c[0][c] = ana.pressure_exact(xc, 0.0);
  }

  if (u.HasComponent("face")) {
    Epetra_MultiVector& u_f = *u.ViewComponent("face", false);

    for (int f = 0; f != nfaces; ++f) {
      const AmanziGeometry::Point& xf = mesh->getFaceCentroid(f);
      u_f[0][f] = ana.pressure_exact(xf, 0.0);
    }
  }

  if (u.HasComponent("boundary_face")) {
    int nboundary_faces =
      mesh->getNumEntities(AmanziMesh::Entity_kind::BOUNDARY_FACE, AmanziMesh::Parallel_kind::OWNED);
    Epetra_MultiVector& u_f = *u.ViewComponent("boundary_face", false);

    for (int bf = 0; bf != nboundary_faces; ++bf) {
      int f = mesh->getMap(AmanziMesh::Entity_kind::FACE,false).LID(mesh->getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE,false).GID(bf));

      const AmanziGeometry::Point& xf = mesh->getFaceCentroid(f);
      u_f[0][f] = ana.pressure_exact(xf, 0.0);
    }
  }

  Teuchos::ParameterList pc_list;
  pc_list.set("preconditioning method", "boomer amg");
  pc_list.sublist("boomer amg parameters").set("tolerance", 0.0);
  pc_list.set("iterative method", "gmres");
  op->global_operator()->set_inverse_parameters(pc_list);
  op->global_operator()->InitializeInverse();
  op->global_operator()->ComputeInverse();

  if (write_matrix) {
    std::stringstream fname;
    fname << "matrix_" << nx;
    if (discretization == "fv: default")
      fname << "_fv";
    else
      fname << "_mfd";
    fname << ".dat";
    EpetraExt::RowMatrixToMatlabFile(fname.str().c_str(), *op->global_operator()->A());
  }

  error.PutScalar(0.0);
  int ierr = op->global_operator()->ApplyInverse(*op->global_operator()->rhs(), error);
  CHECK(ierr >= 0);
  CHECK(op->global_operator()->num_itrs() < 10);

  error.Update(-1.0, u, 1.0);

  double pnorm(0.0);
  u.ViewComponent("cell", false)->Norm2(&pnorm);

  double l2(0.0);
  error.ViewComponent("cell", false)->Norm2(&l2);
  l2 /= pnorm;
  double linf(0.0);
  error.NormInf(&linf);

  if (comm->MyPID() == 0) {
    printf("[%4d, %6.12e, %6.12e],\n", (int)round(log2(nx)), log2(l2), log2(linf));
  }
  return std::make_pair(log2(l2), log2(linf));
}


void
MeanConvergenceRate(const std::vector<std::pair<double, double>>& l2s)
{
  double mean_dl2 = 0.0;
  double mean_dlinf = 0.0;
  int size = 0;
  for (int i = 1; i != l2s.size(); ++i) {
    mean_dl2 += (l2s[i].first - l2s[i - 1].first);
    mean_dlinf += (l2s[i].second - l2s[i - 1].second);
    size++;
  }

  double rate2 = -mean_dl2 / size;
  double rateinf = -mean_dlinf / size;
  std::cout << " Mean convergence rate (l2, linf) = " << rate2 << ", " << rateinf << std::endl;
  CHECK(rate2 > 1.9);
  CHECK(rateinf > 1.8);
}


void
RunForwardTest(const std::string& discretization)
{
  std::cout << "Convergence for discretization: " << discretization << std::endl;
  std::cout << "x = np.array([";
  std::vector<std::pair<double, double>> l2s;
  for (int i = 2; i <= 65; i *= 2) {
    std::pair<double, double> l2 = RunForwardProblem(discretization, i, i);
    l2s.push_back(l2);
  }
  std::cout << "])" << std::endl;

  MeanConvergenceRate(l2s);
}


void
RunInverseTest(const std::string& discretization)
{
  std::cout << "Convergence for discretization: " << discretization << std::endl;
  std::cout << "x = np.array([";
  std::vector<std::pair<double, double>> l2s;
  for (int i = 2; i <= 65; i *= 2) {
    std::pair<double, double> l2 = RunInverseProblem(discretization, i, i, false);
    l2s.push_back(l2);
  }
  std::cout << "])" << std::endl;

  MeanConvergenceRate(l2s);
}


TEST(OPERATOR_DIFFUSION_CONVERGENCE_FV)
{
  RunForwardTest("fv: default");
  RunInverseTest("fv: default");
}


TEST(OPERATOR_DIFFUSION_CONVERGENCE_MFD)
{
  RunForwardTest("mfd: default");
  RunInverseTest("mfd: default");
}

// just writes a typical matrix to file for debugging
// TEST(OPERATOR_DIFFUSION_WRITE_MFD) {
//   RunInverseProblem("mfd: default", 2, 2, true);
// }
