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
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

// Amanzi
#include "MeshFactory.hh"
#include "OutputXDMF.hh"
#include "Tensor.hh"
#include "WhetStoneDefs.hh"

// Amanzi::Operators
#include "Analytic00.hh"
#include "Analytic07.hh"
#include "Operator.hh"
#include "OperatorDefs.hh"
#include "PDE_DiffusionFactory.hh"
#include "ReconstructionCellLinear.hh"

using namespace Teuchos;
using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;
using namespace Amanzi::Operators;

// given a point in mesh1, it returns face in mesh2
int
find_face(const Point& xf1, const Point& ray, const Mesh& mesh2)
{
  std::string rgn("MIDDLE2");
  AmanziMesh::cEntity_ID_View block;
  AmanziMesh::cDouble_View vofs;

  Kokkos::tie(block, vofs) =
    mesh2.getSetEntitiesAndVolumeFractions(rgn, Entity_kind::FACE, Parallel_kind::OWNED);
  int nblock = block.size();

  // check of ray hits interior of a face
  for (int n = 0; n < nblock; ++n) {
    int f2 = block[n];
    const auto& xf2 = mesh2.getFaceCentroid(f2);
    const auto& normal = mesh2.getFaceNormal(f2);

    double angle = ray * normal;
    if (std::fabs(angle) < 1e-12) AMANZI_ASSERT(false);
    double s = ((xf1 - xf2) * normal) / angle;
    auto xf1_proj = xf1 - s * ray;

    auto coords2 = mesh2.getFaceCoordinates(f2);
    if ((coords2[0] - xf1_proj) * (coords2[1] - xf1_proj) <= 0.0) return f2;
  }

  // check for closest face
  int f2min(-1);
  double dist_min(1e+50);

  for (int n = 0; n < nblock; ++n) {
    int f2 = block[n];
    const auto& xf2 = mesh2.getFaceCentroid(f2);
    const auto& normal = mesh2.getFaceNormal(f2);

    double angle = ray * normal;
    double s = ((xf1 - xf2) * normal) / angle;
    auto xf1_proj = xf1 - s * ray;

    auto coords2 = mesh2.getFaceCoordinates(f2);
    s = std::min(norm(coords2[0] - xf1_proj), norm(coords2[1] - xf1_proj));
    if (s < dist_min) {
      dist_min = s;
      f2min = f2;
    }
  }

  return f2min;
}


// conservative remap of fluxes from mesh1 to mesh2 using particles
const int REFINE = 10;

std::map<int, double>
interpolate_fluxes(const Mesh& mesh1, const Mesh& mesh2, const Epetra_MultiVector& flux1)
{
  std::map<int, double> flux2;

  std::string rgn1("MIDDLE1"), rgn2("MIDDLE2");
  AmanziMesh::cEntity_ID_View block1, block2;
  AmanziMesh::cDouble_View vofs1, vofs2;

  Kokkos::tie(block1, vofs1) =
    mesh1.getSetEntitiesAndVolumeFractions(rgn1, Entity_kind::FACE, Parallel_kind::OWNED);

  Kokkos::tie(block2, vofs2) =
    mesh2.getSetEntitiesAndVolumeFractions(rgn2, Entity_kind::FACE, Parallel_kind::OWNED);

  int nblock1 = block1.size();
  int nblock2 = block2.size();
  double s;

  for (int n = 0; n < nblock1; ++n) {
    int f1 = block1[n];
    const auto& normal1 = mesh1.getFaceNormal(f1);
    auto coords1 = mesh1.getFaceCoordinates(f1);

    for (int k = 0; k < REFINE; ++k) {
      s = double(k + 0.5) / REFINE;
      auto xa1 = coords1[0] * (1.0 - s) + coords1[1] * s;

      int f2 = find_face(xa1, normal1, mesh2);
      s = flux2[f2];
      flux2[f2] += (flux1[0][f1] / REFINE) / mesh2.getFaceArea(f2);
    }
  }

  return flux2;
}


/* *****************************************************************
* TBW.
* **************************************************************** */
void
RunTest(int n1, int n2, double gap)
{
  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0) std::cout << "\nTest: Diffusion Schwarz solver on two meshes" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_diffusion_schwarz.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  std::vector<double> middle2({ 1.0 + gap, 0.0 });
  plist->sublist("regions")
    .sublist("MIDDLE2")
    .sublist("region: plane")
    .set<Teuchos::Array<double>>("point", middle2);

  // create an SIMPLE mesh framework
  ParameterList region_list = plist->sublist("regions");
  auto gm = Teuchos::rcp(new GeometricModel(2, region_list, *comm));

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));

  double mid1(1.0), mid2(1.0 + gap);
  RCP<Mesh> mesh1 = meshfactory.create(0.0, 0.0, mid1, 1.0, n1, n1);
  RCP<Mesh> mesh2 = meshfactory.create(mid2, 0.0, 2.0, 1.0, n2, n2);

  int nnodes1_owned = mesh1->getNumEntities(Entity_kind::NODE, Parallel_kind::OWNED);
  int ncells1_owned = mesh1->getNumEntities(Entity_kind::CELL, Parallel_kind::OWNED);
  int nfaces1_wghost = mesh1->getNumEntities(Entity_kind::FACE, Parallel_kind::ALL);

  int nnodes2_owned = mesh2->getNumEntities(Entity_kind::NODE, Parallel_kind::OWNED);
  int ncells2_owned = mesh2->getNumEntities(Entity_kind::CELL, Parallel_kind::OWNED);
  int nfaces2_wghost = mesh2->getNumEntities(Entity_kind::FACE, Parallel_kind::ALL);

  // extract right side of domain 1
  // NOTE: we have to initialize sets using a uniform mesh
  AmanziMesh::cEntity_ID_View block;
  AmanziMesh::cDouble_View vofs;
  std::set<int> side1, side2;

  Kokkos::tie(block, vofs) =
    mesh1->getSetEntitiesAndVolumeFractions("MIDDLE1", Entity_kind::FACE, Parallel_kind::OWNED);
  for (int n = 0; n < block.size(); ++n) side1.insert(block[n]);

  Kokkos::tie(block, vofs) =
    mesh2->getSetEntitiesAndVolumeFractions("MIDDLE2", Entity_kind::FACE, Parallel_kind::OWNED);
  for (int n = 0; n < block.size(); ++n) side2.insert(block[n]);

  // -- modify meshes
  for (int n = 0; n < nnodes1_owned; ++n) {
    auto xp = mesh1->getNodeCoordinate(n);
    double s = std::sin(1.5 * xp[0]) * std::sin(3 * M_PI * xp[1]);
    xp[0] += 0.03 * s;
    xp[1] += 0.02 * s;
    mesh1->setNodeCoordinate(n, xp);
  }
  mesh1->recacheGeometry();

  for (int n = 0; n < nnodes2_owned; ++n) {
    auto xp = mesh2->getNodeCoordinate(n);
    if (fabs(xp[0] - 2.0) > 1e-6) {
      double s = std::sin(1.5 * xp[0]) * std::sin(3 * M_PI * xp[1]);
      xp[0] += 0.02 * s;
      xp[1] += 0.02 * s;
      mesh2->setNodeCoordinate(n, xp);
    }
  }
  mesh2->recacheGeometry();

  // create boundary data
  // Analytic00 ana1(mesh1, 2), ana2(mesh2, 2);
  Analytic07 ana1(mesh1), ana2(mesh2);

  // create solutions
  auto cvs1 = Teuchos::rcp(new CompositeVectorSpace());
  cvs1->SetMesh(mesh1)
    ->SetGhosted(true)
    ->AddComponent("cell", Entity_kind::CELL, 1)
    ->AddComponent("face", Entity_kind::FACE, 1);

  auto sol1 = Teuchos::rcp(new CompositeVector(*cvs1));
  auto flux1 = Teuchos::rcp(new CompositeVector(*cvs1));
  sol1->PutScalar(0.0);

  auto cvs2 = Teuchos::rcp(new CompositeVectorSpace());
  cvs2->SetMesh(mesh2)
    ->SetGhosted(true)
    ->AddComponent("cell", Entity_kind::CELL, 1)
    ->AddComponent("face", Entity_kind::FACE, 1);

  auto sol2 = Teuchos::rcp(new CompositeVector(*cvs2));
  sol2->PutScalar(0.0);

  auto& sol1_f = *sol1->ViewComponent("face");
  auto& sol2_f = *sol2->ViewComponent("face");
  auto& flux1_f = *flux1->ViewComponent("face");

  auto sol2_c = sol2->ViewComponent("cell", true);

  // create gradients
  Teuchos::ParameterList grad_list = plist->sublist("reconstruction");
  auto gradient2 = Teuchos::rcp(new Operators::ReconstructionCellLinear(mesh2));
  gradient2->Init(grad_list);

  // ------------------
  // Schwarz iterations
  // ------------------
  double lambda(1.0);
  std::map<int, double> interface;

  for (int loop = 0; loop < 100; ++loop) {
    // Dirichlet solver
    gradient2->Compute(sol2_c);

    auto bc1 = Teuchos::rcp(new BCs(mesh1, Entity_kind::FACE, WhetStone::DOF_Type::SCALAR));
    std::vector<int>& bc1_model = bc1->bc_model();
    std::vector<double>& bc1_value = bc1->bc_value();

    double err(0.0);
    for (int f = 0; f < nfaces1_wghost; ++f) {
      const Point& xf = mesh1->getFaceCentroid(f);
      const Point& normal = mesh1->getFaceNormal(f);

      if (fabs(xf[0]) < 1e-6) {
        bc1_model[f] = OPERATOR_BC_DIRICHLET;
        bc1_value[f] = ana1.pressure_exact(xf, 0.0);
        // } else if (fabs(xf[0] - mid1) < 1e-6) {
      } else if (side1.find(f) != side1.end()) {
        int f2 = find_face(xf, normal, *mesh2);
        int c2 = getFaceOnBoundaryInternalCell(*mesh2, f2);
        double val = gradient2->getValue(c2, xf);

        bc1_model[f] = OPERATOR_BC_DIRICHLET;
        bc1_value[f] = 0.1 * interface[f] + 0.9 * val;
        err += std::pow(interface[f] - val, 2) * mesh1->getFaceArea(f);
        interface[f] = val;
      } else if (fabs(xf[1]) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6) {
        bc1_model[f] = OPERATOR_BC_DIRICHLET;
        bc1_value[f] = ana1.pressure_exact(xf, 0.0);
      }
    }
    if (loop > 0 && err < 1e-12) break;

    Teuchos::ParameterList olist = plist->sublist("PK operators").sublist("diffusion operator");
    Operators::PDE_DiffusionFactory opfactory1(olist, mesh1);
    // opfactory.SetVariableTensorCoefficient(K);

    Teuchos::RCP<Operators::PDE_Diffusion> op1 = opfactory1.Create();
    op1->SetBCs(bc1, bc1);

    Teuchos::RCP<Operator> global_op1 = op1->global_operator();
    global_op1->Init();

    auto& rhs1 = *global_op1->rhs()->ViewComponent("cell");
    for (int c = 0; c < ncells1_owned; c++) {
      const Point& xc = mesh1->getCellCentroid(c);
      rhs1[0][c] = ana1.source_exact(xc, 0.0) * mesh1->getCellVolume(c);
    }

    op1->UpdateMatrices(Teuchos::null, Teuchos::null);
    op1->ApplyBCs(true, true, true);

    global_op1->set_inverse_parameters(
      "Hypre AMG", plist->sublist("preconditioners"), "PCG", plist->sublist("solvers"));
    global_op1->InitializeInverse();

    global_op1->ApplyInverse(*global_op1->rhs(), *sol1);
    op1->UpdateFlux(sol1.ptr(), flux1.ptr());

    // Neumann solver
    auto flux2 = interpolate_fluxes(*mesh1, *mesh2, flux1_f);

    auto bc2 = Teuchos::rcp(new BCs(mesh2, Entity_kind::FACE, WhetStone::DOF_Type::SCALAR));
    std::vector<int>& bc2_model = bc2->bc_model();
    std::vector<double>& bc2_value = bc2->bc_value();
    std::vector<double>& bc2_mixed = bc2->bc_mixed();

    for (int f = 0; f < nfaces2_wghost; ++f) {
      const Point& xf = mesh2->getFaceCentroid(f);
      // if (fabs(xf[0] - mid2) < 1e-6) {
      if (side2.find(f) != side2.end()) {
        double val = flux2[f];
        bc2_model[f] = OPERATOR_BC_MIXED;
        bc2_value[f] = -val - lambda * sol2_f[0][f];
        bc2_mixed[f] = lambda;
      } else if (fabs(xf[0] - 2.0) < 1e-6) {
        bc2_model[f] = OPERATOR_BC_DIRICHLET;
        bc2_value[f] = ana2.pressure_exact(xf, 0.0);
      } else if (fabs(xf[1]) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6) {
        bc2_model[f] = OPERATOR_BC_DIRICHLET;
        bc2_value[f] = ana2.pressure_exact(xf, 0.0);
      }
    }

    Operators::PDE_DiffusionFactory opfactory2(olist, mesh2);
    Teuchos::RCP<Operators::PDE_Diffusion> op2 = opfactory2.Create();
    op2->SetBCs(bc2, bc2);

    Teuchos::RCP<Operator> global_op2 = op2->global_operator();
    global_op2->Init();

    auto& rhs2 = *global_op2->rhs()->ViewComponent("cell");
    for (int c = 0; c < ncells2_owned; c++) {
      const Point& xc = mesh2->getCellCentroid(c);
      rhs2[0][c] = ana2.source_exact(xc, 0.0) * mesh2->getCellVolume(c);
    }

    op2->UpdateMatrices(Teuchos::null, Teuchos::null);
    op2->ApplyBCs(true, true, true);

    global_op2->set_inverse_parameters(
      "Hypre AMG", plist->sublist("preconditioners"), "PCG", plist->sublist("solvers"));
    global_op2->InitializeInverse();

    global_op2->ApplyInverse(*global_op2->rhs(), *sol2);

    // statistics
    if (MyPID == 0) {
      std::cout << "loop: " << loop << " err=" << std::sqrt(err)
                << " itrs: " << global_op1->num_itrs() << " " << global_op2->num_itrs();
      if (global_op1->returned_code() != 1 || global_op2->returned_code() != 1) {
        std::cout << "||r||=" << global_op1->residual() << " " << global_op2->residual()
                  << " code=" << global_op1->returned_code() << " " << global_op2->returned_code();
      }
      std::cout << std::endl;
    }
  }

  // calculate error in potential
  double pnorm, l2_err, inf_err;
  Epetra_MultiVector& p = *sol1->ViewComponent("cell");

  ana1.ComputeCellError(p, 0.0, pnorm, l2_err, inf_err);
  CHECK(l2_err < 3e-2);

  if (MyPID == 0) {
    l2_err /= pnorm;
    printf("L2(p)=%9.6f  Inf(p)=%9.6f  #dofs=%d\n", l2_err, inf_err, ncells1_owned + ncells2_owned);
  }

  // initialize I/O
  const auto& field1_c = *sol1->ViewComponent("cell");
  const auto& field2_c = *sol2->ViewComponent("cell");

  Teuchos::ParameterList iolist;
  iolist.set<std::string>("file name base", "plot1");
  OutputXDMF io1(iolist, mesh1, true, false);

  io1.InitializeCycle(0.0, 1, "");
  io1.WriteVector(*field1_c(0), "solution1", AmanziMesh::Entity_kind::CELL);
  io1.FinalizeCycle();

  iolist.set<std::string>("file name base", "plot2");
  OutputXDMF io2(iolist, mesh2, true, false);

  io2.InitializeCycle(0.0, 1, "");
  io2.WriteVector(*field2_c(0), "solution2", AmanziMesh::Entity_kind::CELL);
  io2.FinalizeCycle();
}


TEST(DIRICHLET_NEUMANN)
{
  RunTest(19, 13, 0.0025); // gap = 0.0025;
  // RunTest(39, 25, 0.00125);
}
