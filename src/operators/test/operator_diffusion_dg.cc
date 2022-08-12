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
#include "GMVMesh.hh"
#include "CompositeVector.hh"
#include "MeshFactory.hh"
#include "NumericalIntegration.hh"
#include "Tensor.hh"

// Operators
#include "AnalyticDG01.hh"
#include "AnalyticDG02.hh"

#include "OperatorDefs.hh"
#include "PDE_DiffusionDG.hh"
#include "Verification.hh"


/* *****************************************************************
* This test diffusion solver with full tensor and source term.
* **************************************************************** */
class MyFunction : public Amanzi::WhetStone::WhetStoneFunction {
 public:
  MyFunction(AnalyticDG02* ana) : ana_(ana) {};
  ~MyFunction() {};

  virtual double Value(const Amanzi::AmanziGeometry::Point& x) const override { return (ana_->Tensor(x, 0.0))(0, 0); }
  AnalyticDG02* ana_;
};

void OperatorDiffusionDG(std::string solver_name,
                         std::string dg_basis = "regularized",
                         int dim = 2, int numi_order = 0) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0) std::cout << "\nTest: "<< dim << "D elliptic problem, dG method, solver: " 
                            << solver_name << ", basis=" << dg_basis 
                            << ", quadrature=" << numi_order << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_diffusion.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  // create a mesh framework
  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({Framework::MSTK, Framework::STK}));
  RCP<const Mesh> mesh;
  if (dim == 2) {
    // mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 1, 2);
    mesh = meshfactory.create("test/median7x8_filtered.exo");
    // mesh = meshfactory.create("test/triangular8_clockwise.exo");
  } else {
    mesh = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 2.0, 2, 3, 3, true, true);
  }

  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int ncells_wghost = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  // modify diffusion coefficient
  int d = mesh->space_dimension();
  auto Kc = std::make_shared<std::vector<WhetStone::Tensor> >();
  auto Kc_poly = std::make_shared<std::vector<WhetStone::MatrixPolynomial> >();
  auto Kc_func = std::make_shared<std::vector<WhetStone::WhetStoneFunction*> >();
  auto Kf = std::make_shared<std::vector<double> >();

  AnalyticDG02 ana(mesh, 2, false);
  MyFunction func(&ana);

  for (int c = 0; c < ncells_wghost; c++) {
    const Point& xc = mesh->cell_centroid(c);
    const WhetStone::Tensor& Ktmp = ana.Tensor(xc, 0.0);
    Kc->push_back(Ktmp);

    WhetStone::MatrixPolynomial Kpoly(d, d, d, 0);
    for (int i = 0; i < d; ++i) 
      for (int j = 0; j < d; ++j) Kpoly(i, j)(0) = Ktmp(i, j);
    Kpoly.set_origin(xc);
    Kc_poly->push_back(Kpoly);

    Kc_func->push_back(&func);
  }

  for (int f = 0; f < nfaces_wghost; f++) {
    double area = mesh->face_area(f);
    Kf->push_back(40.0 / area);
  }

  // create boundary data. We use full Taylor expansion of boundary data in
  // the vicinity of domain boundary.
  ParameterList op_list = plist.sublist("PK operator").sublist("diffusion operator dg");
  int order = op_list.sublist("schema").get<int>("method order");
  int nk = WhetStone::PolynomialSpaceDimension(dim, order);

  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(mesh, AmanziMesh::FACE, WhetStone::DOF_Type::VECTOR));
  std::vector<int>& bc_model = bc->bc_model();
  std::vector<std::vector<double> >& bc_value = bc->bc_value_vector(nk);

  WhetStone::Polynomial coefs;
  WhetStone::DenseVector data;

  const auto& fmap = mesh->face_map(true);
  const auto& bmap = mesh->exterior_face_map(true);
  for (int bf = 0; bf < bmap.NumMyElements(); ++bf) {
    int f = fmap.LID(bmap.GID(bf));
    const Point& xf = mesh->face_centroid(f);

    if (fabs(xf[0]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6 ||
        fabs(xf[1]) < 1e-6) {
      bc_model[f] = OPERATOR_BC_DIRICHLET;

      ana.SolutionTaylor(xf, 0.0, coefs);
      for (int i = 0; i < coefs.size(); ++i) {
        bc_value[f][i] = coefs(i);
      }
    } else {
      // bc_model[f] = OPERATOR_BC_NEUMANN;
      bc_model[f] = OPERATOR_BC_DIRICHLET;

      ana.SolutionTaylor(xf, 0.0, coefs);
      for (int i = 0; i < coefs.size(); ++i) {
        bc_value[f][i] = coefs(i);
      }
    }
  }

  // create diffusion operator 
  // -- primary term
  op_list.set<int>("quadrature order", numi_order);
  op_list.sublist("schema").set<std::string>("dg basis", dg_basis);
  Teuchos::RCP<PDE_DiffusionDG> op = Teuchos::rcp(new PDE_DiffusionDG(op_list, mesh));
  auto global_op = op->global_operator();
  const WhetStone::DG_Modal& dg = op->dg();

  // -- boundary conditions
  op->SetBCs(bc, bc);
  const CompositeVectorSpace& cvs = op->global_operator()->DomainMap();

  // create source and add it to the operator
  CompositeVector src(cvs);
  Epetra_MultiVector& src_c = *src.ViewComponent("cell");

  WhetStone::Polynomial pc(dim, order);
  WhetStone::NumericalIntegration numi(mesh);

  for (int c = 0; c < ncells; ++c) {
    const Point& xc = mesh->cell_centroid(c);

    ana.SourceTaylor(xc, 0.0, coefs);
    coefs.set_origin(xc);

    // -- calculate moments in natural basis
    data.Reshape(pc.size());
    for (auto it = pc.begin(); it < pc.end(); ++it) {
      int n = it.PolynomialPosition();

      WhetStone::Polynomial cmono(dim, it.multi_index(), 1.0);
      cmono.set_origin(xc);      
      WhetStone::Polynomial tmp = coefs * cmono;      

      data(n) = numi.IntegratePolynomialCell(c, tmp);
    } 

    // -- convert moment to my basis 
    dg.cell_basis(c).LinearFormNaturalToMy(data);
    for (int n = 0; n < pc.size(); ++n) {
      src_c[n][c] = data(n);
    }
  }

  // populate the diffusion operator
  if (numi_order == 0) 
    op->Setup(Kc, Kf);
  else 
    op->Setup(Kc_poly, Kf);
    // op->Setup(Kc_func, Kf);
  op->UpdateMatrices(Teuchos::null, Teuchos::null);

  // update the source term
  *global_op->rhs() = src;

  // apply BCs (primary=true, eliminate=true) and assemble
  op->ApplyBCs(true, true, true);

  // create preconditoner using the base operator class
  global_op->set_inverse_parameters("Hypre AMG", plist.sublist("preconditioners"));
  global_op->InitializeInverse();
  global_op->ComputeInverse();

  // Test SPD properties of the matrix.
  VerificationCV ver(global_op);
  ver.CheckMatrixSPD(false, true, 1);

  // create preconditoner with iterative method
  global_op->set_inverse_parameters("Hypre AMG", plist.sublist("preconditioners"), solver_name, plist.sublist("solvers"));
  global_op->InitializeInverse();
  global_op->ComputeInverse();
  
  CompositeVector& rhs = *global_op->rhs();
  CompositeVector solution(rhs);
  solution.PutScalar(0.0);

  global_op->ApplyInverse(rhs, solution);

  ver.CheckResidual(solution, 1.0e-11);

  if (MyPID == 0) {
    std::cout << "pressure solver (pcg): ||r||=" << global_op->residual() 
              << " itr=" << global_op->num_itrs()
              << " code=" << global_op->returned_code() 
              << " dofs=" << global_op->A()->NumGlobalRows() << std::endl;

    // visualization
    const Epetra_MultiVector& p = *solution.ViewComponent("cell");
    GMV::open_data_file(*mesh, (std::string)"operators.gmv");
    GMV::start_data();
    GMV::write_cell_data(p, 0, "solution");
    GMV::write_cell_data(p, 1, "gradx");
    GMV::write_cell_data(p, 2, "grady");
    GMV::close_data_file();
  }

  CHECK(global_op->num_itrs() < 200);

  // compute pressure error
  solution.ScatterMasterToGhosted();
  Epetra_MultiVector& p = *solution.ViewComponent("cell", false);

  double pnorm, pl2_err, pinf_err, pl2_mean, pinf_mean, pl2_int;
  ana.ComputeCellError(dg, p, 0.0, pnorm, pl2_err, pinf_err, pl2_mean, pinf_mean, pl2_int);

  if (MyPID == 0) {
    printf("Mean:     L2(p)=%9.6f  Inf(p)=%9.6f  itr=%3d\n", pl2_mean, pinf_mean, global_op->num_itrs());
    printf("Total:    L2(p)=%9.6f  Inf(p)=%9.6f\n", pl2_err, pinf_err);
    printf("Integral: L2(p)=%9.6f\n", pl2_int);

    CHECK(pl2_err < 1e-10);
  }
}

TEST(OPERATOR_DIFFUSION_DG) {
  OperatorDiffusionDG("AztecOO CG", "orthonormalized");
  OperatorDiffusionDG("AztecOO CG", "normalized");
  OperatorDiffusionDG("AztecOO CG");
  OperatorDiffusionDG("Amesos1");
  OperatorDiffusionDG("Amesos2: basker");
  OperatorDiffusionDG("Amesos2: superludist");

  int order = 2;
  OperatorDiffusionDG("AztecOO CG", "orthonormalized", 2, order);
}
