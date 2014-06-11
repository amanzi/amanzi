#include "UnitTest++.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_RCP.hpp"
#include "EpetraExt_RowMatrixOut.h"

#include "MeshFactory.hh"
#include "Mesh.hh"

#include "MatrixMFD.hh"
#include "MatrixMFD_Coupled.hh"

#include "LinearOperatorGMRES.hh"

using namespace Amanzi;

struct mfd {
  Epetra_MpiComm *comm;
  Teuchos::RCP<AmanziMesh::Mesh> mesh;
  Teuchos::RCP<Teuchos::ParameterList> plist;
  Teuchos::RCP<Operators::MatrixMFD> A;
  Teuchos::RCP<Operators::MatrixMFD> B;
  Teuchos::RCP<Operators::MatrixMFD_Coupled> C;
  Teuchos::RCP<CompositeVector> xA,xB;
  Teuchos::RCP<CompositeVector> bA,bB;
  Teuchos::RCP<TreeVector> x,b;
  Teuchos::RCP<TreeVector> offdiag;

  std::vector<Operators::MatrixBC> bc_markers;
  std::vector<double> bc_values;

  mfd() {
    comm = new Epetra_MpiComm(MPI_COMM_WORLD);

    plist = Teuchos::rcp(new Teuchos::ParameterList());
    Teuchos::updateParametersFromXmlFile("test-mesh.xml",plist.ptr());

    AmanziMesh::MeshFactory factory(comm);
    AmanziMesh::FrameworkPreference prefs(factory.preference());
    prefs.clear();
    prefs.push_back(AmanziMesh::MSTK);
    factory.preference(prefs);

    // create the meshes
    AmanziGeometry::GeometricModel gm(3, plist->sublist("Regions"), comm);
    mesh = factory.create(plist->sublist("Mesh").sublist("Generate Mesh"), &gm);

    // Boundary conditions
    int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
    bc_markers.resize(nfaces, Operators::MATRIX_BC_NULL);
    bc_values.resize(nfaces, 0.);
  }

  void createMFD(std::string method, std::string pc,
		 Teuchos::Ptr<CompositeVector> kr) {
    Teuchos::ParameterList plist;
    plist.set("MFD method", method);
    plist.set("dump Schur complement", true);
    plist.sublist("preconditioner").set("preconditioner type", pc);
    plist.sublist("preconditioner").sublist("boomer amg parameters")
      .set("cycle applications", 100);
    plist.sublist("preconditioner").sublist("boomer amg parameters")
      .set("tolerance", 1.e-14);
    plist.sublist("preconditioner").sublist("boomer amg parameters")
      .set("verbosity", 3);
    plist.sublist("preconditioner").sublist("boomer amg parameters")
      .set("number of functions", 2);
    plist.sublist("consistent face solver")
      .set("iterative method", "gmres");
    plist.sublist("consistent face solver").sublist("gmres parameters")
      .set("error tolerance", 1.e-14);
    plist.sublist("consistent face solver").sublist("gmres parameters")
      .set("maximum number of iterations", 100);
    plist.sublist("consistent face solver").sublist("preconditioner")
      .set("preconditioner type", "block ilu");

    A = Teuchos::rcp(new Operators::MatrixMFD(plist, mesh));
    B = Teuchos::rcp(new Operators::MatrixMFD(plist, mesh));
    C = Teuchos::rcp(new Operators::MatrixMFD_Coupled(plist, mesh));
    C->SetSubBlocks(A,B);
    C->SymbolicAssembleGlobalMatrices();

    // -- mass block
    A->set_symmetric(false);
    A->SymbolicAssembleGlobalMatrices();
    A->CreateMFDmassMatrices(Teuchos::null);
    A->CreateMFDstiffnessMatrices(kr);
    A->CreateMFDrhsVectors();
    A->ApplyBoundaryConditions(bc_markers, bc_values);

    // -- energy block
    B->set_symmetric(false);
    B->SymbolicAssembleGlobalMatrices();
    B->CreateMFDmassMatrices(Teuchos::null);
    B->CreateMFDstiffnessMatrices(Teuchos::null);
    B->CreateMFDrhsVectors();
    B->ApplyBoundaryConditions(bc_markers, bc_values);

    // -- vectors
    x = Teuchos::rcp(new TreeVector(C->DomainMap()));
    b = Teuchos::rcp(new TreeVector(C->RangeMap()));
    offdiag = Teuchos::rcp(new TreeVector(C->RangeMap()));
    offdiag->PutScalar(0.);

    xA = x->SubVector(0)->Data();
    xB = x->SubVector(1)->Data();
    bA = b->SubVector(0)->Data();
    bB = b->SubVector(1)->Data();

    // set up of block system
    C->SetOffDiagonals(offdiag->SubVector(0)->Data()->ViewComponent("cell",false),
		       offdiag->SubVector(1)->Data()->ViewComponent("cell",false));

  }


  double value(AmanziGeometry::Point x) {
    return x[0] + x[1] + x[2];
  }

  void setDirichletLinear() {
    for (int f=0; f!=bc_markers.size(); ++f) {
      AmanziMesh::Entity_ID_List cells;
      mesh->face_get_cells(f, AmanziMesh::USED, &cells);
      if (cells.size() == 1) {
	bc_markers[f] = Operators::MATRIX_BC_DIRICHLET;
	bc_values[f] = value(mesh->face_centroid(f));
      }	
    }
  }

  void setSolution(const Teuchos::Ptr<TreeVector>& x) {
    setSolution(x->SubVector(0)->Data().ptr());
    setSolution(x->SubVector(1)->Data().ptr());
  }

  void setSolution(const Teuchos::Ptr<CompositeVector>& x) {
    if (x->HasComponent("cell")) {
      Epetra_MultiVector& x_c = *x->ViewComponent("cell",false);
      for (int c=0; c!=x_c.MyLength(); ++c) {
	x_c[0][c] = value(mesh->cell_centroid(c));
      }
    }

    if (x->HasComponent("face")) {
      Epetra_MultiVector& x_f = *x->ViewComponent("face",false);
      for (int f=0; f!=x_f.MyLength(); ++f) {
	x_f[0][f] = value(mesh->face_centroid(f));
      }
    }
  }
};

TEST_FIXTURE(mfd, ApplyConstantTwoPoint) {
  Teuchos::RCP<CompositeVector> kr = Teuchos::null;
  createMFD("two point flux approximation", "boomer amg", kr.ptr());

  x->PutScalar(1.);
  C->Apply(*x, *b);

  double norm;
  b->Norm2(&norm);
  std::cout << "norm = " <<  norm << std::endl;
  CHECK_CLOSE(0., norm, 1.e-8);
}

TEST_FIXTURE(mfd, ApplyConstantTwoPointKr) {
  CompositeVectorSpace kr_sp;
  kr_sp.SetMesh(mesh)->SetGhosted()->SetComponent("face",AmanziMesh::FACE,1);

  Teuchos::RCP<CompositeVector> kr = 
    Teuchos::rcp(new CompositeVector(kr_sp));
  kr->PutScalar(0.5);
  createMFD("two point flux approximation", "boomer amg", kr.ptr());

  x->PutScalar(1.);
  C->Apply(*x, *b);

  double norm;
  b->Norm2(&norm);
  std::cout << "norm = " <<  norm << std::endl;
  CHECK_CLOSE(0., norm, 1.e-8);
}

TEST_FIXTURE(mfd, SymmetryRandomTwoPointKr) {
  bc_markers[0] = Operators::MATRIX_BC_DIRICHLET;
  createMFD("two point flux approximation", "block ilu", Teuchos::null);

  Epetra_MultiVector& bA_c = *bA->ViewComponent("cell",false);
  std::srand(0);
  for (int c=0; c!=bA_c.MyLength(); ++c) {
    bA_c[0][c] = (std::rand() % 1000) / 1000.0;
  }
  Epetra_MultiVector& bA_f = *bA->ViewComponent("face",false);
  for (int f=0; f!=bA_f.MyLength(); ++f) {
    bA_f[0][f] = (std::rand() % 1000) / 1000.0;
  }
  Epetra_MultiVector& bB_c = *bB->ViewComponent("cell",false);
  std::srand(0);
  for (int c=0; c!=bB_c.MyLength(); ++c) {
    bB_c[0][c] = (std::rand() % 1000) / 1000.0;
  }
  Epetra_MultiVector& bB_f = *bB->ViewComponent("face",false);
  for (int f=0; f!=bB_f.MyLength(); ++f) {
    bB_f[0][f] = (std::rand() % 1000) / 1000.0;
  }

  TreeVector r(*b);

  // check data starts out T/p symmetric
  r = *b;
  r.SubVector(0)->Update(1., *r.SubVector(1), -1.);
  double norm = 0.;
  r.SubVector(0)->Norm2(&norm);
  std::cout << "norm on IC = " << norm << std::endl;
  CHECK_CLOSE(0., norm, 1.e-8);

  // check forward op maintains T/p symmetry
  x->PutScalar(0.);
  int ierr = C->Apply(*b,*x);
  CHECK(!ierr);

  r = *x;
  r.SubVector(0)->Update(1., *r.SubVector(1), -1.);
  norm = 0.;
  r.SubVector(0)->Norm2(&norm);
  std::cout << "norm on forward op = " << norm << std::endl;
  CHECK_CLOSE(0., norm, 1.e-8);

  // check inverse operator maintains T/p symmetry
  x->PutScalar(0.);
  C->ApplyInverse(*b, *x);

  r = *x;
  r.SubVector(0)->Update(1., *r.SubVector(1), -1.);
  norm = 0.;
  r.SubVector(0)->Norm2(&norm);
  std::cout << "norm on inverse = " << norm << std::endl;
  CHECK_CLOSE(0., norm, 1.e-8);
}


TEST_FIXTURE(mfd, ApplyRandomTwoPointKr) {
  CompositeVectorSpace kr_sp;
  kr_sp.SetMesh(mesh)->SetGhosted()->SetComponent("face",AmanziMesh::FACE,1);

  Teuchos::RCP<CompositeVector> kr = 
    Teuchos::rcp(new CompositeVector(kr_sp));
  kr->PutScalar(0.5);

  bc_markers[0] = Operators::MATRIX_BC_DIRICHLET;
  createMFD("two point flux approximation", "boomer amg", Teuchos::null);

  Epetra_MultiVector& bA_c = *bA->ViewComponent("cell",false);
  std::srand(0);
  for (int c=0; c!=bA_c.MyLength(); ++c) {
    bA_c[0][c] = (std::rand() % 1000) / 1000.0;
  }
  Epetra_MultiVector& bA_f = *bA->ViewComponent("face",false);
  for (int f=0; f!=bA_f.MyLength(); ++f) {
    bA_f[0][f] = (std::rand() % 1000) / 1000.0;
  }
  Epetra_MultiVector& bB_c = *bB->ViewComponent("cell",false);
  std::srand(0);
  for (int c=0; c!=bB_c.MyLength(); ++c) {
    bB_c[0][c] = (std::rand() % 1000) / 1000.0;
  }
  Epetra_MultiVector& bB_f = *bB->ViewComponent("face",false);
  for (int f=0; f!=bB_f.MyLength(); ++f) {
    bB_f[0][f] = (std::rand() % 1000) / 1000.0;
  }

  TreeVector r(*b);
  r = *b;

  // need a true solver, not just PC
  Teuchos::ParameterList solver_list;
  solver_list.set("error tolerance", 1.e-10);
  
  AmanziSolvers::LinearOperatorGMRES<TreeMatrix,
		      TreeVector,TreeVectorSpace> solver(C,C);
  solver.Init(solver_list);

  // test A * A^1 * r - r == 0
  x->PutScalar(0.);

  // int ierr = solver.ApplyInverse(*b, *x);
  // CHECK(!ierr);
  // CHECK(solver.num_itrs() <= 3);

  int ierr = C->ApplyInverse(*b,*x);
  CHECK(!ierr);

  b->PutScalar(0.);
  C->Apply(*x, *b);
  b->Update(-1., r, 1.);
  std::cout << "A(Inv(b))-b = " << std::endl;
  b->Print(std::cout);

  double norm = 0.;
  b->Norm2(&norm);
  std::cout << "norm = " <<  norm << std::endl;
  CHECK_CLOSE(0., norm, 1.e-8);
}


TEST_FIXTURE(mfd, ApplyInverseRandomTwoPointKr) {
  CompositeVectorSpace kr_sp;
  kr_sp.SetMesh(mesh)->SetGhosted()->SetComponent("face",AmanziMesh::FACE,1);

  Teuchos::RCP<CompositeVector> kr = 
    Teuchos::rcp(new CompositeVector(kr_sp));
  kr->PutScalar(0.5);

  bc_markers[0] = Operators::MATRIX_BC_DIRICHLET;
  createMFD("two point flux approximation", "boomer amg", Teuchos::null);

  Epetra_MultiVector& xA_c = *xA->ViewComponent("cell",false);
  std::srand(0);
  for (int c=0; c!=xA_c.MyLength(); ++c) {
    xA_c[0][c] = (std::rand() % 1000) / 1000.0;
  }
  Epetra_MultiVector& xA_f = *xA->ViewComponent("face",false);
  for (int f=0; f!=xA_f.MyLength(); ++f) {
    xA_f[0][f] = (std::rand() % 1000) / 1000.0;
  }
  Epetra_MultiVector& xB_c = *xB->ViewComponent("cell",false);
  std::srand(0);
  for (int c=0; c!=xB_c.MyLength(); ++c) {
    xB_c[0][c] = (std::rand() % 1000) / 1000.0;
  }
  Epetra_MultiVector& xB_f = *xB->ViewComponent("face",false);
  for (int f=0; f!=xB_f.MyLength(); ++f) {
    xB_f[0][f] = (std::rand() % 1000) / 1000.0;
  }

  TreeVector r(*x);
  r = *x;

  // need a true solver, not just PC
  Teuchos::ParameterList solver_list;
  solver_list.set("error tolerance", 1.e-10);
  
  AmanziSolvers::LinearOperatorGMRES<TreeMatrix,
		      TreeVector,TreeVectorSpace> solver(C,C);
  solver.Init(solver_list);

  // test A * A^-1 * r - r == 0
  b->PutScalar(0.);
  C->Apply(*x, *b);
  x->PutScalar(0.);

  // int ierr = solver.ApplyInverse(*b, *x);
  // CHECK(!ierr);
  // CHECK(solver.num_itrs() <= 3);

  int ierr = C->ApplyInverse(*b,*x);
  CHECK(!ierr);

  x->Update(-1., r, 1.);
  std::cout << "Inv(A(x))-x = " << std::endl;
  x->Print(std::cout);

  double norm = 0.;
  x->Norm2(&norm);
  std::cout << "norm = " <<  norm << std::endl;
  CHECK_CLOSE(0., norm, 1.e-8);
}

