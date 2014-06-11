#include "UnitTest++.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_RCP.hpp"
#include "EpetraExt_RowMatrixOut.h"
#include "Epetra_Util.h"

#include "MeshFactory.hh"
#include "Mesh.hh"

#include "MatrixMFD_Surf.hh"
#include "MatrixMFD_TPFA_ScaledConstraint.hh"

#include "LinearOperatorNKA.hh"

using namespace Amanzi;

struct mfd {
  Epetra_MpiComm *comm;
  Teuchos::RCP<AmanziMesh::Mesh> mesh;
  Teuchos::RCP<AmanziMesh::Mesh> surf_mesh;
  Teuchos::RCP<Teuchos::ParameterList> plist;
  Teuchos::RCP<Operators::MatrixMFD_Surf> A;
  Teuchos::RCP<Operators::MatrixMFD_TPFA_ScaledConstraint> As;
  Teuchos::RCP<CompositeVector> x,b;
  Teuchos::RCP<CompositeVector> xs,bs;
  std::vector<Operators::MatrixBC> bc_markers;
  std::vector<double> bc_values;
  std::vector<Operators::MatrixBC> surf_bc_markers;
  std::vector<double> surf_bc_values;

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
    std::vector<std::string> surface_sets(1,"3D surface domain");
    surf_mesh = factory.create(&*mesh, surface_sets, AmanziMesh::FACE, true, false);

    // Boundary conditions
    int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
    bc_markers.resize(nfaces, Operators::MATRIX_BC_NULL);
    bc_values.resize(nfaces, 0.);

    int nfaces_surf = surf_mesh->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
    surf_bc_markers.resize(nfaces_surf, Operators::MATRIX_BC_NULL);
    surf_bc_values.resize(nfaces_surf, 0.);
  }

  void createMFD(std::string method, Teuchos::Ptr<CompositeVector> kr) {
    Teuchos::ParameterList plist;
    plist.set("MFD method", method);
    plist.sublist("preconditioner").set("preconditioner type", "boomer amg");
    plist.sublist("preconditioner").sublist("boomer amg parameters")
      .set("number of cycles", 100);
    plist.sublist("preconditioner").sublist("boomer amg parameters")
      .set("tolerance", 1.e-10);
    plist.sublist("preconditioner").sublist("boomer amg parameters")
      .set("verbosity", 0);
    plist.sublist("consistent face solver")
      .set("iterative method", "nka");
    plist.sublist("consistent face solver").sublist("nka parameters")
      .set("error tolerance", 1.e-10);
    plist.sublist("consistent face solver").sublist("nka parameters")
      .set("maximum number of iterations", 100);
    plist.sublist("consistent face solver").sublist("preconditioner")
      .set("preconditioner type", "block ilu");

    // -- TPFA on surface
    As = Teuchos::rcp(new Operators::MatrixMFD_TPFA_ScaledConstraint(plist, surf_mesh));
    As->set_symmetric(false);
    As->SymbolicAssembleGlobalMatrices();
    As->CreateMFDmassMatrices(Teuchos::null);
    As->CreateMFDstiffnessMatrices(Teuchos::null);
    As->CreateMFDrhsVectors();
    As->ApplyBoundaryConditions(surf_bc_markers, surf_bc_values);
    
    // -- combined on subsurf
    A = Teuchos::rcp(new Operators::MatrixMFD_Surf(plist, mesh));
    A->SetSurfaceOperator(As);
    A->set_symmetric(false);
    A->SymbolicAssembleGlobalMatrices();
    A->CreateMFDmassMatrices(Teuchos::null);
    A->CreateMFDstiffnessMatrices(kr);
    A->CreateMFDrhsVectors();
    A->ApplyBoundaryConditions(bc_markers, bc_values);
    
    x = Teuchos::rcp(new CompositeVector(A->DomainMap(), true));
    xs = Teuchos::rcp(new CompositeVector(As->DomainMap(), true));
    b = Teuchos::rcp(new CompositeVector(A->RangeMap(), true));
    bs = Teuchos::rcp(new CompositeVector(As->RangeMap(), true));
  }

  double value(AmanziGeometry::Point x) {
    return x[0] + x[1];
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

  void setDirichletSurfLinear() {
    for (int f=0; f!=surf_bc_markers.size(); ++f) {
      AmanziMesh::Entity_ID_List cells;
      surf_mesh->face_get_cells(f, AmanziMesh::USED, &cells);
      if (cells.size() == 1) {
	surf_bc_markers[f] = Operators::MATRIX_BC_DIRICHLET;
	surf_bc_values[f] = value(surf_mesh->face_centroid(f));
      }	
    }
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

  void setSolutionSurf(const Teuchos::Ptr<CompositeVector>& x) {
    if (x->HasComponent("cell")) {
      Epetra_MultiVector& x_c = *x->ViewComponent("cell",false);
      for (int c=0; c!=x_c.MyLength(); ++c) {
	x_c[0][c] = value(surf_mesh->cell_centroid(c));
      }
    }

    if (x->HasComponent("face")) {
      Epetra_MultiVector& x_f = *x->ViewComponent("face",false);
      for (int f=0; f!=x_f.MyLength(); ++f) {
	x_f[0][f] = value(surf_mesh->face_centroid(f));
      }
    }
  }
};

TEST_FIXTURE(mfd, ApplyConstantTwoPoint) {
  Teuchos::RCP<CompositeVector> kr = Teuchos::null;
  createMFD("two point flux approximation", kr.ptr());

  x->PutScalar(1.);
  A->Apply(*x, *b);

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
  createMFD("two point flux approximation", kr.ptr());

  x->PutScalar(1.);
  A->Apply(*x, *b);

  double norm;
  b->Norm2(&norm);
  std::cout << "norm = " <<  norm << std::endl;
  CHECK_CLOSE(0., norm, 1.e-8);
}

TEST_FIXTURE(mfd, ApplyLinearTwoPointKr) {
  CompositeVectorSpace kr_sp;
  kr_sp.SetMesh(mesh)->SetGhosted()->SetComponent("face",AmanziMesh::FACE,1);

  Teuchos::RCP<CompositeVector> kr = 
    Teuchos::rcp(new CompositeVector(kr_sp));
  kr->PutScalar(0.5);

  setDirichletLinear();
  setDirichletSurfLinear();
  createMFD("two point flux approximation", kr.ptr());
  setSolution(x.ptr());
  setSolutionSurf(xs.ptr());

  // test As * x - b == 0
  As->ComputeResidual(*xs, bs.ptr());
  double norm_surf;
  bs->Norm2(&norm_surf);
  std::cout << "norm surf = " << norm_surf << std::endl;
  CHECK_CLOSE(0., norm_surf, 1.e-8);

  // test Ax - b == 0
  A->ComputeResidual(*x, b.ptr());
  double norm;
  b->Norm2(&norm);
  std::cout << "norm = " <<  norm << std::endl;
  CHECK_CLOSE(0., norm, 1.e-8);
}

TEST_FIXTURE(mfd, ApplyInverseLinearTwoPointKr) {
  CompositeVectorSpace kr_sp;
  kr_sp.SetMesh(mesh)->SetGhosted()->SetComponent("face",AmanziMesh::FACE,1);

  Teuchos::RCP<CompositeVector> kr = 
    Teuchos::rcp(new CompositeVector(kr_sp));
  kr->PutScalar(0.5);

  setDirichletLinear();
  setDirichletSurfLinear();
  createMFD("two point flux approximation", kr.ptr());
  setSolution(x.ptr());

  // test A^1*b - x == 0
  b->PutScalar(0.);
  int ierr = A->ApplyInverse(*A->rhs(), *b);
  CHECK(!ierr);
  b->Update(-1., *x, 1.);

  double norm = 0.;
  b->Norm2(&norm);
  std::cout << "norm = " <<  norm << std::endl;
  CHECK_CLOSE(0., norm, 1.e-8);
}


TEST_FIXTURE(mfd, ConsistentFaceLinearTwoPointKr) {
  CompositeVectorSpace kr_sp;
  kr_sp.SetMesh(mesh)->SetGhosted()->SetComponent("face",AmanziMesh::FACE,1);

  Teuchos::RCP<CompositeVector> kr = 
    Teuchos::rcp(new CompositeVector(kr_sp));
  kr->PutScalar(0.5);

  setDirichletLinear();
  setDirichletSurfLinear();
  createMFD("two point flux approximation", kr.ptr());
  setSolution(x.ptr());

  // test Aff^-1 * (rhs - Afc*x_c) - x_f == 0
  *b = *x;
  b->ViewComponent("face",false)->PutScalar(0.);
  A->UpdateConsistentFaceConstraints(b.ptr());
  b->Update(-1., *x, 1.);

  double norm = 0.;
  b->Norm2(&norm);
  std::cout << "norm = " <<  norm << std::endl;
  CHECK_CLOSE(0., norm, 1.e-8);
}


TEST_FIXTURE(mfd, ApplyRandomTwoPointKr) {
  CompositeVectorSpace kr_sp;
  kr_sp.SetMesh(mesh)->SetGhosted()->SetComponent("face",AmanziMesh::FACE,1);

  Teuchos::RCP<CompositeVector> kr = 
    Teuchos::rcp(new CompositeVector(kr_sp));
  kr->PutScalar(0.5);

  bc_markers[0] = Operators::MATRIX_BC_DIRICHLET;
  createMFD("two point flux approximation", kr.ptr());

  Epetra_MultiVector& b_c = *b->ViewComponent("cell",false);
  Epetra_Util rand;
  rand.SetSeed(1);
  for (int c=0; c!=b_c.MyLength(); ++c) {
    b_c[0][c] = rand.RandomDouble();
  }
  Epetra_MultiVector& b_f = *b->ViewComponent("face",false);
  for (int f=0; f!=b_f.MyLength(); ++f) {
    b_f[0][f] = rand.RandomDouble();
  }
  CompositeVector r(*b);
  r = *b;

  // need a true solver, not just PC
  Teuchos::ParameterList solver_list;
  solver_list.set("error tolerance", 1.e-10);
  solver_list.sublist("VerboseObject").set("Verbosity Level", "extreme");
  
  AmanziSolvers::LinearOperatorNKA<CompositeMatrix,
		      CompositeVector,CompositeVectorSpace> solver(A,A);
  solver.Init(solver_list);

  // test A * A^1 * r - r == 0
  x->PutScalar(0.);

  int ierr = solver.ApplyInverse(*b, *x);
  CHECK(!ierr);
  CHECK(solver.num_itrs() <= 3);

  //int ierr = A->ApplyInverse(*b,*x);
  //CHECK(!ierr);

  b->PutScalar(0.);
  A->Apply(*x, *b);
  b->Update(-1., r, 1.);

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
  createMFD("two point flux approximation", kr.ptr());

  Epetra_MultiVector& x_c = *x->ViewComponent("cell",false);
  Epetra_Util rand;
  rand.SetSeed(1);
  for (int c=0; c!=x_c.MyLength(); ++c) {
    x_c[0][c] = rand.RandomDouble();
  }
  Epetra_MultiVector& x_f = *x->ViewComponent("face",false);
  for (int f=0; f!=x_f.MyLength(); ++f) {
    x_f[0][f] = rand.RandomDouble();
  }
  
  CompositeVector r(*x);
  r = *x;

  // need a true solver, not just PC
  Teuchos::ParameterList solver_list;
  solver_list.set("error tolerance", 1.e-10);
  
  AmanziSolvers::LinearOperatorGMRES<CompositeMatrix,
		      CompositeVector,CompositeVectorSpace> solver(A,A);
  solver.Init(solver_list);

  // test A * A^-1 * r - r == 0
  b->PutScalar(0.);
  A->Apply(*x, *b);
  x->PutScalar(0.);

  int ierr = solver.ApplyInverse(*b, *x);
  //  int ierr = A->ApplyInverse(*b, *x);
  CHECK(!ierr);
  CHECK(solver.num_itrs() <= 3);

  //int ierr = A->ApplyInverse(*b,*x);
  //CHECK(!ierr);

  x->Update(-1., r, 1.);

  double norm = 0.;
  x->Norm2(&norm);
  std::cout << "norm = " <<  norm << std::endl;
  CHECK_CLOSE(0., norm, 1.e-8);
}

TEST_FIXTURE(mfd, ApplyInverseRandomTwoPointKrRandom) {
  CompositeVectorSpace kr_sp;
  kr_sp.SetMesh(mesh)->SetGhosted()->SetComponent("face",AmanziMesh::FACE,1);

  Teuchos::RCP<CompositeVector> kr = 
    Teuchos::rcp(new CompositeVector(kr_sp));
  Epetra_Util rand;
  rand.SetSeed(2);
  Epetra_MultiVector& kr_f = *kr->ViewComponent("face",false);
  for (int f=0; f!=kr_f.MyLength(); ++f) {
    kr_f[0][f] = std::abs(rand.RandomDouble());
  }


  bc_markers[0] = Operators::MATRIX_BC_DIRICHLET;
  createMFD("two point flux approximation", kr.ptr());
  Epetra_Util rand2;
  rand2.SetSeed(1);
  Epetra_MultiVector& x_c = *x->ViewComponent("cell",false);

  for (int c=0; c!=x_c.MyLength(); ++c) {
    x_c[0][c] = rand2.RandomDouble();
  }
  Epetra_MultiVector& x_f = *x->ViewComponent("face",false);
  for (int f=0; f!=x_f.MyLength(); ++f) {
    x_f[0][f] = rand2.RandomDouble();
  }
  
  CompositeVector r(*x);
  r = *x;

  // need a true solver, not just PC
  Teuchos::ParameterList solver_list;
  solver_list.set("error tolerance", 1.e-10);
  
  AmanziSolvers::LinearOperatorGMRES<CompositeMatrix,
		      CompositeVector,CompositeVectorSpace> solver(A,A);
  solver.Init(solver_list);

  // test A * A^-1 * r - r == 0
  b->PutScalar(0.);
  A->Apply(*x, *b);
  x->PutScalar(0.);

  int ierr = solver.ApplyInverse(*b, *x);
  CHECK(!ierr);
  CHECK(solver.num_itrs() <= 3);

  //int ierr = A->ApplyInverse(*b,*x);
  //CHECK(!ierr);

  x->Update(-1., r, 1.);

  double norm = 0.;
  x->Norm2(&norm);
  std::cout << "norm = " <<  norm << std::endl;
  CHECK_CLOSE(0., norm, 1.e-8);
}



