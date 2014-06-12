#include "UnitTest++.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_RCP.hpp"
#include "EpetraExt_RowMatrixOut.h"

#include "MeshFactory.hh"
#include "Mesh.hh"

#include "MatrixMFD_TPFA_ScaledConstraint.hh"

#include "LinearOperatorGMRES.hh"

using namespace Amanzi;

struct mfd {
  Epetra_MpiComm *comm;
  Teuchos::RCP<AmanziMesh::Mesh> mesh;
  Teuchos::RCP<Teuchos::ParameterList> plist;
  Teuchos::RCP<Operators::MatrixMFD_TPFA_ScaledConstraint> A;
  Teuchos::RCP<CompositeVector> x,b;
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

  void createMFD(std::string method, Teuchos::Ptr<CompositeVector> kr) {
    Teuchos::ParameterList plist;
    plist.set("MFD method", method);
    plist.sublist("preconditioner").set("preconditioner type", "boomer amg");
    plist.sublist("preconditioner").sublist("boomer amg parameters")
      .set("cycle applications", 100);
    plist.sublist("preconditioner").sublist("boomer amg parameters")
      .set("tolerance", 1.e-14);
    plist.sublist("preconditioner").sublist("boomer amg parameters")
      .set("verbosity", 0);
    plist.sublist("consistent face solver")
      .set("iterative method", "gmres");
    plist.sublist("consistent face solver").sublist("gmres parameters")
      .set("error tolerance", 1.e-10);
    plist.sublist("consistent face solver").sublist("gmres parameters")
      .set("maximum number of iterations", 100);
    plist.sublist("consistent face solver").sublist("preconditioner")
      .set("preconditioner type", "block ilu");

    A = Teuchos::rcp(new Operators::MatrixMFD_TPFA_ScaledConstraint(plist, mesh));

    // -- TPF on surface
    A->set_symmetric(false);
    A->SymbolicAssembleGlobalMatrices();
    A->CreateMFDmassMatrices(Teuchos::null);
    A->CreateMFDstiffnessMatrices(kr);
    A->CreateMFDrhsVectors();
    A->ApplyBoundaryConditions(bc_markers, bc_values);
    
    x = Teuchos::rcp(new CompositeVector(A->DomainMap(), true));
    b = Teuchos::rcp(new CompositeVector(A->RangeMap(), true));
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

TEST_FIXTURE(mfd, ApplyConstantTPFA_ScaledConstraint) {
  Teuchos::RCP<CompositeVector> kr = Teuchos::null;
  createMFD("two point flux approximation", kr.ptr());

  x->PutScalar(1.);
  A->Apply(*x, *b);

  double norm;
  b->NormInf(&norm);
  std::cout << "norm = " <<  norm << std::endl;
  CHECK_CLOSE(0., norm, 1.e-8);
}

TEST_FIXTURE(mfd, ApplyConstantTPFA_ScaledConstraintKr) {
  CompositeVectorSpace kr_sp;
  kr_sp.SetMesh(mesh)->SetGhosted()->SetComponent("face",AmanziMesh::FACE,1);

  Teuchos::RCP<CompositeVector> kr = 
    Teuchos::rcp(new CompositeVector(kr_sp));
  kr->PutScalar(0.5);
  createMFD("two point flux approximation", kr.ptr());

  x->PutScalar(1.);
  A->Apply(*x, *b);

  double norm;
  b->NormInf(&norm);
  std::cout << "norm = " <<  norm << std::endl;
  CHECK_CLOSE(0., norm, 1.e-8);
}

TEST_FIXTURE(mfd, ApplyLinearTPFA_ScaledConstraintKr) {
  CompositeVectorSpace kr_sp;
  kr_sp.SetMesh(mesh)->SetGhosted()->SetComponent("face",AmanziMesh::FACE,1);

  Teuchos::RCP<CompositeVector> kr = 
    Teuchos::rcp(new CompositeVector(kr_sp));
  kr->PutScalar(0.5);

  setDirichletLinear();
  createMFD("two point flux approximation", kr.ptr());
  setSolution(x.ptr());

  // test Ax - b == 0
  A->ComputeResidual(*x, b.ptr());
  double norm;
  b->NormInf(&norm);
  std::cout << "norm = " <<  norm << std::endl;
  CHECK_CLOSE(0., norm, 1.e-8);
}

TEST_FIXTURE(mfd, ApplyInverseLinearTPFA_ScaledConstraintKr) {
  CompositeVectorSpace kr_sp;
  kr_sp.SetMesh(mesh)->SetGhosted()->SetComponent("face",AmanziMesh::FACE,1);

  Teuchos::RCP<CompositeVector> kr = 
    Teuchos::rcp(new CompositeVector(kr_sp));
  kr->PutScalar(0.5);

  setDirichletLinear();
  createMFD("two point flux approximation", kr.ptr());
  setSolution(x.ptr());

  // test A^1*b - x == 0
  b->PutScalar(0.);
  int ierr = A->ApplyInverse(*A->rhs(), *b);
  CHECK(!ierr);
  b->Update(-1., *x, 1.);

  double norm = 0.;
  b->NormInf(&norm);
  std::cout << "norm = " <<  norm << std::endl;
  CHECK_CLOSE(0., norm, 1.e-8);
}


TEST_FIXTURE(mfd, ConsistentFaceLinearTPFA_ScaledConstraintKr) {
  CompositeVectorSpace kr_sp;
  kr_sp.SetMesh(mesh)->SetGhosted()->SetComponent("face",AmanziMesh::FACE,1);

  Teuchos::RCP<CompositeVector> kr = 
    Teuchos::rcp(new CompositeVector(kr_sp));
  kr->PutScalar(0.5);

  setDirichletLinear();
  createMFD("two point flux approximation", kr.ptr());
  setSolution(x.ptr());

  // test Aff^-1 * (rhs - Afc*x_c) - x_f == 0
  *b = *x;
  b->ViewComponent("face",false)->PutScalar(0.);
  A->UpdateConsistentFaceConstraints(b.ptr());
  b->Update(-1., *x, 1.);

  double norm = 0.;
  b->NormInf(&norm);
  std::cout << "norm = " <<  norm << std::endl;
  CHECK_CLOSE(0., norm, 1.e-8);
}


TEST_FIXTURE(mfd, ApplyRandomTPFA_ScaledConstraintKr) {
  CompositeVectorSpace kr_sp;
  kr_sp.SetMesh(mesh)->SetGhosted()->SetComponent("face",AmanziMesh::FACE,1);

  Teuchos::RCP<CompositeVector> kr = 
    Teuchos::rcp(new CompositeVector(kr_sp));
  kr->PutScalar(0.5);

  bc_markers[0] = Operators::MATRIX_BC_DIRICHLET;
  createMFD("two point flux approximation", kr.ptr());

  Epetra_MultiVector& b_c = *b->ViewComponent("cell",false);
  std::srand(0);
  for (int c=0; c!=b_c.MyLength(); ++c) {
    b_c[0][c] = (std::rand() % 1000) / 1000.0;
  }
  Epetra_MultiVector& b_f = *b->ViewComponent("face",false);
  for (int f=0; f!=b_f.MyLength(); ++f) {
    b_f[0][f] = (std::rand() % 1000) / 1000.0;
  }
  CompositeVector r(*b);
  r = *b;

  // test A * A^1 * r - r == 0
  x->PutScalar(0.);

  int ierr = A->ApplyInverse(*b,*x);
  CHECK(!ierr);

  b->PutScalar(0.);
  A->Apply(*x, *b);
  b->Update(-1., r, 1.);

  double norm = 0.;
  b->NormInf(&norm);
  std::cout << "norm = " <<  norm << std::endl;
  CHECK_CLOSE(0., norm, 1.e-8);
}


TEST_FIXTURE(mfd, ApplyInverseRandomTPFA_ScaledConstraintKr) {
  CompositeVectorSpace kr_sp;
  kr_sp.SetMesh(mesh)->SetGhosted()->SetComponent("face",AmanziMesh::FACE,1);

  Teuchos::RCP<CompositeVector> kr = 
    Teuchos::rcp(new CompositeVector(kr_sp));
  kr->PutScalar(0.5);

  bc_markers[0] = Operators::MATRIX_BC_DIRICHLET;
  createMFD("two point flux approximation", kr.ptr());

  Epetra_MultiVector& x_c = *x->ViewComponent("cell",false);
  std::srand(0);
  for (int c=0; c!=x_c.MyLength(); ++c) {
    x_c[0][c] = (std::rand() % 1000) / 1000.0;
  }
  Epetra_MultiVector& x_f = *x->ViewComponent("face",false);
  for (int f=0; f!=x_f.MyLength(); ++f) {
    x_f[0][f] = (std::rand() % 1000) / 1000.0;
  }
  
  CompositeVector r(*x);
  r = *x;

  // test A * A^-1 * r - r == 0
  b->PutScalar(0.);
  A->Apply(*x, *b);
  x->PutScalar(0.);

  int ierr = A->ApplyInverse(*b,*x);
  CHECK(!ierr);

  x->Update(-1., r, 1.);

  double norm = 0.;
  x->NormInf(&norm);
  std::cout << "norm = " <<  norm << std::endl;
  CHECK_CLOSE(0., norm, 1.e-8);
}



