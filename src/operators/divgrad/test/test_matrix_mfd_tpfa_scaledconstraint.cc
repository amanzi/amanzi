#include "UnitTest++.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_RCP.hpp"
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_CrsMatrixIn.h"
#include "Epetra_SerialComm.h"

#include "MeshFactory.hh"
#include "Mesh.hh"

#include "MatrixMFD_TPFA_ScaledConstraint.hh"

#include "LinearOperatorGMRES.hh"

using namespace Amanzi;

struct mfd {
  Epetra_MpiComm *comm;
  Teuchos::RCP<AmanziGeometry::GeometricModel> gm;
  Teuchos::RCP<AmanziMesh::Mesh> mesh;
  Teuchos::RCP<Teuchos::ParameterList> plist;
  Teuchos::RCP<Operators::MatrixMFD_TPFA_ScaledConstraint> A;
  Teuchos::RCP<CompositeVector> x,b;
  std::vector<Operators::MatrixBC> bc_markers;
  std::vector<double> bc_values;

  mfd() {
    comm = new Epetra_MpiComm(MPI_COMM_WORLD);

    plist = Teuchos::rcp(new Teuchos::ParameterList());
    Teuchos::updateParametersFromXmlFile("test/test-mesh.xml",plist.ptr());

    AmanziMesh::MeshFactory factory(comm);
    AmanziMesh::FrameworkPreference prefs(factory.preference());
    prefs.clear();
    prefs.push_back(AmanziMesh::MSTK);
    factory.preference(prefs);

    // create the meshes
    Teuchos::ParameterList& regionlist = plist->sublist("Regions");
    gm = Teuchos::rcp(new AmanziGeometry::GeometricModel(3, regionlist, comm));
    mesh = factory.create(plist->sublist("Mesh").sublist("Generate Mesh"), &*gm);

    // Boundary conditions
    int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);
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
      mesh->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
      if (cells.size() == 1) {
	bc_markers[f] = Operators::MATRIX_BC_DIRICHLET;
	bc_values[f] = value(mesh->face_centroid(f));
      }	
    }
  }

  void setDirichletOne() {
    AmanziMesh::Entity_ID_List bottom;
    mesh->get_set_entities("bottom side", AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL, &bottom);
    for (int f=0; f!=bottom.size(); ++f) {
      bc_markers[bottom[f]] = Operators::MATRIX_BC_DIRICHLET;
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



TEST_FIXTURE(mfd, AssembleRandomTPFASCNormed) {
  CompositeVectorSpace kr_sp;
  kr_sp.SetMesh(mesh)->SetGhosted()->SetComponent("face",AmanziMesh::FACE,1);

  Teuchos::RCP<CompositeVector> kr = 
    Teuchos::rcp(new CompositeVector(kr_sp));
  Epetra_MultiVector& kr_f = *kr->ViewComponent("face",false);
  for (int f=0; f!=kr_f.MyLength(); ++f) {
    AmanziGeometry::Point fc = mesh->face_centroid(f);
    kr_f[0][f] = std::sqrt(std::abs(fc[0]) + std::abs(fc[1]) + std::abs(fc[2]));
  }

  setDirichletOne();
  createMFD("two point flux approximation", kr.ptr());
  
  Epetra_MultiVector& x_c = *x->ViewComponent("cell",false);
  for (int c=0; c!=x_c.MyLength(); ++c) {
    AmanziGeometry::Point fc = mesh->cell_centroid(c);
    x_c[0][c] = 3. * std::pow(fc[0],2) - 1.123 * std::sqrt(std::abs(fc[1])) + std::pow(fc[2],3);
  }

  b->PutScalar(0.);
  A->Schur()->Apply(x_c, *b->ViewComponent("cell",false));
  double norm_Acc(0.);
  b->ViewComponent("cell",false)->Norm2(&norm_Acc);
  std::cout << std::setprecision(15) << "norm = " << norm_Acc << std::endl;

  CHECK_CLOSE(3.54171767185051, norm_Acc, 1.e-8);

  // const Epetra_Map& fmap = mesh->face_map(false);
  // const Epetra_Map& fmap_ghosted = mesh->face_map(true);
  // Epetra_Vector Aff_diag(fmap);
  // A->Aff()->ExtractDiagonalCopy(Aff_diag);

  // for (int f=0; f!=fmap_ghosted.NumMyElements(); ++f) {
  //   AmanziGeometry::Point fc = mesh->face_centroid(f);
  //   if (std::abs(fc[0] - 0.5) < 1.e-8 && std::abs(fc[1] - 0.166666666666666666) < 1.e-8 &&
  //       std::abs(fc[2] - 0.3333333333333) < 1.e-8) {
  //     std::cout << "We is here!" << std::endl;
  //     AmanziMesh::Entity_ID_List cells;
  //     mesh->face_get_cells(f, AmanziMesh::Parallel_type::OWNED, &cells);
  //     std::cout << "On proc: " << comm->MyPID() << "face gid: " << fmap_ghosted.GID(f) << " is in " << cells.size() << " cells." << std::endl;
  //   }
  // }

  // for (int pid=0; pid!=comm->NumProc(); ++pid) {
  //   if (pid == comm->MyPID()) {
  //     for (int f=0; f!=fmap.NumMyElements(); ++f) {
  //       AmanziGeometry::Point fc = mesh->face_centroid(f);
  //       std::cout << "GID: " << fmap.GID(f) << ", Centroid: " << fc << " val=" << Aff_diag[f] << " b=" << b_f[0][f] << std::endl;
  //     }
  //   }
  //   comm->Barrier();
  // }

  // dump matrices for later comparison
  std::stringstream filename_Sff;
  filename_Sff << "test/MatrixMFD_TPFA_SC_Scc_np" << comm->NumProc() << ".txt";
  EpetraExt::RowMatrixToMatlabFile(filename_Sff.str().c_str(), *A->Schur());

  if (comm->MyPID() == 0) {
    Epetra_SerialComm mycomm;
    
    // load matrices for comparison
    std::stringstream filename_Sff_ref;
    filename_Sff_ref << "test/MatrixMFD_TPFA_SC_Scc_ref_np" << comm->NumProc() << ".txt";
    Epetra_CrsMatrix* Sref;
    EpetraExt::MatlabFileToCrsMatrix(filename_Sff_ref.str().c_str(), mycomm, Sref);

    // load matrices for comparison
    std::stringstream filename_Sff_test;
    filename_Sff_test << "test/MatrixMFD_TPFA_SC_Scc_np" << comm->NumProc() << ".txt";
    Epetra_CrsMatrix* Stest;
    EpetraExt::MatlabFileToCrsMatrix(filename_Sff_test.str().c_str(), mycomm, Stest);
    
    // compare
    int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
    for (int c=0; c!=ncells; ++c) {
      std::vector<int> inds(20);
      std::vector<double> vals(20);
      std::vector<int> inds_ref(20);
      std::vector<double> vals_ref(20);
      int num_entries;

      Stest->ExtractMyRowCopy(c, 20, num_entries, &vals[0], &inds[0]);
      inds.resize(num_entries);
      vals.resize(num_entries);
      
      Sref->ExtractMyRowCopy(c, 20, num_entries, &vals_ref[0], &inds_ref[0]);
      inds_ref.resize(num_entries);
      vals_ref.resize(num_entries);

      CHECK(inds_ref == inds);
      CHECK(vals_ref == vals);
    }

    delete Sref;
    delete Stest;
  }
  
}

