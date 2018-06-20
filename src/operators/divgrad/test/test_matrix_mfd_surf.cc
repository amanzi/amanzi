#include "UnitTest++.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_RCP.hpp"
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_CrsMatrixIn.h"
#include "Epetra_SerialComm.h"
#include "Epetra_Util.h"

#include "MeshFactory.hh"
#include "Mesh.hh"

#include "MatrixMFD_Surf.hh"
#include "MatrixMFD_TPFA.hh"

#include "LinearOperatorNKA.hh"

using namespace Amanzi;

struct mfd {
  Epetra_MpiComm *comm;
  Teuchos::RCP<AmanziGeometry::GeometricModel> gm;
  Teuchos::RCP<AmanziMesh::Mesh> mesh;
  Teuchos::RCP<AmanziMesh::Mesh> surf_mesh;
  Teuchos::RCP<Teuchos::ParameterList> plist;
  Teuchos::RCP<Operators::MatrixMFD_Surf> A;
  Teuchos::RCP<Operators::MatrixMFD_TPFA> As;
  Teuchos::RCP<CompositeVector> x,b;
  Teuchos::RCP<CompositeVector> xs,bs;
  std::vector<Operators::MatrixBC> bc_markers;
  std::vector<double> bc_values;
  std::vector<Operators::MatrixBC> surf_bc_markers;
  std::vector<double> surf_bc_values;

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

    std::vector<std::string> surface_sets(1,"3D surface domain");
    surf_mesh = factory.create(&*mesh, surface_sets, AmanziMesh::FACE, true, false);
    
    // Boundary conditions
    int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);
    bc_markers.resize(nfaces, Operators::MATRIX_BC_NULL);
    bc_values.resize(nfaces, 0.);

    int nfaces_surf = surf_mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);
    surf_bc_markers.resize(nfaces_surf, Operators::MATRIX_BC_NULL);
    surf_bc_values.resize(nfaces_surf, 0.);
  }

  void createMFD(std::string method, Teuchos::Ptr<CompositeVector> kr) {
    Teuchos::ParameterList plist;
    plist.set("MFD method", method);
    plist.set("dump Schur complement", true);
    plist.sublist("preconditioner").set("preconditioner type", "boomer amg");
    plist.sublist("preconditioner").sublist("boomer amg parameters")
      .set("cycle applications", 100);
    plist.sublist("preconditioner").sublist("boomer amg parameters")
      .set("tolerance", 1.e-14);
    plist.sublist("preconditioner").sublist("boomer amg parameters")
      .set("verbosity", 0);
    plist.sublist("consistent face solver")
      .set("iterative method", "nka");
    plist.sublist("consistent face solver").sublist("nka parameters")
      .set("error tolerance", 1.e-14);
    plist.sublist("consistent face solver").sublist("nka parameters")
      .set("maximum number of iterations", 100);
    plist.sublist("consistent face solver").sublist("preconditioner")
      .set("preconditioner type", "block ilu");

    // -- TPFA on surface
    As = Teuchos::rcp(new Operators::MatrixMFD_TPFA(plist, surf_mesh));
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

  void communicateBCs() {
    CompositeVectorSpace space;
    space.SetMesh(mesh)->SetGhosted()->SetComponent("face", AmanziMesh::FACE, 2);
    CompositeVector bcs(space);

    {
      Epetra_MultiVector& bcs_f = *bcs.ViewComponent("face",false);
      for (int f=0; f!=bcs_f.MyLength(); ++f) {
        if (bc_markers[f] == Operators::MATRIX_BC_DIRICHLET) {
          bcs_f[0][f] = 1.0;
          bcs_f[1][f] = bc_values[f];
        }	
      }
    }
    bcs.ScatterMasterToGhosted("face");

    const Epetra_MultiVector& bcs_f = *bcs.ViewComponent("face",true);
    for (int f=0; f!=bcs_f.MyLength(); ++f) {
      if (bcs_f[0][f] > 0) {
        bc_markers[f] = Operators::MATRIX_BC_DIRICHLET;
        bc_values[f] = bcs_f[1][f];
      }
    }
  }

  void communicateBCsSurf() {
    CompositeVectorSpace space;
    space.SetMesh(surf_mesh)->SetGhosted()->SetComponent("face", AmanziMesh::FACE, 2);
    CompositeVector bcs(space);

    {
      Epetra_MultiVector& bcs_f = *bcs.ViewComponent("face",false);
      for (int f=0; f!=bcs_f.MyLength(); ++f) {
        if (surf_bc_markers[f] == Operators::MATRIX_BC_DIRICHLET) {
          bcs_f[0][f] = 1.0;
          bcs_f[1][f] = surf_bc_values[f];
        }	
      }
    }
    bcs.ScatterMasterToGhosted("face");

    const Epetra_MultiVector& bcs_f = *bcs.ViewComponent("face",true);
    for (int f=0; f!=bcs_f.MyLength(); ++f) {
      if (bcs_f[0][f] > 0) {
        surf_bc_markers[f] = Operators::MATRIX_BC_DIRICHLET;
        surf_bc_values[f] = bcs_f[1][f];
      }
    }
  }

  void setDirichletLinear() {
    int nfaces = mesh->num_entities(AmanziMesh::FACE,AmanziMesh::Parallel_type::OWNED);
    for (int f=0; f!=nfaces; ++f) {
      AmanziMesh::Entity_ID_List cells;
      mesh->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
      if (cells.size() == 1) {
        bc_markers[f] = Operators::MATRIX_BC_DIRICHLET;
        bc_values[f] = value(mesh->face_centroid(f));
      }	
    }
    communicateBCs();
  }

  void setDirichletSurfLinear() {
    int nfaces = surf_mesh->num_entities(AmanziMesh::FACE,AmanziMesh::Parallel_type::OWNED);
    for (int f=0; f!=nfaces; ++f) {
      AmanziMesh::Entity_ID_List cells;
      surf_mesh->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
      if (cells.size() == 1) {
        surf_bc_markers[f] = Operators::MATRIX_BC_DIRICHLET;
        surf_bc_values[f] = value(surf_mesh->face_centroid(f));
      }	
    }
    communicateBCsSurf();
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
  createMFD("two point flux approximation", Teuchos::null);

  x->PutScalar(1.);
  A->Apply(*x, *b);

  double norm;
  b->NormInf(&norm);
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
  b->NormInf(&norm);
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
  bs->NormInf(&norm_surf);
  std::cout << "norm surf = " << norm_surf << std::endl;
  CHECK_CLOSE(0., norm_surf, 1.e-8);

  // test Ax - b == 0
  A->ComputeResidual(*x, b.ptr());
  double norm;
  b->NormInf(&norm);
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
  b->NormInf(&norm);
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
  b->NormInf(&norm);
  std::cout << "norm = " <<  norm << std::endl;
  CHECK_CLOSE(0., norm, 1.e-8);
}


TEST_FIXTURE(mfd, ApplyRandomTwoPointKr) {
  CompositeVectorSpace kr_sp;
  kr_sp.SetMesh(mesh)->SetGhosted()->SetComponent("face",AmanziMesh::FACE,1);

  Teuchos::RCP<CompositeVector> kr = 
    Teuchos::rcp(new CompositeVector(kr_sp));
  kr->PutScalar(0.5);

  setDirichletOne();
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


TEST_FIXTURE(mfd, ApplyInverseRandomTwoPointKr) {
  CompositeVectorSpace kr_sp;
  kr_sp.SetMesh(mesh)->SetGhosted()->SetComponent("face",AmanziMesh::FACE,1);

  Teuchos::RCP<CompositeVector> kr = 
    Teuchos::rcp(new CompositeVector(kr_sp));
  kr->PutScalar(0.5);

  setDirichletOne();
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

  // test A * A^-1 * r - r == 0
  b->PutScalar(0.);
  A->Apply(*x, *b);
  x->PutScalar(0.);

  int ierr = A->ApplyInverse(*b, *x);
  CHECK(!ierr);

  x->Update(-1., r, 1.);

  double norm = 0.;
  x->NormInf(&norm);
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
    kr_f[0][f] = std::max(std::abs(rand.RandomDouble()), 0.01);
  }


  setDirichletOne();
  createMFD("two point flux approximation", kr.ptr());
  //createMFD("two point flux approximation", Teuchos::null);
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

TEST_FIXTURE(mfd, AssembleRandomSurfNormed) {
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
  
  Epetra_MultiVector& x_f = *x->ViewComponent("face",false);
  for (int f=0; f!=kr_f.MyLength(); ++f) {
    AmanziGeometry::Point fc = mesh->face_centroid(f);
    x_f[0][f] = 3. * std::pow(fc[0],2) - 1.123 * std::sqrt(std::abs(fc[1])) + std::pow(fc[2],3);
  }

  b->PutScalar(0.);
  A->Aff()->Apply(x_f, *b->ViewComponent("face",false));
  double norm_Aff(0.);
  b->ViewComponent("face",false)->Norm2(&norm_Aff);

  b->PutScalar(0.);
  A->Schur()->Apply(x_f, *b->ViewComponent("face",false));
  double norm_schur(0.);

  Epetra_MultiVector& b_f = *b->ViewComponent("face",false);
  b_f.Norm2(&norm_schur);
  std::cout << std::setprecision(15) << "norms = " << norm_Aff << ", " << norm_schur << std::endl;

  CHECK_CLOSE(15.7418101076706, norm_Aff, 1.e-8);
  CHECK_CLOSE(5.13793706832201, norm_schur, 1.e-8);

  const Epetra_Map& fmap = mesh->face_map(false);
  const Epetra_Map& fmap_ghosted = mesh->face_map(true);
  Epetra_Vector Aff_diag(fmap);
  A->Aff()->ExtractDiagonalCopy(Aff_diag);

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
  std::stringstream filename_Aff;
  filename_Aff << "test/MatrixMFD_Surf_Aff_np" << comm->NumProc() << ".txt";
  EpetraExt::RowMatrixToMatlabFile(filename_Aff.str().c_str(), *A->Aff());

  std::stringstream filename_Sff;
  filename_Sff << "test/MatrixMFD_Surf_Sff_np" << comm->NumProc() << ".txt";
  EpetraExt::RowMatrixToMatlabFile(filename_Sff.str().c_str(), *A->Schur());

  if (comm->MyPID() == 0) {
    Epetra_SerialComm mycomm;
    
    // load matrices for comparison
    std::stringstream filename_Aff_ref;
    filename_Aff_ref << "test/MatrixMFD_Surf_Aff_ref_np" << comm->NumProc() << ".txt";
    Epetra_CrsMatrix* Aref;
    EpetraExt::MatlabFileToCrsMatrix(filename_Aff_ref.str().c_str(), mycomm, Aref);

    std::stringstream filename_Sff_ref;
    filename_Sff_ref << "test/MatrixMFD_Surf_Sff_ref_np" << comm->NumProc() << ".txt";
    Epetra_CrsMatrix* Sref;
    EpetraExt::MatlabFileToCrsMatrix(filename_Sff_ref.str().c_str(), mycomm, Sref);

    // load matrices for comparison
    std::stringstream filename_Aff_test;
    filename_Aff_test << "test/MatrixMFD_Surf_Aff_np" << comm->NumProc() << ".txt";
    Epetra_CrsMatrix* Atest;
    EpetraExt::MatlabFileToCrsMatrix(filename_Aff_test.str().c_str(), mycomm, Atest);

    std::stringstream filename_Sff_test;
    filename_Sff_test << "test/MatrixMFD_Surf_Sff_np" << comm->NumProc() << ".txt";
    Epetra_CrsMatrix* Stest;
    EpetraExt::MatlabFileToCrsMatrix(filename_Sff_test.str().c_str(), mycomm, Stest);
    
    // compare
    for (int f=0; f!=fmap.NumGlobalElements(); ++f) {
      std::vector<int> inds(20);
      std::vector<double> vals(20);
      std::vector<int> inds_ref(20);
      std::vector<double> vals_ref(20);
      int num_entries;

      Atest->ExtractMyRowCopy(f, 20, num_entries, &vals[0], &inds[0]);
      inds.resize(num_entries);
      vals.resize(num_entries);
    
      Aref->ExtractMyRowCopy(f, 20, num_entries, &vals_ref[0], &inds_ref[0]);
      inds_ref.resize(num_entries);
      vals_ref.resize(num_entries);

      CHECK(inds_ref == inds);
      CHECK(vals_ref == vals);

      Stest->ExtractMyRowCopy(f, 20, num_entries, &vals[0], &inds[0]);
      inds.resize(num_entries);
      vals.resize(num_entries);
      
      Sref->ExtractMyRowCopy(f, 20, num_entries, &vals_ref[0], &inds_ref[0]);
      inds_ref.resize(num_entries);
      vals_ref.resize(num_entries);

      CHECK(inds_ref == inds);
      CHECK(vals_ref == vals);
    }

    delete Aref;
    delete Sref;
    delete Atest;
    delete Stest;
  }
  
}



