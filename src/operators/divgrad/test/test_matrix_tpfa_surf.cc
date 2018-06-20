#include "UnitTest++.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_RCP.hpp"
#include "EpetraExt_RowMatrixOut.h"
#include "Epetra_Util.h"

#include "MeshFactory.hh"
#include "Mesh.hh"

#include "Matrix_TPFA_Surf.hh"
#include "Matrix_TPFA.hh"
#include "MatrixMFD_TPFA.hh"

#include "LinearOperatorNKA.hh"

using namespace Amanzi;

struct mfd {
  Epetra_MpiComm *comm;
  Teuchos::RCP<AmanziMesh::Mesh> mesh;
  Teuchos::RCP<AmanziMesh::Mesh> surf_mesh;
  Teuchos::RCP<Teuchos::ParameterList> plist;
  Teuchos::RCP<Operators::Matrix_TPFA_Surf> A;
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
    std::cout << "On proc: " << surf_mesh->get_comm()->MyPID() << " there are " << surf_mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED) << " owned cells." << std::endl;
    
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
    if (method == "finite volume"){
      //plist.set("TPFA use cells only", true);
      plist.set("TPFA", false);
      plist.set("FV", true);
    }
    plist.set("dump Schur complement", true);
    // plist.sublist("preconditioner").set("preconditioner type", "boomer amg");
    // plist.sublist("preconditioner").sublist("boomer amg parameters")
    //   .set("cycle applications", 100);
    // plist.sublist("preconditioner").sublist("boomer amg parameters")
    //   .set("tolerance", 1.e-14);
    // plist.sublist("preconditioner").sublist("boomer amg parameters")
    //   .set("verbosity", 0);
    // plist.sublist("consistent face solver")
    //   .set("iterative method", "nka");
    // plist.sublist("consistent face solver").sublist("nka parameters")
    //   .set("error tolerance", 1.e-14);
    // plist.sublist("consistent face solver").sublist("nka parameters")
    //   .set("maximum number of iterations", 100);
    // plist.sublist("consistent face solver").sublist("preconditioner")
    //   .set("preconditioner type", "block ilu");
    plist.sublist("preconditioner").set("preconditioner type", "ml");
    plist.sublist("preconditioner").sublist("ml parameters")
      .set("cycle applications", 5);


    Teuchos::ParameterList plist_surf;
    plist_surf.set("MFD method", "two point flux approximation");
    plist_surf.set("dump Schur complement", true);
    plist_surf.sublist("preconditioner").set("preconditioner type", "boomer amg");
    plist_surf.sublist("preconditioner").sublist("boomer amg parameters")
      .set("cycle applications", 100);
    plist_surf.sublist("preconditioner").sublist("boomer amg parameters")
      .set("tolerance", 1.e-14);
    plist_surf.sublist("preconditioner").sublist("boomer amg parameters")
      .set("verbosity", 0);
    plist_surf.sublist("consistent face solver")
      .set("iterative method", "nka");
    plist_surf.sublist("consistent face solver").sublist("nka parameters")
      .set("error tolerance", 1.e-14);
    plist_surf.sublist("consistent face solver").sublist("nka parameters")
      .set("maximum number of iterations", 100);
    plist_surf.sublist("consistent face solver").sublist("preconditioner")
      .set("preconditioner type", "block ilu");

    // -- TPFA on surface
    As = Teuchos::rcp(new Operators::MatrixMFD_TPFA(plist_surf, surf_mesh));
    As->set_symmetric(false);
    As->SymbolicAssembleGlobalMatrices();
    As->CreateMFDmassMatrices(Teuchos::null);
    As->CreateMFDstiffnessMatrices(Teuchos::null);
    As->CreateMFDrhsVectors();
    As->ApplyBoundaryConditions(surf_bc_markers, surf_bc_values);
    
    // -- FV combined on subsurf
    A = Teuchos::rcp(new Operators::Matrix_TPFA_Surf(plist, mesh));
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
    
  void setSolution(const Teuchos::Ptr<CompositeVector>& x) {
    if (x->HasComponent("cell")) {
      Epetra_MultiVector& x_c = *x->ViewComponent("cell",false);
      for (int c=0; c!=x_c.MyLength(); ++c) {
	x_c[0][c] = value(mesh->cell_centroid(c));
      }
    }

    const Epetra_Map& fb_map = mesh->exterior_face_map(false);
    const Epetra_Map& f_map = mesh->face_map(false);


    if (x->HasComponent("boundary_face")) {
      Epetra_MultiVector& x_f = *x->ViewComponent("boundary_face",false);
      for (int f=0; f!=x_f.MyLength(); ++f) {
	int f_gid = fb_map.GID(f);
	int f_lid = f_map.LID(f_gid);
	x_f[0][f] = value(mesh->face_centroid(f_lid));
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
  createMFD("finite volume", Teuchos::null);

  x->PutScalar(1.);
  A->Apply(*x, *b);

  double norm;
  b->NormInf(&norm);
  //Epetra_MultiVector& b_cell = *b->ViewComponent("cell", false);
  //Epetra_MultiVector& b_bf = *b->ViewComponent("boundary_face", false);
  //std::cout << b_cell<<"\n";
  //std::cout << b_bf<<"\n";
  std::cout << "norm = " <<  norm << std::endl;
  CHECK_CLOSE(0., norm, 1.e-8);
}

TEST_FIXTURE(mfd, ApplyConstantTwoPointKr) {
  CompositeVectorSpace kr_sp;
  kr_sp.SetMesh(mesh)->SetGhosted()->SetComponent("face",AmanziMesh::FACE,1);

  Teuchos::RCP<CompositeVector> kr = 
    Teuchos::rcp(new CompositeVector(kr_sp));
  kr->PutScalar(0.5);
  createMFD("finite volume", kr.ptr());

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
  createMFD("finite volume", kr.ptr());
  setSolution(x.ptr());
  setSolutionSurf(xs.ptr());

  // test As * x - b == 0
  As->ComputeNegativeResidual(*xs, bs.ptr());
  double norm_surf;
  bs->NormInf(&norm_surf);
  std::cout << "norm surf = " << norm_surf << std::endl;
  CHECK_CLOSE(0., norm_surf, 1.e-8);

  // test Ax - b == 0
  A->ComputeNegativeResidual(*x, b.ptr());
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
  createMFD("finite volume", kr.ptr());
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


// TEST_FIXTURE(mfd, ConsistentFaceLinearTwoPointKr) {
//   CompositeVectorSpace kr_sp;
//   kr_sp.SetMesh(mesh)->SetGhosted()->SetComponent("face",AmanziMesh::FACE,1);

//   Teuchos::RCP<CompositeVector> kr = 
//     Teuchos::rcp(new CompositeVector(kr_sp));
//   kr->PutScalar(0.5);

//   setDirichletLinear();
//   setDirichletSurfLinear();
//   createMFD("finite volume", kr.ptr());
//   setSolution(x.ptr());

//   // test Aff^-1 * (rhs - Afc*x_c) - x_f == 0
//   *b = *x;
//   b->ViewComponent("face",false)->PutScalar(0.);
//   A->UpdateConsistentFaceConstraints(b.ptr());
//   b->Update(-1., *x, 1.);

//   double norm = 0.;
//   b->NormInf(&norm);
//   std::cout << "norm = " <<  norm << std::endl;
//   CHECK_CLOSE(0., norm, 1.e-8);
// }


TEST_FIXTURE(mfd, ApplyRandomTwoPointKr) {
  CompositeVectorSpace kr_sp;
  kr_sp.SetMesh(mesh)->SetGhosted()->SetComponent("face",AmanziMesh::FACE,1);

  Teuchos::RCP<CompositeVector> kr = 
    Teuchos::rcp(new CompositeVector(kr_sp));
  kr->PutScalar(0.5);

  if (mesh->get_comm()->MyPID() == 0)
    bc_markers[0] = Operators::MATRIX_BC_DIRICHLET;
  communicateBCs();
  createMFD("finite volume", kr.ptr());

  Epetra_MultiVector& b_c = *b->ViewComponent("cell",false);
  Epetra_Util rand;
  rand.SetSeed(1);
  for (int c=0; c!=b_c.MyLength(); ++c) {
    b_c[0][c] = rand.RandomDouble();
  }
  Epetra_MultiVector& b_f = *b->ViewComponent("boundary_face",false);
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

  if (mesh->get_comm()->MyPID() == 0)
    bc_markers[0] = Operators::MATRIX_BC_DIRICHLET;
  communicateBCs();
  createMFD("finite volume", kr.ptr());

  Epetra_MultiVector& x_c = *x->ViewComponent("cell",false);
  Epetra_Util rand;
  rand.SetSeed(1);
  for (int c=0; c!=x_c.MyLength(); ++c) {
    x_c[0][c] = rand.RandomDouble();
  }
  Epetra_MultiVector& x_f = *x->ViewComponent("boundary_face",false);
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


  if (mesh->get_comm()->MyPID() == 0)
    bc_markers[0] = Operators::MATRIX_BC_DIRICHLET;
  communicateBCs();
  //  createMFD("finite volume", kr.ptr());
  createMFD("finite volume", Teuchos::null);
  Epetra_Util rand2;
  rand2.SetSeed(1);
  Epetra_MultiVector& x_c = *x->ViewComponent("cell",false);

  for (int c=0; c!=x_c.MyLength(); ++c) {
    x_c[0][c] = rand2.RandomDouble();
  }
  Epetra_MultiVector& x_f = *x->ViewComponent("boundary_face",false);
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



