#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_RCP.hpp"
#include "EpetraExt_RowMatrixOut.h"

#include "MeshFactory.hh"
#include "Mesh.hh"

#include "matrix_mfd_tpfa.hh"
#include "matrix_mfd_surf.hh"

using namespace Amanzi;

struct test_mfd {
  Epetra_MpiComm *comm;
  Teuchos::RCP<AmanziMesh::Mesh> mesh;
  Teuchos::RCP<AmanziMesh::Mesh> surf_mesh;
  Teuchos::RCP<Teuchos::ParameterList> plist;

  Teuchos::RCP<Operators::MatrixMFD_TPFA> Atpf;
  Teuchos::RCP<Operators::MatrixMFD_Surf> As;

  test_mfd() {
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

    // -- Boundary conditions
    int nfaces_surf = surf_mesh->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
    std::vector<Operators::Matrix_bc> bc_surf_markers(nfaces_surf,
            Operators::MATRIX_BC_NULL);
    std::vector<double> bc_surf_values(nfaces_surf, 0.);

    int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
    std::vector<Operators::Matrix_bc> bc_markers(nfaces, Operators::MATRIX_BC_NULL);
    std::vector<double> bc_values(nfaces, 0.);

    // create the matrices
    Teuchos::ParameterList mfd_plist;
    mfd_plist.set("MFD method", "optimized");

    // set a subsurface face to Dirichlet -- this is the bottom and should be
    // trivially handled.
    bc_markers[0] = Operators::MATRIX_BC_DIRICHLET;
    bc_values[0] = 1.;

    // set a surface BC edge to Dirichlet
    bc_surf_markers[0] = Operators::MATRIX_BC_DIRICHLET;
    bc_surf_values[0] = 100.;

    // -- TPF on surface
    Atpf = Teuchos::rcp(new Operators::MatrixMFD_TPFA(mfd_plist, surf_mesh));
    Atpf->SetSymmetryProperty(false);
    Atpf->SymbolicAssembleGlobalMatrices();
    Atpf->CreateMFDmassMatrices(Teuchos::null);
    Atpf->CreateMFDstiffnessMatrices(Teuchos::null);
    Atpf->CreateMFDrhsVectors();
    Atpf->ApplyBoundaryConditions(bc_surf_markers, bc_surf_values);
    Atpf->AssembleGlobalMatrices();

    // -- combined on domain
    As = Teuchos::rcp(new Operators::MatrixMFD_Surf(mfd_plist, mesh, surf_mesh));
    As->set_surface_A(Atpf);
    As->SetSymmetryProperty(false);
    As->SymbolicAssembleGlobalMatrices();
    As->CreateMFDmassMatrices(Teuchos::null);
    As->CreateMFDstiffnessMatrices(Teuchos::null);
    As->CreateMFDrhsVectors();
    As->ApplyBoundaryConditions(bc_surf_markers, bc_surf_values);
    As->AssembleGlobalMatrices();


    // dump the schur complement
    Teuchos::RCP<const Epetra_FECrsMatrix> Spp = Atpf->TPFA();
    EpetraExt::RowMatrixToMatlabFile("TPFA.txt", *Spp);
    Teuchos::RCP<const Epetra_FECrsMatrix> Aff = As->Aff();
    EpetraExt::RowMatrixToMatlabFile("Aff.txt", *Aff);

  }

};
