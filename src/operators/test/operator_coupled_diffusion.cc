/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Operators

  Tests nonlinear, coupled diffusion problem.
*/


#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "EpetraExt_RowMatrixOut.h"
#include "UnitTest++.h"

// Amanzi
#include "Mesh_MSTK.hh"
#include "Mesh_Algorithms.hh"
#include "Tensor.hh"
#include "WhetStoneDefs.hh"

// Operators
#include "BCs.hh"
#include "OperatorDefs.hh"
#include "PDE_DiffusionFactory.hh"
#include "TreeOperator.hh"

// test classes
#include "AnalyticNonlinearCoupled00.hh"
#include "AnalyticNonlinearCoupledBase.hh"

using namespace Amanzi;

double
BoundaryFaceGetValue(Operators::BCs& bc, const CompositeVector& u, int f)
{
  if (bc.bc_model()[f] == Operators::OPERATOR_BC_DIRICHLET) {
    return bc.bc_value()[f];
  } else {
    if (u.HasComponent("face")) {
      return u("face", 0, f);
    } else if (u.HasComponent("boundary_face")) {
      int bf = AmanziMesh::getFaceOnBoundaryBoundaryFace(*u.Mesh(), f);
      return u("boundary_face", 0, bf);
    } else {
      return u("cell", 0, AmanziMesh::getFaceOnBoundaryInternalCell(*u.Mesh(), f));
    }
  }
}


struct Problem {
 public:
  Problem(const Teuchos::RCP<AmanziMesh::Mesh>& mesh_,
          const Teuchos::RCP<AnalyticNonlinearCoupledBase>& ana_,
          const std::string& discretization_)
    : mesh(mesh_), ana(ana_), discretization(discretization_), comm(mesh_->getComm())
  {}

  ~Problem() {}

  void Setup() { MakeBCs(); }

  void FillCoefs(const CompositeVector& u, const CompositeVector& v)
  {
    int ncells = mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);
    int nfaces = mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::OWNED);

    // kr0,1 on faces, always upwinded
    Epetra_MultiVector& kr0_f = *kr0->ViewComponent("face", true);
    Epetra_MultiVector& kr1_f = *kr1->ViewComponent("face", true);

    const Epetra_MultiVector& u_c = *u.ViewComponent("cell", true);
    const Epetra_MultiVector& v_c = *v.ViewComponent("cell", true);

    for (int f = 0; f != nfaces; ++f) {
      AmanziMesh::Entity_ID_List cells;
      cells = mesh->getFaceCells(f, AmanziMesh::Parallel_type::ALL);

      if (cells.size() == 2) {
        if (u_c[0][cells[0]] > u_c[0][cells[1]]) {
          kr0_f[0][f] = ana->ScalarCoefficient00(u_c[0][cells[0]], v_c[0][cells[0]]);
        } else {
          kr0_f[0][f] = ana->ScalarCoefficient00(u_c[0][cells[1]], v_c[0][cells[1]]);
        }

        if (v_c[0][cells[0]] > v_c[0][cells[1]]) {
          kr1_f[0][f] = ana->ScalarCoefficient11(u_c[0][cells[0]], v_c[0][cells[0]]);
        } else {
          kr1_f[0][f] = ana->ScalarCoefficient11(u_c[0][cells[1]], v_c[0][cells[1]]);
        }
      } else {
        double u_boundary = BoundaryFaceGetValue(*bc0, u, f);
        double v_boundary = BoundaryFaceGetValue(*bc1, v, f);
        if (u_c[0][cells[0]] > u_boundary) {
          kr0_f[0][f] = ana->ScalarCoefficient00(u_c[0][cells[0]], v_c[0][cells[0]]);
        } else {
          kr0_f[0][f] = ana->ScalarCoefficient00(u_boundary, v_boundary);
        }

        if (v_c[0][cells[0]] > v_boundary) {
          kr1_f[0][f] = ana->ScalarCoefficient11(u_c[0][cells[0]], v_c[0][cells[0]]);
        } else {
          kr1_f[0][f] = ana->ScalarCoefficient11(u_boundary, v_boundary);
        }
      }
    }

    // derivatives
    if (kr0_u->HasComponent("cell")) {
      // FV, only need derivs on cells
      Epetra_MultiVector& kr0_u_c = *kr0_u->ViewComponent("cell", false);
      Epetra_MultiVector& kr0_v_c = *kr0_v->ViewComponent("cell", false);
      Epetra_MultiVector& kr1_u_c = *kr1_u->ViewComponent("cell", false);
      Epetra_MultiVector& kr1_v_c = *kr1_v->ViewComponent("cell", false);

      const Epetra_MultiVector& u_c = *u.ViewComponent("cell", false);
      const Epetra_MultiVector& v_c = *v.ViewComponent("cell", false);

      for (int c = 0; c != ncells; ++c) {
        kr0_u_c[0][c] = ana->DScalarCoefficient00D0(u_c[0][c], v_c[0][c]);
        kr0_v_c[0][c] = ana->DScalarCoefficient00D1(u_c[0][c], v_c[0][c]);
        kr1_u_c[0][c] = ana->DScalarCoefficient11D0(u_c[0][c], v_c[0][c]);
        kr1_v_c[0][c] = ana->DScalarCoefficient11D1(u_c[0][c], v_c[0][c]);
      }
    } else {
      // upwind derivatives
      Epetra_MultiVector& kr0_u_f = *kr0_u->ViewComponent("face", false);
      Epetra_MultiVector& kr0_v_f = *kr0_v->ViewComponent("face", false);
      Epetra_MultiVector& kr1_u_f = *kr1_u->ViewComponent("face", false);
      Epetra_MultiVector& kr1_v_f = *kr1_v->ViewComponent("face", false);

      for (int f = 0; f != nfaces; ++f) {
        AmanziMesh::Entity_ID_List cells;
        cells = mesh->getFaceCells(f, AmanziMesh::Parallel_type::ALL);

        if (cells.size() == 2) {
          if (u_c[0][cells[0]] > u_c[0][cells[1]]) {
            kr0_u_f[0][f] = ana->DScalarCoefficient00D0(u_c[0][cells[0]], v_c[0][cells[0]]);
            kr0_v_f[0][f] = ana->DScalarCoefficient00D1(u_c[0][cells[0]], v_c[0][cells[0]]);
          } else {
            kr0_u_f[0][f] = ana->DScalarCoefficient00D0(u_c[0][cells[1]], v_c[0][cells[1]]);
            kr0_v_f[0][f] = ana->DScalarCoefficient00D1(u_c[0][cells[1]], v_c[0][cells[1]]);
          }

          if (v_c[0][cells[0]] > v_c[0][cells[1]]) {
            kr1_u_f[0][f] = ana->DScalarCoefficient11D0(u_c[0][cells[0]], v_c[0][cells[0]]);
            kr1_v_f[0][f] = ana->DScalarCoefficient11D1(u_c[0][cells[0]], v_c[0][cells[0]]);
          } else {
            kr1_u_f[0][f] = ana->DScalarCoefficient11D0(u_c[0][cells[1]], v_c[0][cells[1]]);
            kr1_v_f[0][f] = ana->DScalarCoefficient11D1(u_c[0][cells[1]], v_c[0][cells[1]]);
          }

        } else {
          double u_boundary = BoundaryFaceGetValue(*bc0, u, f);
          double v_boundary = BoundaryFaceGetValue(*bc1, v, f);
          if (u_c[0][cells[0]] > u_boundary) {
            kr0_u_f[0][f] = ana->DScalarCoefficient00D0(u_c[0][cells[0]], v_c[0][cells[0]]);
            kr0_v_f[0][f] = ana->DScalarCoefficient00D1(u_c[0][cells[0]], v_c[0][cells[0]]);
          } else {
            kr0_u_f[0][f] = ana->DScalarCoefficient00D0(u_boundary, v_boundary);
            kr0_v_f[0][f] = ana->DScalarCoefficient00D1(u_boundary, v_boundary);
          }

          if (v_c[0][cells[0]] > v_boundary) {
            kr1_u_f[0][f] = ana->DScalarCoefficient11D0(u_c[0][cells[0]], v_c[0][cells[0]]);
            kr1_v_f[0][f] = ana->DScalarCoefficient11D1(u_c[0][cells[0]], v_c[0][cells[0]]);
          } else {
            kr1_u_f[0][f] = ana->DScalarCoefficient11D0(u_boundary, v_boundary);
            kr1_v_f[0][f] = ana->DScalarCoefficient11D1(u_boundary, v_boundary);
          }
        }
      }
    }

    kr0->ScatterMasterToGhosted();
    kr1->ScatterMasterToGhosted();
    kr0_u->ScatterMasterToGhosted();
    kr0_v->ScatterMasterToGhosted();
    kr1_u->ScatterMasterToGhosted();
    kr1_v->ScatterMasterToGhosted();
  }

  void MakeBCs()
  {
    bc0 = Teuchos::rcp(new Operators::BCs(mesh, AmanziMesh::Entity_kind::FACE, WhetStone::DOF_Type::SCALAR));
    std::vector<int>& bc_model0 = bc0->bc_model();
    std::vector<double>& bc_value0 = bc0->bc_value();

    bc1 = Teuchos::rcp(new Operators::BCs(mesh, AmanziMesh::Entity_kind::FACE, WhetStone::DOF_Type::SCALAR));
    std::vector<int>& bc_model1 = bc1->bc_model();
    std::vector<double>& bc_value1 = bc1->bc_value();

    int nfaces_wghost = mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::ALL);
    for (int f = 0; f < nfaces_wghost; f++) {
      const AmanziGeometry::Point& xf = mesh->getFaceCentroid(f);

      if (fabs(xf[0]) < 1e-12) {
        double area = mesh->getFaceArea(f);
        int dir = 0;
        int c = AmanziMesh::getFaceOnBoundaryInternalCell(*mesh, f);
        const AmanziGeometry::Point& normal = mesh->getFaceNormal(f, &dir);

        bc_model0[f] = Operators::OPERATOR_BC_NEUMANN;
        bc_value0[f] = ana->velocity_exact0(xf, 0.0) * normal / area;
        bc_model1[f] = Operators::OPERATOR_BC_NEUMANN;
        bc_value1[f] = ana->velocity_exact1(xf, 0.0) * normal / area;

      } else if (fabs(xf[1]) < 1e-12) {
        double area = mesh->getFaceArea(f);
        int dir = 0;
        int c = AmanziMesh::getFaceOnBoundaryInternalCell(*mesh, f);
        const AmanziGeometry::Point& normal = mesh->getFaceNormal(f, &dir);

        // y = 0 boundary
        bc_model0[f] = Operators::OPERATOR_BC_DIRICHLET;
        bc_value0[f] = ana->exact0(xf, 0.0);
        bc_model1[f] = Operators::OPERATOR_BC_NEUMANN;
        bc_value1[f] = ana->velocity_exact1(xf, 0.0) * normal / area;

      } else if (fabs(xf[0] - 1.0) < 1e-12) {
        double area = mesh->getFaceArea(f);
        int dir = 0;
        int c = AmanziMesh::getFaceOnBoundaryInternalCell(*mesh, f);
        const AmanziGeometry::Point& normal = mesh->getFaceNormal(f, &dir);

        // x = 1 boundary
        bc_model0[f] = Operators::OPERATOR_BC_NEUMANN;
        bc_value0[f] = ana->velocity_exact0(xf, 0.0) * normal / area;
        bc_model1[f] = Operators::OPERATOR_BC_DIRICHLET;
        bc_value1[f] = ana->exact1(xf, 0.0);

      } else if (fabs(xf[1] - 1.0) < 1e-12) {
        // y = 1 boudaries
        bc_model0[f] = Operators::OPERATOR_BC_DIRICHLET;
        bc_value0[f] = ana->exact0(xf, 0.0);
        bc_model1[f] = Operators::OPERATOR_BC_DIRICHLET;
        bc_value1[f] = ana->exact1(xf, 0.0);
      }
    }
  }

  void MakeTensorCoefs()
  {
    int ncells = mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);

    // tensor list
    K00 = Teuchos::rcp(new std::vector<WhetStone::Tensor>());
    for (int c = 0; c < ncells; c++) {
      const AmanziGeometry::Point& xc = mesh->getCellCentroid(c);
      const WhetStone::Tensor& Kc = ana->Tensor00(xc, 0.0);
      K00->push_back(Kc);
    }

    K11 = Teuchos::rcp(new std::vector<WhetStone::Tensor>());
    for (int c = 0; c < ncells; c++) {
      const AmanziGeometry::Point& xc = mesh->getCellCentroid(c);
      const WhetStone::Tensor& Kc = ana->Tensor11(xc, 0.0);
      K11->push_back(Kc);
    }
  }

  void MakeScalarTensorCoefs()
  {
    int ncells = mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);

    // tensor list
    K00 = Teuchos::rcp(new std::vector<WhetStone::Tensor>());
    K11 = Teuchos::rcp(new std::vector<WhetStone::Tensor>());

    for (int c = 0; c < ncells; c++) {
      const AmanziGeometry::Point& xc = mesh->getCellCentroid(c);

      double u = ana->exact0(xc, 0);
      double v = ana->exact1(xc, 0);

      WhetStone::Tensor Kc0 = ana->Tensor00(xc, 0.0);
      Kc0 *= ana->ScalarCoefficient00(u, v);
      K00->push_back(Kc0);

      WhetStone::Tensor Kc1 = ana->Tensor11(xc, 0.0);
      Kc1 *= ana->ScalarCoefficient11(u, v);
      K11->push_back(Kc1);
    }
  }

  void MakeScalarCoefSpace()
  {
    CompositeVectorSpace kr_space;
    kr_space.SetMesh(mesh)->SetGhosted(true)->SetComponent("face", AmanziMesh::Entity_kind::FACE, 1);

    kr0 = Teuchos::rcp(new CompositeVector(kr_space));
    kr1 = Teuchos::rcp(new CompositeVector(kr_space));

    if (discretization == "fv: default") {
      CompositeVectorSpace dkr_space;
      dkr_space.SetMesh(mesh)->SetGhosted(true)->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
      kr0_u = Teuchos::rcp(new CompositeVector(dkr_space));
      kr0_v = Teuchos::rcp(new CompositeVector(dkr_space));
      kr1_u = Teuchos::rcp(new CompositeVector(dkr_space));
      kr1_v = Teuchos::rcp(new CompositeVector(dkr_space));
    } else {
      kr0_u = Teuchos::rcp(new CompositeVector(kr_space));
      kr0_v = Teuchos::rcp(new CompositeVector(kr_space));
      kr1_u = Teuchos::rcp(new CompositeVector(kr_space));
      kr1_v = Teuchos::rcp(new CompositeVector(kr_space));
    }
  }

  void FillSolution(CompositeVector& u, CompositeVector& v)
  {
    int ncells = mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);
    Epetra_MultiVector& u_c = *u.ViewComponent("cell", false);
    Epetra_MultiVector& v_c = *v.ViewComponent("cell", false);

    for (int c = 0; c != ncells; ++c) {
      const AmanziGeometry::Point& xc = mesh->getCellCentroid(c);
      u_c[0][c] = ana->exact0(xc, 0);
      v_c[0][c] = ana->exact1(xc, 0);
    }

    if (u.HasComponent("face")) {
      int nfaces = mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::OWNED);
      Epetra_MultiVector& u_f = *u.ViewComponent("face", false);
      Epetra_MultiVector& v_f = *v.ViewComponent("face", false);

      for (int f = 0; f != nfaces; ++f) {
        const AmanziGeometry::Point& xf = mesh->getFaceCentroid(f);
        u_f[0][f] = ana->exact0(xf, 0);
        v_f[0][f] = ana->exact1(xf, 0);
      }
    }

    if (u.HasComponent("boundary_face")) {
      int nboundary_faces =
        mesh->getNumEntities(AmanziMesh::Entity_kind::BOUNDARY_FACE, AmanziMesh::Parallel_type::OWNED);
      Epetra_MultiVector& u_f = *u.ViewComponent("boundary_face", false);
      Epetra_MultiVector& v_f = *v.ViewComponent("boundary_face", false);

      for (int bf = 0; bf != nboundary_faces; ++bf) {
        int f = mesh->getMap(AmanziMesh::Entity_kind::FACE,false).LID(mesh->getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE,false).GID(bf));

        const AmanziGeometry::Point& xf = mesh->getFaceCentroid(f);
        u_f[0][f] = ana->exact0(xf, 0);
        v_f[0][f] = ana->exact1(xf, 0);
      }
    }
  }

  void CreateTreeVectorSpace()
  {
    tvs = Teuchos::rcp(new TreeVectorSpace());
    tvs->PushBack(
      Teuchos::rcp(new TreeVectorSpace(Teuchos::rcpFromRef(op00->global_operator()->DomainMap()))));
    tvs->PushBack(
      Teuchos::rcp(new TreeVectorSpace(Teuchos::rcpFromRef(op11->global_operator()->DomainMap()))));
  }

  void CreateOperator()
  {
    op = Teuchos::rcp(new Operators::TreeOperator(tvs));
    op->set_operator_block(0, 0, op00->global_operator());
    op->set_operator_block(1, 1, op11->global_operator());
  }

  void CreatePC()
  {
    pc = Teuchos::rcp(new Operators::TreeOperator(tvs));
    pc->set_operator_block(0, 0, pc00->global_operator());
    pc->set_operator_block(1, 1, pc11->global_operator());

    if (pc01 != Teuchos::null) {
      pc->set_operator_block(0, 1, pc01->global_operator());
      pc->set_operator_block(1, 0, pc10->global_operator());
    }
  }


  void CreateBlockOperators(bool upwind)
  {
    Operators::PDE_DiffusionFactory fac;
    Teuchos::RCP<const AmanziMesh::Mesh> mesh_c = mesh;

    Teuchos::ParameterList op_list00;
    op_list00.set("discretization primary", discretization);
    if (upwind) {
      op_list00.set("nonlinear coefficient", "upwind: face");
    } else {
      op_list00.set("nonlinear coefficient", "none");
    }
    op00 = fac.Create(op_list00, mesh_c);
    op00->SetBCs(bc0, bc0);
    op00->SetTensorCoefficient(K00);

    Teuchos::ParameterList op_list11;
    op_list11.set("discretization primary", discretization);
    if (upwind) {
      op_list11.set("nonlinear coefficient", "upwind: face");
    } else {
      op_list00.set("nonlinear coefficient", "none");
    }

    op11 = fac.Create(op_list11, mesh_c);
    op11->SetBCs(bc1, bc1);
    op11->SetTensorCoefficient(K11);
  }


  void CreateBlockPCs(bool jac_ondiag, bool jac_offdiag, bool upwind)
  {
    Operators::PDE_DiffusionFactory fac;
    Teuchos::RCP<const AmanziMesh::Mesh> mesh_c = mesh;

    Teuchos::ParameterList op_list00;
    op_list00.set("discretization primary", discretization);
    if (upwind) {
      op_list00.set("nonlinear coefficient", "upwind: face");
    } else {
      op_list00.set("nonlinear coefficient", "none");
    }
    if (jac_ondiag) {
      if (discretization == "fv: default") {
        op_list00.set("Newton correction", "true Jacobian");
      } else {
        op_list00.set("Newton correction", "approximate Jacobian");
      }
    }
    pc00 = fac.Create(op_list00, mesh_c);
    pc00->SetBCs(bc0, bc0);
    pc00->SetTensorCoefficient(K00);

    Teuchos::ParameterList op_list11;
    op_list11.set("discretization primary", discretization);
    if (upwind) {
      op_list11.set("nonlinear coefficient", "upwind: face");
    } else {
      op_list11.set("nonlinear coefficient", "none");
    }
    if (jac_ondiag) {
      if (discretization == "fv: default") {
        op_list11.set("Newton correction", "true Jacobian");
      } else {
        op_list11.set("Newton correction", "approximate Jacobian");
      }
    }
    pc11 = fac.Create(op_list11, mesh_c);
    pc11->SetBCs(bc1, bc1);
    pc11->SetTensorCoefficient(K11);


    if (jac_offdiag) {
      Teuchos::ParameterList op_list01;
      op_list01.set("discretization primary", discretization);
      if (upwind) {
        op_list01.set("nonlinear coefficient", "upwind: face");
      } else {
        op_list01.set("nonlinear coefficient", "none");
      }
      op_list01.set("exclude primary terms", true);
      if (discretization == "fv: default") {
        op_list01.set("Newton correction", "true Jacobian");
      } else {
        op_list01.set("Newton correction", "approximate Jacobian");
      }
      pc01 = fac.Create(op_list01, mesh_c);
      pc01->SetBCs(bc0, bc1);
      pc01->SetTensorCoefficient(K00);

      Teuchos::ParameterList op_list10;
      op_list10.set("discretization primary", discretization);
      if (upwind) {
        op_list10.set("nonlinear coefficient", "upwind: face");
      } else {
        op_list10.set("nonlinear coefficient", "none");
      }
      op_list10.set("exclude primary terms", true);
      if (discretization == "fv: default") {
        op_list10.set("Newton correction", "true Jacobian");
      } else {
        op_list10.set("Newton correction", "approximate Jacobian");
      }
      pc10 = fac.Create(op_list10, mesh_c);
      pc10->SetBCs(bc1, bc0);
      pc10->SetTensorCoefficient(K11);
    }
  }

  void CreateSources()
  {
    CompositeVectorSpace src_space;
    src_space.SetMesh(mesh)->SetGhosted(true)->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    f0 = Teuchos::rcp(new CompositeVector(src_space));
    f1 = Teuchos::rcp(new CompositeVector(src_space));

    int ncells = mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);
    Epetra_MultiVector& f0_c = *f0->ViewComponent("cell", false);
    Epetra_MultiVector& f1_c = *f1->ViewComponent("cell", false);

    for (int c = 0; c != ncells; ++c) {
      const AmanziGeometry::Point& xc = mesh->getCellCentroid(c);
      f0_c[0][c] = ana->source_exact0(xc, 0);
      f1_c[0][c] = ana->source_exact1(xc, 0);
    }
  }

  void CreateFluxes()
  {
    CompositeVectorSpace q_space;
    q_space.SetMesh(mesh)->SetGhosted(true)->SetComponent("face", AmanziMesh::Entity_kind::FACE, 1);
    q0 = Teuchos::rcp(new CompositeVector(q_space));
    q1 = Teuchos::rcp(new CompositeVector(q_space));
  }


  void Report(CompositeVector& u, CompositeVector& v)
  {
    std::cout << "Problem with "
              << mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::ALL) << " on "
              << mesh->getComm()->NumProc() << " cores with discretization \"" << discretization
              << "\"" << std::endl;
    std::cout << "Solution:" << std::endl;
    std::cout << "U = ";
    u.Print(std::cout);
    std::cout << "V = ";
    v.Print(std::cout);
    std::cout << "-------------------------" << std::endl;
    std::cout << "Kr0 = ";
    kr0->Print(std::cout);
    std::cout << "Kr1 = ";
    kr1->Print(std::cout);
    std::cout << "-------------------------" << std::endl;
  }

  void ReportError(const std::string& filename,
                   CompositeVector& err0,
                   CompositeVector& err1,
                   bool faces = false)
  {
    std::ofstream fid;
    fid.open(filename.c_str());

    int ncells = mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);
    Epetra_MultiVector& e0_c = *err0.ViewComponent("cell", false);
    Epetra_MultiVector& e1_c = *err1.ViewComponent("cell", false);

    fid.precision(14);
    fid << std::scientific;

    for (int c = 0; c != ncells; ++c) {
      const AmanziGeometry::Point& xc = mesh->getCellCentroid(c);
      fid << xc[0] << " " << xc[1] << " " << e0_c[0][c] << " " << e1_c[0][c] << std::endl;
    }

    if (faces && err0.HasComponent("face")) {
      int nfaces = mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::OWNED);
      Epetra_MultiVector& e0_f = *err0.ViewComponent("face", false);
      Epetra_MultiVector& e1_f = *err1.ViewComponent("face", false);

      for (int c = 0; c != nfaces; ++c) {
        const AmanziGeometry::Point& xc = mesh->getFaceCentroid(c);
        fid << xc[0] << " " << xc[1] << " " << e0_f[0][c] << " " << e1_f[0][c] << std::endl;
      }
    }
    fid.close();
  }

 public:
  Teuchos::RCP<Operators::BCs> bc0;
  Teuchos::RCP<Operators::BCs> bc1;

  Teuchos::RCP<std::vector<WhetStone::Tensor>> K00;
  Teuchos::RCP<std::vector<WhetStone::Tensor>> K11;

  Teuchos::RCP<CompositeVector> kr0;
  Teuchos::RCP<CompositeVector> kr1;
  Teuchos::RCP<CompositeVector> kr0_u;
  Teuchos::RCP<CompositeVector> kr0_v;
  Teuchos::RCP<CompositeVector> kr1_u;
  Teuchos::RCP<CompositeVector> kr1_v;

  Teuchos::RCP<CompositeVector> f0;
  Teuchos::RCP<CompositeVector> f1;

  Teuchos::RCP<CompositeVector> q0;
  Teuchos::RCP<CompositeVector> q1;

  Teuchos::RCP<Operators::PDE_Diffusion> op00;
  Teuchos::RCP<Operators::PDE_Diffusion> op11;

  Teuchos::RCP<Operators::PDE_Diffusion> pc00;
  Teuchos::RCP<Operators::PDE_Diffusion> pc11;
  Teuchos::RCP<Operators::PDE_Diffusion> pc01;
  Teuchos::RCP<Operators::PDE_Diffusion> pc10;

  Teuchos::RCP<Operators::TreeOperator> op;
  Teuchos::RCP<Operators::TreeOperator> pc;

  Teuchos::RCP<TreeVectorSpace> tvs;

  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh;
  Teuchos::RCP<AnalyticNonlinearCoupledBase> ana;
  std::string discretization;

  Comm_ptr_type comm;
};


Teuchos::RCP<Problem>
getProblem(const std::string& discretization, bool upwind, int nx, int ny)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();

  // create a mesh
  auto mesh_mstk = Teuchos::rcp(new Mesh_MSTK(0., 0., 1., 1., nx, ny, comm));
  auto mesh = Teuchos::rcp(new Mesh(mesh_mstk, Teuchos::null)); 

  // create the analytic solution
  Teuchos::RCP<AnalyticNonlinearCoupledBase> ana =
    Teuchos::rcp(new AnalyticNonlinearCoupled00(mesh));

  // create the problem
  Teuchos::RCP<Problem> problem = Teuchos::rcp(new Problem(mesh, ana, discretization));
  problem->MakeBCs();
  if (upwind) {
    problem->MakeTensorCoefs();
    problem->MakeScalarCoefSpace();
  } else {
    problem->MakeScalarTensorCoefs();
  }
  return problem;
}


// -----------------------------------------------------------------------------
// Creates the linear operator with the correct coefficients, checks the error
// in the forward application of the Operator to the true solution relative to
// the rhs.
//
// Returns the l2, linf errors as a pair.
// -----------------------------------------------------------------------------
std::pair<double, double>
RunForwardProblem(const std::string& discretization, bool upwind, int nx, int ny)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Teuchos::RCP<Problem> problem = getProblem(discretization, upwind, nx, ny);

  // test the forward solution
  problem->CreateBlockOperators(upwind);
  problem->CreateTreeVectorSpace();
  problem->CreateOperator();
  problem->CreateSources();

  // get u,v
  Teuchos::RCP<CompositeVector> u =
    Teuchos::rcp(new CompositeVector(problem->op00->global_operator()->DomainMap()));
  Teuchos::RCP<CompositeVector> v =
    Teuchos::rcp(new CompositeVector(problem->op11->global_operator()->DomainMap()));
  problem->FillSolution(*u, *v);

  // get coefficients, fill matrices
  if (upwind) problem->FillCoefs(*u, *v);
  problem->op00->SetScalarCoefficient(problem->kr0, Teuchos::null);
  problem->op00->UpdateMatrices(Teuchos::null, Teuchos::null);
  problem->op00->global_operator()->UpdateRHS(*problem->f0, false);
  problem->op00->ApplyBCs(true, true, true);

  problem->op11->SetScalarCoefficient(problem->kr1, Teuchos::null);
  problem->op11->UpdateMatrices(Teuchos::null, Teuchos::null);
  problem->op11->global_operator()->UpdateRHS(*problem->f1, false);
  problem->op11->ApplyBCs(true, true, true);

  // apply!
  TreeVector B(*problem->tvs);
  *B.SubVector(0)->Data() = *problem->op00->global_operator()->rhs();
  *B.SubVector(1)->Data() = *problem->op11->global_operator()->rhs();

  TreeVector X(*problem->tvs);
  *X.SubVector(0)->Data() = *u;
  *X.SubVector(1)->Data() = *v;

  TreeVector AX(*problem->tvs);
  problem->op->Apply(X, AX);
  AX.Update(-1., B, 1.);

  double error0_l2, error1_l2;
  double error0_linf, error1_linf;
  double unorm, vnorm;

  B.SubVector(0)->Data()->ViewComponent("cell")->Norm2(&unorm);
  B.SubVector(0)->Data()->ViewComponent("cell")->Norm2(&vnorm);

  AX.SubVector(0)->Data()->ViewComponent("cell")->Norm2(&error0_l2);
  AX.SubVector(1)->Data()->ViewComponent("cell")->Norm2(&error1_l2);
  AX.SubVector(0)->Data()->ViewComponent("cell")->NormInf(&error0_linf);
  AX.SubVector(1)->Data()->ViewComponent("cell")->NormInf(&error1_linf);

  double error_l2 =
    sqrt((error0_l2 * error0_l2 + error1_l2 * error1_l2) / (unorm * unorm + vnorm * vnorm));
  double error_linf = std::max(error0_linf, error1_linf);

  if (problem->comm->MyPID() == 0) {
    printf("[%4d, %6.12e, %6.12e],\n", (int)round(log2(nx)), log2(error_l2), log2(error_linf));
  }

  return std::make_pair(log2(error_l2), log2(error_linf));
}


// -----------------------------------------------------------------------------
// Creates the linear operator with the correct coefficients, checks the error
// in the forward application of the Operator to the true solution relative to
// the rhs.
//
// Returns the l2, linf errors as a pair.
// -----------------------------------------------------------------------------
std::pair<double, double>
RunForwardProblem_Assembled(const std::string& discretization, bool upwind, int nx, int ny)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Teuchos::RCP<Problem> problem = getProblem(discretization, upwind, nx, ny);

  // test the forward solution
  problem->CreateBlockOperators(upwind);
  problem->CreateTreeVectorSpace();
  problem->CreateOperator();
  problem->CreateSources();

  // get u,v
  Teuchos::RCP<CompositeVector> u =
    Teuchos::rcp(new CompositeVector(problem->op00->global_operator()->DomainMap()));
  Teuchos::RCP<CompositeVector> v =
    Teuchos::rcp(new CompositeVector(problem->op11->global_operator()->DomainMap()));
  problem->FillSolution(*u, *v);

  // get coefficients, fill matrices
  if (upwind) problem->FillCoefs(*u, *v);
  problem->op00->SetScalarCoefficient(problem->kr0, Teuchos::null);
  problem->op00->UpdateMatrices(Teuchos::null, Teuchos::null);
  problem->op00->global_operator()->UpdateRHS(*problem->f0, false);
  problem->op00->ApplyBCs(true, true, true);

  problem->op11->SetScalarCoefficient(problem->kr1, Teuchos::null);
  problem->op11->UpdateMatrices(Teuchos::null, Teuchos::null);
  problem->op11->global_operator()->UpdateRHS(*problem->f1, false);
  problem->op11->ApplyBCs(true, true, true);

  // create vector storage
  TreeVector B(*problem->tvs);
  *B.SubVector(0)->Data() = *problem->op00->global_operator()->rhs();
  *B.SubVector(1)->Data() = *problem->op11->global_operator()->rhs();

  TreeVector X(*problem->tvs);
  *X.SubVector(0)->Data() = *u;
  *X.SubVector(1)->Data() = *v;

  TreeVector AX(*problem->tvs);

  // apply
  problem->op->SymbolicAssembleMatrix();
  problem->op->AssembleMatrix();
  problem->op->ApplyAssembled(X, AX, 0.0);

  // subtract off true RHS
  AX.Update(-1., B, 1.);

  // norms
  double error_l2, error_linf;
  double unorm;

  B.Norm2(&unorm);
  AX.Norm2(&error_l2);
  AX.NormInf(&error_linf);

  error_l2 /= unorm;

  if (problem->comm->MyPID() == 0) {
    printf("[%4d, %6.12e, %6.12e],\n", (int)round(log2(nx)), log2(error_l2), log2(error_linf));
  }

  return std::make_pair(log2(error_l2), log2(error_linf));
}


// -----------------------------------------------------------------------------
// Creates the linear operator with the correct coefficients, checks the error
// in the inverse application of the Operator to the rhs relative to the true
// solution.
//
// Returns the l2, linf errors as a pair.
// -----------------------------------------------------------------------------
std::pair<double, double>
RunInverseProblem(const std::string& discretization, bool upwind, int nx, int ny, bool write_file)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Teuchos::RCP<Problem> problem = getProblem(discretization, upwind, nx, ny);

  // test the inverse problem
  problem->CreateBlockOperators(upwind);
  problem->CreateTreeVectorSpace();
  problem->CreateOperator();
  problem->CreateSources();

  // get u,v
  Teuchos::RCP<CompositeVector> u =
    Teuchos::rcp(new CompositeVector(problem->op00->global_operator()->DomainMap()));
  Teuchos::RCP<CompositeVector> v =
    Teuchos::rcp(new CompositeVector(problem->op11->global_operator()->DomainMap()));
  problem->FillSolution(*u, *v);

  // get coefficients, fill matrices
  if (upwind) problem->FillCoefs(*u, *v);
  problem->op00->SetScalarCoefficient(problem->kr0, Teuchos::null);
  problem->op00->UpdateMatrices(Teuchos::null, Teuchos::null);
  problem->op00->global_operator()->UpdateRHS(*problem->f0, false);
  problem->op00->ApplyBCs(true, true, true);

  problem->op11->SetScalarCoefficient(problem->kr1, Teuchos::null);
  problem->op11->UpdateMatrices(Teuchos::null, Teuchos::null);
  problem->op11->global_operator()->UpdateRHS(*problem->f1, false);
  problem->op11->ApplyBCs(true, true, true);

  // create vector storage
  TreeVector B(*problem->tvs);
  *B.SubVector(0)->Data() = *problem->op00->global_operator()->rhs();
  *B.SubVector(1)->Data() = *problem->op11->global_operator()->rhs();

  TreeVector X(*problem->tvs);
  TreeVector AX(*problem->tvs);

  problem->op->SymbolicAssembleMatrix();
  problem->op->AssembleMatrix();

  if (write_file) {
    std::stringstream fname;
    fname << "matrix_" << nx;
    if (discretization == "fv: default")
      fname << "_fv";
    else
      fname << "_mfd";
    fname << ".dat";
    EpetraExt::RowMatrixToMatlabFile(fname.str().c_str(), *problem->op->A());
  }


  Teuchos::ParameterList pc_list;
  pc_list.set("preconditioning method", "boomer amg");
  pc_list.sublist("boomer amg parameters").set("tolerance", 0.0);
  pc_list.sublist("boomer amg parameters").set("verbosity", 0);
  //  pc_list.sublist("boomer amg parameters").set("number of functions", 2);
  pc_list.set("iterative method", "pcg");
  pc_list.sublist("pcg parameters").set("maximum number of iterations", 200);
  pc_list.sublist("verbose object").set("verbosity level", "medium");
  problem->op->set_inverse_parameters(pc_list);
  problem->op->InitializeInverse();
  problem->op->ComputeInverse();

  X.PutScalar(0.);
  int ierr = problem->op->ApplyInverse(B, X);
  CHECK(ierr >= 0);
  CHECK(problem->op->num_itrs() < 100);

  // subtract off true solution
  X.SubVector(0)->Data()->Update(-1., *u, 1.);
  X.SubVector(1)->Data()->Update(-1., *v, 1.);

  // write error if requested
  if (write_file) {
    std::stringstream fname;
    fname << "error_" << nx;
    if (discretization == "fv: default")
      fname << "_fv";
    else
      fname << "_mfd";
    fname << ".dat";
    problem->ReportError(fname.str(), *X.SubVector(0)->Data(), *X.SubVector(1)->Data());
  }

  // norms
  double error0_l2, error1_l2;
  double error0_linf, error1_linf;
  double unorm, vnorm;

  u->ViewComponent("cell", false)->Norm2(&unorm);
  v->ViewComponent("cell", false)->Norm2(&vnorm);

  X.SubVector(0)->Data()->ViewComponent("cell", false)->Norm2(&error0_l2);
  X.SubVector(1)->Data()->ViewComponent("cell", false)->Norm2(&error1_l2);
  X.SubVector(0)->Data()->ViewComponent("cell", false)->NormInf(&error0_linf);
  X.SubVector(1)->Data()->ViewComponent("cell", false)->NormInf(&error1_linf);

  double error_l2 = sqrt(pow(error0_l2 / unorm, 2) + pow(error1_l2 / vnorm, 2));
  double error_linf = std::max(error0_linf, error1_linf);

  if (problem->comm->MyPID() == 0) {
    printf("[%4d, %6.12e, %6.12e],\n", (int)round(log2(nx)), log2(error_l2), log2(error_linf));
  }

  return std::make_pair(log2(error_l2), log2(error_linf));
}


// -----------------------------------------------------------------------------
// Solves the full nonlinear problem
//
// Returns the l2, linf errors as a pair.
// -----------------------------------------------------------------------------
std::pair<double, double>
RunNonlinearProblem(const std::string& discretization,
                    const std::string& jacobian,
                    bool upwind,
                    int nx,
                    int ny,
                    bool write_file,
                    double damping)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Teuchos::RCP<Problem> problem = getProblem(discretization, upwind, nx, ny);

  // test the forward solution
  problem->CreateBlockOperators(upwind);
  if (jacobian == "none") {
    problem->CreateBlockPCs(false, false, upwind);
  } else if (jacobian == "diagonal") {
    problem->CreateBlockPCs(true, false, upwind);
  } else if (jacobian == "full") {
    problem->CreateBlockPCs(true, true, upwind);
  } else {
    AMANZI_ASSERT(0);
  }

  problem->CreateTreeVectorSpace();
  problem->CreateOperator();
  if (write_file) problem->op->SymbolicAssembleMatrix(); // debug, write assembled matrix
  problem->CreatePC();
  problem->CreateSources();
  problem->CreateFluxes();

  // simple newton iteration
  TreeVector X(*problem->tvs);
  X.PutScalar(0.2);
  Teuchos::RCP<CompositeVector> u = X.SubVector(0)->Data();
  Teuchos::RCP<CompositeVector> v = X.SubVector(1)->Data();
  TreeVector R(*problem->tvs);
  TreeVector DX(*problem->tvs);

  double error = 1.;
  double l2error = 1.;
  int nits = 0;
  while (error > 1.e-8 && nits <= 100) {
    // update the data
    if (upwind) problem->FillCoefs(*u, *v);

    // assemble the forward operator
    problem->op00->global_operator()->Init();
    problem->op00->SetScalarCoefficient(problem->kr0, Teuchos::null);
    problem->op00->UpdateMatrices(Teuchos::null, Teuchos::null);
    problem->op00->global_operator()->UpdateRHS(*problem->f0, false);
    problem->op00->ApplyBCs(true, true, true);
    problem->op11->global_operator()->Init();
    problem->op11->SetScalarCoefficient(problem->kr1, Teuchos::null);
    problem->op11->UpdateMatrices(Teuchos::null, Teuchos::null);
    problem->op11->global_operator()->UpdateRHS(*problem->f1, false);
    problem->op11->ApplyBCs(true, true, true);

    // write forward operator to file for debugging
    if (write_file) {
      problem->op->AssembleMatrix();

      std::stringstream fname;
      fname << "files/operator_" << nits << "_" << nx;
      if (discretization == "fv: default")
        fname << "_fv";
      else
        fname << "_mfd";
      fname << ".dat";
      EpetraExt::RowMatrixToMatlabFile(fname.str().c_str(), *problem->op->A());
    }

    // apply the forward operator to get the residual
    // -- NOTE: update the interface to improve how we access TreeOp's RHS
    problem->op->Apply(X, R);
    R.SubVector(0)->Data()->Update(-1., *problem->op00->global_operator()->rhs(), 1);
    R.SubVector(1)->Data()->Update(-1., *problem->op11->global_operator()->rhs(), 1);

    R.NormInf(&error);
    R.Norm2(&l2error);
    double unorm;
    X.Norm2(&unorm);
    // std::cout << "Iteration: " << nits << " residual: Linf = " << error
    //           << ", L2 = " << l2error / unorm << std::endl;

    if (write_file) {
      std::stringstream suffix;
      suffix << nits << "_" << nx;
      if (discretization == "fv: default")
        suffix << "_fv";
      else
        suffix << "_mfd";
      suffix << ".dat";

      std::string fname = "files/residual" + suffix.str();
      problem->ReportError(fname, *R.SubVector(0)->Data(), *R.SubVector(1)->Data(), true);

      fname = "files/solution" + suffix.str();
      problem->ReportError(fname, *X.SubVector(0)->Data(), *X.SubVector(1)->Data(), true);

      if (nits > 0) {
        fname = "files/correction" + suffix.str();
        problem->ReportError(fname, *DX.SubVector(0)->Data(), *DX.SubVector(1)->Data(), true);
      }
    }


    if (error <= 1.e-8) continue;

    // update fluxes
    problem->op00->UpdateFlux(u.ptr(), problem->q0.ptr());
    problem->op11->UpdateFlux(v.ptr(), problem->q1.ptr());

    // assemble the preconditioner
    problem->pc00->global_operator()->Init();
    problem->pc00->SetScalarCoefficient(problem->kr0, problem->kr0_u);
    problem->pc00->UpdateMatrices(problem->q0.ptr(), u.ptr());
    problem->pc00->UpdateMatricesNewtonCorrection(problem->q0.ptr(), u.ptr());
    problem->pc00->ApplyBCs(true, true, true);

    problem->pc11->global_operator()->Init();
    problem->pc11->SetScalarCoefficient(problem->kr1, problem->kr1_v);
    problem->pc11->UpdateMatrices(problem->q1.ptr(), v.ptr());
    problem->pc11->UpdateMatricesNewtonCorrection(problem->q1.ptr(), v.ptr());
    problem->pc11->ApplyBCs(true, true, true);

    if (problem->pc01 != Teuchos::null) {
      problem->pc01->global_operator()->Init();
      problem->pc01->SetScalarCoefficient(problem->kr0, problem->kr0_v);
      problem->pc01->UpdateMatrices(problem->q0.ptr(), u.ptr());
      problem->pc01->UpdateMatricesNewtonCorrection(problem->q0.ptr(), u.ptr());
      problem->pc01->ApplyBCs(false, true, false);

      problem->pc10->global_operator()->Init();
      problem->pc10->SetScalarCoefficient(problem->kr1, problem->kr1_u);
      problem->pc10->UpdateMatrices(problem->q1.ptr(), v.ptr());
      problem->pc10->UpdateMatricesNewtonCorrection(problem->q1.ptr(), v.ptr());
      problem->pc10->ApplyBCs(false, true, false);
    }
    problem->pc->AssembleMatrix();

    // write preconditioner to file for debugging
    if (write_file) {
      std::stringstream fname;
      fname << "files/preconditioner_" << nits << "_" << nx;
      if (discretization == "fv: default")
        fname << "_fv";
      else
        fname << "_mfd";
      fname << ".dat";
      EpetraExt::RowMatrixToMatlabFile(fname.str().c_str(), *problem->pc->A());
    }

    Teuchos::ParameterList pc_list;
    pc_list.set("preconditioning method", "boomer amg");
    pc_list.sublist("boomer amg parameters").set("tolerance", 0.0);
    pc_list.sublist("boomer amg parameters").set("number of functions", 2);
    pc_list.set("iterative method", "gmres");
    pc_list.sublist("verbose object").set("verbosity level", "medium");
    problem->op->set_inverse_parameters(pc_list);
    problem->op->InitializeInverse();
    problem->op->ComputeInverse();

    // invert the preconditioner to get a correction
    DX.PutScalar(0.);
    int converged_reason = problem->op->ApplyInverse(R, DX);
    CHECK(converged_reason >= 0);

    // apply the correction
    X.Update(-damping, DX, 1.);
    nits++;
  }

  // check convergence
  CHECK(nits <= 100);

  // check error
  TreeVector U(X);
  problem->FillSolution(*U.SubVector(0)->Data(), *U.SubVector(1)->Data());
  X.Update(-1., U, 1.);

  if (write_file) {
    std::stringstream fname;
    fname << "error_" << nx;
    if (discretization == "fv: default")
      fname << "_fv";
    else
      fname << "_mfd";
    fname << ".dat";
    problem->ReportError(fname.str(), *X.SubVector(0)->Data(), *X.SubVector(1)->Data());
  }

  double error0_l2, error1_l2;
  double error0_linf, error1_linf;
  double unorm, vnorm;

  U.SubVector(0)->Data()->ViewComponent("cell", false)->Norm2(&unorm);
  U.SubVector(0)->Data()->ViewComponent("cell", false)->Norm2(&vnorm);

  X.SubVector(0)->Data()->ViewComponent("cell", false)->Norm2(&error0_l2);
  X.SubVector(1)->Data()->ViewComponent("cell", false)->Norm2(&error1_l2);
  X.SubVector(0)->Data()->ViewComponent("cell", false)->NormInf(&error0_linf);
  X.SubVector(1)->Data()->ViewComponent("cell", false)->NormInf(&error1_linf);

  double error_l2 = sqrt(pow(error0_l2 / unorm, 2) + pow(error1_l2 / vnorm, 2));
  double error_linf = std::max(error0_linf, error1_linf);

  if (problem->comm->MyPID() == 0) {
    printf("[%4d, %6.12e, %6.12e],\n", (int)round(log2(nx)), log2(error_l2), log2(error_linf));
  }

  return std::make_pair(log2(error_l2), log2(error_linf));
}


// -----------------------------------------------------------------------------
// Creates the linear operator with the correct coefficients, checks the error
// in the inverse application of the Operator to the rhs relative to the true
// solution.  Uses the block-diagonal operator by inverting each system
// independently.  This is more of an internal debugging check of the Analytic
// Solution than an actual test of the Operators class.
//
// Returns the l2, linf errors as a pair.
// -----------------------------------------------------------------------------
std::pair<double, double>
RunInverseProblem_Diag(const std::string& discretization,
                       bool upwind,
                       int nx,
                       int ny,
                       bool write_file)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Teuchos::RCP<Problem> problem = getProblem(discretization, upwind, nx, ny);

  // test the forward solution
  problem->CreateBlockOperators(upwind);
  problem->CreateTreeVectorSpace();
  problem->CreateOperator();
  problem->CreateSources();

  // get u,v
  Teuchos::RCP<CompositeVector> u =
    Teuchos::rcp(new CompositeVector(problem->op00->global_operator()->DomainMap()));
  Teuchos::RCP<CompositeVector> v =
    Teuchos::rcp(new CompositeVector(problem->op11->global_operator()->DomainMap()));
  problem->FillSolution(*u, *v);

  // get coefficients, fill matrices
  //  problem->FillCoefs(*u,*v);
  problem->op00->SetScalarCoefficient(problem->kr0, Teuchos::null);
  problem->op00->UpdateMatrices(Teuchos::null, Teuchos::null);
  problem->op00->global_operator()->UpdateRHS(*problem->f0, false);
  problem->op00->ApplyBCs(true, true, true);

  problem->op11->SetScalarCoefficient(problem->kr1, Teuchos::null);
  problem->op11->UpdateMatrices(Teuchos::null, Teuchos::null);
  problem->op11->global_operator()->UpdateRHS(*problem->f1, false);
  problem->op11->ApplyBCs(true, true, true);

  // assemble, invert
  Teuchos::ParameterList pc_list;
  pc_list.set("preconditioning method", "boomer amg");
  pc_list.set("iterative method", "gmres");
  pc_list.sublist("boomer amg parameters").set("tolerance", 0.);
  pc_list.sublist("boomer amg parameters").set("number of functions", 2);
  problem->op00->global_operator()->set_inverse_parameters(pc_list);
  problem->op00->global_operator()->InitializeInverse();
  problem->op00->global_operator()->ComputeInverse();

  Teuchos::ParameterList pc_list11;
  pc_list11.set("preconditioning method", "boomer amg");
  pc_list11.set("iterative method", "gmres");
  pc_list11.sublist("boomer amg parameters").set("tolerance", 0.);
  pc_list11.sublist("boomer amg parameters").set("number of functions", 2);
  problem->op11->global_operator()->set_inverse_parameters(pc_list11);
  problem->op11->global_operator()->InitializeInverse();
  problem->op11->global_operator()->ComputeInverse();

  TreeVector B(*problem->tvs);
  *B.SubVector(0)->Data() = *problem->op00->global_operator()->rhs();
  *B.SubVector(1)->Data() = *problem->op11->global_operator()->rhs();

  TreeVector X(*problem->tvs);
  X.PutScalar(0.);

  int ierr = problem->op00->global_operator()->ApplyInverse(*B.SubVector(0)->Data(),
                                                            *X.SubVector(0)->Data());
  CHECK(ierr >= 0);

  ierr = problem->op11->global_operator()->ApplyInverse(*B.SubVector(1)->Data(),
                                                        *X.SubVector(1)->Data());
  CHECK(ierr >= 0);

  //  std::cout << "Ran inverse with nitrs = " << lin_op->num_itrs() << std::endl;
  X.SubVector(0)->Data()->Update(-1., *u, 1.);
  X.SubVector(1)->Data()->Update(-1., *v, 1.);

  if (write_file) {
    std::stringstream fname;
    fname << "error_" << nx;
    if (discretization == "fv: default")
      fname << "_fv";
    else
      fname << "_mfd";
    fname << ".dat";
    problem->ReportError(fname.str(), *X.SubVector(0)->Data(), *X.SubVector(1)->Data());
  }

  double error0_l2, error1_l2;
  double error0_linf, error1_linf;
  double unorm, vnorm;

  u->ViewComponent("cell", false)->Norm2(&unorm);
  v->ViewComponent("cell", false)->Norm2(&vnorm);

  X.SubVector(0)->Data()->ViewComponent("cell", false)->Norm2(&error0_l2);
  X.SubVector(1)->Data()->ViewComponent("cell", false)->Norm2(&error1_l2);
  X.SubVector(0)->Data()->ViewComponent("cell", false)->NormInf(&error0_linf);
  X.SubVector(1)->Data()->ViewComponent("cell", false)->NormInf(&error1_linf);

  double error_l2 = sqrt(pow(error0_l2 / unorm, 2) + pow(error1_l2 / vnorm, 2));
  double error_linf = std::max(error0_linf, error1_linf);

  if (problem->comm->MyPID() == 0) {
    printf("[%4d, %6.12e, %6.12e],\n", (int)round(log2(nx)), log2(error_l2), log2(error_linf));
  }

  return std::make_pair(log2(error_l2), log2(error_linf));
}


void
RunForwardTest(const std::string& discretization, bool upwind)
{
  std::cout << std::endl
            << std::endl
            << "============================================================================="
            << std::endl;
  std::cout << "Convergence for forward operator with discretization: " << discretization;
  if (upwind) std::cout << " with upwinding";
  std::cout << std::endl;
  std::cout << "x = np.array([";
  std::vector<std::pair<double, double>> l2s;
  for (int i = 2; i <= 129; i *= 2) {
    std::pair<double, double> l2 = RunForwardProblem_Assembled(discretization, upwind, i, i);
    l2s.push_back(l2);
  }
  std::cout << "])" << std::endl;

  double mean_dl2 = 0.;
  double mean_dlinf = 0.;
  int size = 0;
  std::cout << "  Convergence rates = ";
  // mean via the last few -- let's not base convergence on a 2x2 problem!
  for (int i = 3; i != l2s.size(); ++i) {
    mean_dl2 += (l2s[i].first - l2s[i - 1].first);
    mean_dlinf += (l2s[i].second - l2s[i - 1].second);
    size++;
    std::cout << " (" << l2s[i].first - l2s[i - 1].first << ", "
              << (l2s[i].second - l2s[i - 1].second) << "),";
  }
  std::cout << std::endl;

  double rate2 = -mean_dl2 / size;
  double rateinf = -mean_dlinf / size;
  std::cout << " Mean convergence rate (l2, linf) = " << rate2 << ", " << rateinf << std::endl;

  if (upwind) {
    CHECK(rate2 > 0.8);
    CHECK(rateinf > 0.8);
  } else {
    CHECK(rate2 > 1.5);
    CHECK(rateinf > 1.5);
  }
}


void
RunInverseTest(const std::string& discretization, bool upwind)
{
  std::cout << std::endl
            << std::endl
            << "============================================================================="
            << std::endl;
  std::cout << "Convergence for inverse operator with discretization: " << discretization;
  if (upwind) std::cout << " with upwinding";
  std::cout << std::endl;

  std::cout << "x = np.array([\n";
  std::vector<std::pair<double, double>> l2s;
  for (int i = 4; i <= 128; i *= 2) {
    std::pair<double, double> l2 = RunInverseProblem(discretization, upwind, i, i, false);
    l2s.push_back(l2);
  }
  std::cout << "])" << std::endl;

  double mean_dl2 = 0.;
  double mean_dlinf = 0.;
  int size = 0;
  std::cout << "  Convergence rates = ";

  // mean via the last few -- let's not base convergence on a 2x2 problem!
  for (int i = 1; i != l2s.size(); ++i) {
    mean_dl2 += (l2s[i].first - l2s[i - 1].first);
    mean_dlinf += (l2s[i].second - l2s[i - 1].second);
    size++;
  }
  std::cout << std::endl;

  double rate2 = -mean_dl2 / size;
  double rateinf = -mean_dlinf / size;
  std::cout << " Mean convergence rate (l2, linf) = " << rate2 << ", " << rateinf << std::endl;

  if (upwind) {
    CHECK(rate2 > 0.8);
    CHECK(rateinf > 0.8);
  } else {
    CHECK(rate2 > 1.8);
    CHECK(rateinf > 1.5);
  }
}


void
RunNonlinearTest(const std::string& discretization, const std::string& jacobian)
{
  std::cout << std::endl
            << std::endl
            << "============================================================================="
            << std::endl;
  std::cout << "Convergence for nonlinear operator with discretization: " << discretization
            << " with jacobian: " << jacobian << std::endl;

  std::cout << "x = np.array([";
  std::vector<std::pair<double, double>> l2s;
  for (int i = 2; i <= 65; i *= 2) {
    std::pair<double, double> l2 =
      RunNonlinearProblem(discretization, jacobian, true, i, i, false, 0.5);
    l2s.push_back(l2);
  }
  std::cout << "])" << std::endl;

  double mean_dl2 = 0.;
  double mean_dlinf = 0.;
  int size = 0;
  std::cout << "  Convergence rates = ";

  // mean via the last few -- let's not base convergence on a 2x2 problem!
  for (int i = 3; i != l2s.size(); ++i) {
    mean_dl2 += (l2s[i].first - l2s[i - 1].first);
    mean_dlinf += (l2s[i].second - l2s[i - 1].second);
    std::cout << " (" << l2s[i].first - l2s[i - 1].first << ", "
              << (l2s[i].second - l2s[i - 1].second) << "),";
    size++;
  }
  std::cout << std::endl;

  double rate2 = -mean_dl2 / size;
  double rateinf = -mean_dlinf / size;
  std::cout << " Mean convergence rate (l2, linf) = " << rate2 << ", " << rateinf << std::endl;

  CHECK(rate2 > 0.8);
  CHECK(rateinf > 0.8);
}

// -----------------------------------------------------------------------------
// Full Suite is the following 12 tests
// -----------------------------------------------------------------------------
TEST(OPERATOR_COUPLED_DIFFUSION_FORWARD_HARMONIC_CONVERGENCE_MFD)
{
  RunForwardTest("mfd: default", false);
}
TEST(OPERATOR_COUPLED_DIFFUSION_FORWARD_HARMONIC_CONVERGENCE_FV)
{
  RunForwardTest("fv: default", false);
}

// TEST(OPERATOR_COUPLED_DIFFUSION_FORWARD_UPWIND_CONVERGENCE_MFD) {
//   RunForwardTest("mfd: default", true);
// }

TEST(OPERATOR_COUPLED_DIFFUSION_INVERSE_HARMONIC_CONVERGENCE_MFD)
{
  RunInverseTest("mfd: default", false);
}
TEST(OPERATOR_COUPLED_DIFFUSION_INVERSE_HARMONIC_CONVERGENCE_FV)
{
  RunInverseTest("fv: default", false);
}

TEST(OPERATOR_COUPLED_DIFFUSION_INVERSE_UPWIND_CONVERGENCE_MFD)
{
  RunInverseTest("mfd: default", true);
}
TEST(OPERATOR_COUPLED_DIFFUSION_INVERSE_UPWIND_CONVERGENCE_FV)
{
  RunInverseTest("fv: default", true);
}

// TEST(OPERATOR_COUPLED_DIFFUSION_NONLINEAR_UPWIND_CONVERGENCE_FV) {
//   RunNonlinearTest("fv: default", "none");
// }
// TEST(OPERATOR_COUPLED_DIFFUSION_NONLINEAR_UPWIND_CONVERGENCE_FV_JACOBIAN) {
//   RunNonlinearTest("fv: default", "full");
// }


// TEST(OPERATOR_COUPLED_DIFFUSION_NONLINEAR_UPWIND_CONVERGENCE_MFD) {
//   RunNonlinearTest("mfd: default", "none");
// }
// TEST(OPERATOR_COUPLED_DIFFUSION_NONLINEAR_UPWIND_CONVERGENCE_MFD_JACOBIAN) {
//   RunNonlinearTest("mfd: default", "full");
// }


// -----------------------------------------------------------------------------
// debugging tests, simply run one run and write more info to files
// -----------------------------------------------------------------------------
// TEST(OPERATOR_COUPLED_DIFFUSION_NONLINEAR_MFD_NONE) {
//   RunNonlinearProblem("fv: default", "none", true, 2,2, true, 1.);
// }
// TEST(OPERATOR_COUPLED_DIFFUSION_NONLINEAR_MFD_FULL) {
//   RunNonlinearProblem("fv: default", "full", true, 2,2, true, 1.);
// }
// TEST(OPERATOR_COUPLED_DIFFUSION_WRITE_MATRIX) {
//   RunInverseProblem("mfd: default", true, 2,2, true);
// }
