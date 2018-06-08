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
#include <vector>

// TPLs
#include "Epetra_MultiVector.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "UnitTest++.h"

// Amanzi
#include "DG_Modal.hh"
#include "Explicit_TI_RK.hh"
#include "GMVMesh.hh"
#include "LinearOperatorPCG.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "MeshMapsFactory.hh"
#include "NumericalIntegration.hh"
#include "Tensor.hh"
#include "WhetStone_typedefs.hh"

// Amanzi::Operators
#include "OperatorDefs.hh"
#include "PDE_Abstract.hh"
#include "PDE_AdvectionRiemann.hh"
#include "PDE_Reaction.hh"

#include "AnalyticDG00.hh"
#include "AnalyticDG01.hh"
#include "AnalyticDG04.hh"

namespace Amanzi {

class RemapDG : public Explicit_TI::fnBase<CompositeVector> {
 public:
  RemapDG(const Teuchos::RCP<const AmanziMesh::Mesh> mesh,
          const std::shared_ptr<WhetStone::MeshMaps> maps) : mesh_(mesh), maps_(maps) {
    // mesh data
    ncells_owned_ = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
    ncells_wghost_ = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
    nfaces_owned_ = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
    nfaces_wghost_ = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);
    dim_ = mesh_->space_dimension();

    // face velocities
    velf_vec_.resize(nfaces_wghost_);
    for (int f = 0; f < nfaces_wghost_; ++f) {
      maps_->VelocityFace(f, velf_vec_[f]);
    }

    // cell-baced velocities and Jacobian matrices
    WhetStone::Entity_ID_List faces;
    uc_.resize(ncells_wghost_);
    J_.resize(ncells_wghost_);

    for (int c = 0; c < ncells_wghost_; ++c) {
      mesh_->cell_get_faces(c, &faces);

      std::vector<WhetStone::VectorPolynomial> vvf;
      for (int n = 0; n < faces.size(); ++n) {
        vvf.push_back(velf_vec_[faces[n]]);
      }

      maps->VelocityCell(c, vvf, uc_[c]);
      maps->Jacobian(uc_[c], J_[c]);
      // maps->JacobianCell(c, vvf, J_[c]);
    }

    velf_ = Teuchos::rcp(new std::vector<WhetStone::Polynomial>(nfaces_wghost_));
    velc_ = Teuchos::rcp(new std::vector<WhetStone::VectorPolynomial>(ncells_wghost_));
    jac_ = Teuchos::rcp(new std::vector<WhetStone::VectorPolynomial>(ncells_wghost_));

    tprint_ = 0.0;
    tl2_ = 0.0;
    nfun_ = 0;
    l2norm_ = -1.0;
  };
  ~RemapDG() {};

  void Setup(const Teuchos::Ptr<Operators::PDE_Abstract> op_adv,
             const Teuchos::Ptr<Operators::PDE_AdvectionRiemann> op_flux,
             const Teuchos::Ptr<Operators::PDE_Reaction> op_reac) {
    op_adv_ = op_adv;
    op_flux_ = op_flux;
    op_reac_ = op_reac;
  }

  virtual void Functional(double t, const CompositeVector& u, CompositeVector& f) override {
    UpdateGeometricQuantities_(t);

    // -- populate operators
    op_adv_->SetupPolyVector(velc_);
    op_adv_->UpdateMatrices();

    op_flux_->Setup(velc_, velf_);
    op_flux_->UpdateMatrices(velf_.ptr());

    op_reac_->Setup(jac_);
    op_reac_->UpdateMatrices(Teuchos::null);

    // -- solve the problem with mass matrix
    auto& matrices = op_reac_->local_matrices()->matrices;
    for (int n = 0; n < matrices.size(); ++n) {
      matrices[n].Inverse();
    }
    auto global_reac = op_reac_->global_operator();
    CompositeVector& x = *global_reac->rhs();  // reusing internal memory
    global_reac->Apply(u, x);

    // statistics
    nfun_++;
    Epetra_MultiVector& xc = *x.ViewComponent("cell");
    int nk = xc.NumVectors();
    double xmax[nk], xmin[nk];
    xc.MaxValue(xmax);
    xc.MinValue(xmin);

    if (fabs(tprint_ - t) < 1e-6 && mesh_->get_comm()->MyPID() == 0) {
      printf("t=%8.5f  L2=%9.5g  nfnc=%3d  umax: ", t, l2norm_, nfun_);
      for (int i = 0; i < nk; ++i) printf("%9.5g ", xmax[i]);
      printf("\n");
      // printf("                                    umin: ", t, l2norm_, nfun_);
      // for (int i = 0; i < nk; ++i) printf("%9.5g ", xmin[i]);
      // printf("\n");
      tprint_ += 0.1;
    } 

    // -- calculate right-hand_side
    op_flux_->global_operator()->Apply(x, f);
  }

  double L2Norm(double t, const CompositeVector& p1) {
    if (fabs(tl2_ - t) < 1e-6) {
      CompositeVector p2(p1);

      ChangeVariables(t, p1, p2, false);
      p1.Dot(p2, &l2norm_);
      tl2_ += 0.1;
    }
    return l2norm_;
  }

  // access 
  const std::vector<WhetStone::VectorPolynomial> velf_vec() const { return velf_vec_; }
  const std::vector<WhetStone::VectorPolynomial> uc() const { return uc_; }
  const std::vector<WhetStone::MatrixPolynomial> J() const { return J_; }
  const Teuchos::RCP<std::vector<WhetStone::VectorPolynomial> > jac() const { return jac_; }

  void ChangeVariables(double t, const CompositeVector& p1, CompositeVector& p2, bool flag) {
    UpdateGeometricQuantities_(t);
    op_reac_->Setup(jac_);
    op_reac_->UpdateMatrices(Teuchos::null);

    auto global_reac = op_reac_->global_operator();
    if (flag) {
      global_reac->Apply(p1, p2);
    } else {
      auto& matrices = op_reac_->local_matrices()->matrices;
      for (int n = 0; n < matrices.size(); ++n) {
        matrices[n].Inverse();
      }
      global_reac->Apply(p1, p2);
    }
  }

  // estimate trace error
  void CalculateTraceError(double* errl2, double* errinf) {
    *errl2 = 0.0;
    *errinf = 0.0;

    WhetStone::Entity_ID_List faces;
    std::vector<const WhetStone::Polynomial*> polys(2);
    WhetStone::NumericalIntegration numi(mesh_, false);

    for (int c = 0; c < ncells_owned_; ++c) {
      mesh_->cell_get_faces(c, &faces);
      int nfaces = faces.size();

      for (int n = 0; n < nfaces; ++n) {
        int f = faces[n];

        WhetStone::Polynomial tmp(uc_[c][0]);
        tmp -= velf_vec_[f][0];

        polys[0] = &tmp;
        polys[1] = &tmp;
        double err = numi.IntegratePolynomialsFace(f, polys);
        err /= mesh_->face_area(f);

        *errinf = std::max(*errinf, std::pow(err, 0.5));
        *errl2 += err * mesh_->cell_volume(c) / nfaces;
      }
    }

    double errtmp = *errinf;
    mesh_->get_comm()->MaxAll(&errtmp, errinf, 1);

    errtmp = *errl2;
    mesh_->get_comm()->SumAll(&errtmp, errl2, 1);
    *errl2 = std::pow(*errl2, 0.5);
  }

  // estimate node error
  void CalculateNodeError(double* errl2, double* errinf) {
    *errl2 = 0.0;
    *errinf = 0.0;

    WhetStone::Entity_ID_List faces, nodes;
    AmanziGeometry::Point xv(dim_);

    for (int c = 0; c < ncells_owned_; ++c) {
      mesh_->cell_get_faces(c, &faces);
      int nfaces = faces.size();

      for (int n = 0; n < nfaces; ++n) {
        int f = faces[n];
        mesh_->face_get_nodes(f, &nodes);
        mesh_->node_get_coordinates(nodes[0], &xv);

        WhetStone::Polynomial tmp(uc_[c][0]);
        tmp -= velf_vec_[f][0];

        double err = tmp.Value(xv);
        *errinf = std::max(*errinf, fabs(err));
        *errl2 += err * err * mesh_->cell_volume(c) / nfaces;
      }
    }

    double errtmp = *errinf;
    mesh_->get_comm()->MaxAll(&errtmp, errinf, 1);

    errtmp = *errl2;
    mesh_->get_comm()->SumAll(&errtmp, errl2, 1);
    *errl2 = std::pow(*errl2, 0.5);
  }

  // estimate determinant error
  void CalculateDeterminantError(Teuchos::RCP<const AmanziMesh::Mesh> mesh1,
                                 double* errl2, double* errinf) {
    *errl2 = 0.0;
    *errinf = 0.0;

    WhetStone::NumericalIntegration numi(mesh_, false);

    UpdateGeometricQuantities_(1.0);

    for (int c = 0; c < ncells_owned_; ++c) {
      double tmp = numi.IntegratePolynomialCell(c, (*jac_)[c][0]);
      double err = std::fabs(tmp - mesh1->cell_volume(c));
      err /= mesh_->cell_volume(c); 

      *errinf = std::max(*errinf, err);
      *errl2 += err * err * mesh_->cell_volume(c);
    }

    double errtmp = *errinf;
    mesh_->get_comm()->MaxAll(&errtmp, errinf, 1);

    errtmp = *errl2;
    mesh_->get_comm()->SumAll(&errtmp, errl2, 1);
    *errl2 = std::pow(*errl2, 0.5);
  }

 private:
  void UpdateGeometricQuantities_(double t) {
    WhetStone::VectorPolynomial cn;

    for (int f = 0; f < nfaces_wghost_; ++f) {
      // cn = j J^{-t} N dA
      maps_->NansonFormula(f, t, velf_vec_[f], cn);
      (*velf_)[f] = velf_vec_[f] * cn;
    }

    WhetStone::MatrixPolynomial C;
    for (int c = 0; c < ncells_wghost_; ++c) {
      maps_->Cofactors(t, J_[c], C);
    
      // cell-based pseudo velocity -C^t u 
      int nC = C.size();
      (*velc_)[c].resize(nC);

      int kC = nC / dim_;
      for (int n = 0; n < kC; ++n) {
        int m = n * dim_;
        for (int i = 0; i < dim_; ++i) {
          (*velc_)[c][m + i].Reshape(dim_, 0, true);
          (*velc_)[c][m + i].set_origin(uc_[c][0].origin());

          for (int k = 0; k < dim_; ++k) {
            (*velc_)[c][m + i] -= C[m + k][i] * uc_[c][m + k];
          }
        }
      }

      // determinant of Jacobian
      maps_->Determinant(t, J_[c], (*jac_)[c]);
    }
  }

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  int ncells_owned_, ncells_wghost_, nfaces_owned_, nfaces_wghost_;
  int dim_;

  std::shared_ptr<WhetStone::MeshMaps> maps_;
  std::vector<WhetStone::VectorPolynomial> uc_;
  std::vector<WhetStone::MatrixPolynomial> J_;
  Teuchos::RCP<std::vector<WhetStone::VectorPolynomial> > jac_;

  Teuchos::Ptr<Operators::PDE_Abstract> op_adv_;
  Teuchos::Ptr<Operators::PDE_AdvectionRiemann> op_flux_;
  Teuchos::Ptr<Operators::PDE_Reaction> op_reac_;

  Teuchos::RCP<std::vector<WhetStone::VectorPolynomial> > velc_;
  Teuchos::RCP<std::vector<WhetStone::Polynomial> > velf_;

  std::vector<WhetStone::VectorPolynomial> velf_vec_;

  double tprint_, tl2_;
  int nfun_;
  double l2norm_;
};

} // namespace Amanzi


/* *****************************************************************
* Remap of polynomilas in two dimensions. Explicit time scheme.
* Dual formulation places gradient and jumps on a test function.
***************************************************************** */
void RemapTestsDualRK(int order_p, int order_u,
                      const Amanzi::Explicit_TI::method_t& rk_method,
                      std::string maps_name, std::string file_name,
                      int nx, int ny, int nz, double dt) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  int dim = (nz == 0) ? 2 : 3;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();
  if (MyPID == 0) std::cout << "\nTest: " << dim << "D remap, dual formulation:"
                            << " mesh=" << ((ny == 0) ? file_name : "square")
                            << ", orders: " << order_p << " " << order_u 
                            << ", maps=" << maps_name << std::endl;

  // polynomial space
  WhetStone::Polynomial pp(dim, order_p);
  int nk = pp.size();

  // create initial mesh
  MeshFactory meshfactory(&comm);
  meshfactory.set_partitioner(AmanziMesh::Partitioner_type::ZOLTAN_RCB);
  meshfactory.preference(FrameworkPreference({MSTK}));

  Teuchos::RCP<const Mesh> mesh0;
  if (dim == 2) {
    if (ny != 0) 
      mesh0 = meshfactory(0.0, 0.0, 1.0, 1.0, nx, ny);
    else 
      mesh0 = meshfactory(file_name, Teuchos::null);
  } else {
    mesh0 = meshfactory(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, nx, ny, nz, Teuchos::null, true, true);
  }

  int ncells_owned = mesh0->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int ncells_wghost = mesh0->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
  int nnodes_wghost = mesh0->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::ALL);

  // create second and auxiliary mesh
  Teuchos::RCP<Mesh> mesh1;
  if (dim == 2) {
    if (ny != 0) 
      mesh1 = meshfactory(0.0, 0.0, 1.0, 1.0, nx, ny);
    else 
      mesh1 = meshfactory(file_name, Teuchos::null);
  } else {
    mesh1 = meshfactory(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, nx, ny, nz, Teuchos::null, true, true);
  }

  // deform the second mesh
  AmanziGeometry::Point xv(dim), yv(dim), xref(dim), uv(dim);
  Entity_ID_List nodeids;
  AmanziGeometry::Point_List new_positions, final_positions;

  for (int v = 0; v < nnodes_wghost; ++v) {
    mesh1->node_get_coordinates(v, &xv);
    yv = xv;

    double ds(0.0001);
    for (int i = 0; i < 10000; ++i) {
      if (dim == 2) {
        uv[0] = 0.2 * std::sin(M_PI * xv[0]) * std::cos(M_PI * xv[1]);
        uv[1] =-0.2 * std::cos(M_PI * xv[0]) * std::sin(M_PI * xv[1]);
      } else {
        uv[0] = 0.2 * std::sin(M_PI * xv[0]) * std::cos(M_PI * xv[1]) * std::cos(M_PI * xv[2]);
        uv[1] =-0.1 * std::cos(M_PI * xv[0]) * std::sin(M_PI * xv[1]) * std::cos(M_PI * xv[2]);
        uv[2] =-0.1 * std::cos(M_PI * xv[0]) * std::cos(M_PI * xv[1]) * std::sin(M_PI * xv[2]);
        // uv[0] = 0.2 * std::sin(M_PI * xv[0]) * std::cos(M_PI * xv[1]);
        // uv[1] =-0.2 * std::cos(M_PI * xv[0]) * std::sin(M_PI * xv[1]);
      }
      xv += uv * ds;
    }
    /*
    xv[0] = yv[0] * yv[1] + (1.0 - yv[1]) * std::pow(yv[0], 0.8);
    xv[1] = yv[1] * yv[0] + (1.0 - yv[0]) * std::pow(yv[1], 0.8);

    for (int i = 0; i < 2; ++i) {
      xv[i] = yv[i];
      if (xv[i] > 0.01 && xv[i] < 0.99) 
        xv[i] += sin(3 * yv[1 - i]) * sin(4 * yv[i]) / nx / 8;
    }
    */

    nodeids.push_back(v);
    new_positions.push_back(xv);
  }
  mesh1->deform(nodeids, new_positions, false, &final_positions);

  // little factory of mesh maps
  Teuchos::ParameterList map_list;
  map_list.set<std::string>("method", "CrouzeixRaviart")
          .set<int>("method order", order_u)
          .set<std::string>("projector", "H1 harmonic")
          .set<std::string>("map name", maps_name);
  
  WhetStone::MeshMapsFactory maps_factory;
  auto maps = maps_factory.Create(map_list, mesh0, mesh1);

  // numerical integration
  WhetStone::NumericalIntegration numi(mesh0, false);

  // basic remap algorithm
  RemapDG remap(mesh0, maps);

  double err2, err8;
  remap.CalculateNodeError(&err2, &err8);
  if (MyPID == 0) std::cout << "Nodes: " << err2 << " " << err8 << std::endl;
  remap.CalculateTraceError(&err2, &err8);
  if (MyPID == 0) std::cout << "Trace: " << err2 << " " << err8 << std::endl;
  remap.CalculateDeterminantError(mesh1, &err2, &err8);
  if (MyPID == 0) std::cout << "Det:   " << err2 << " " << err8 << std::endl;

  // create and initialize cell-based field 
  CompositeVectorSpace cvs1;
  cvs1.SetMesh(mesh0)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, nk);
  Teuchos::RCP<CompositeVector> p1 = Teuchos::rcp(new CompositeVector(cvs1));
  Epetra_MultiVector& p1c = *p1->ViewComponent("cell", true);

  // we need dg to use correct scaling of basis functions
  WhetStone::DG_Modal dg(order_p, mesh0, "orthonormalized");
  AnalyticDG04 ana(mesh0, order_p);
  ana.InitialGuess(dg, p1c, 1.0);

  // initial mass
  double mass0(0.0);
  for (int c = 0; c < ncells_owned; c++) {
    WhetStone::DenseVector data(nk);
    for (int i = 0; i < nk; ++i) {
      data(i) = p1c[i][c];
    }
    auto poly = dg.cell_basis(c).CalculatePolynomial(mesh0, c, order_p, data);

    numi.ChangeBasisRegularizedToNatural(c, poly);
    mass0 += numi.IntegratePolynomialCell(c, poly);
  }
  double mass_tmp(mass0);
  mesh0->get_comm()->SumAll(&mass_tmp, &mass0, 1);

  // allocate memory for global variables
  CompositeVectorSpace cvs2;
  cvs2.SetMesh(mesh1)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, nk);
  CompositeVector p2(cvs2);
  Epetra_MultiVector& p2c = *p2.ViewComponent("cell");

  // create flux operator
  Teuchos::ParameterList plist;
  plist.set<std::string>("method", "dg modal")
       .set<std::string>("dg basis", "orthonormalized")
       .set<std::string>("matrix type", "flux")
       .set<std::string>("flux formula", "downwind")
       .set<int>("method order", order_p)
       .set<bool>("jump operator on test function", true);

  plist.sublist("schema domain")
      .set<std::string>("base", "face")
      .set<Teuchos::Array<std::string> >("location", std::vector<std::string>({"cell"}))
      .set<Teuchos::Array<std::string> >("type", std::vector<std::string>({"scalar"}))
      .set<Teuchos::Array<int> >("number", std::vector<int>({nk}));

  plist.sublist("schema range") = plist.sublist("schema domain");

  Teuchos::RCP<PDE_AdvectionRiemann> op_flux = Teuchos::rcp(new PDE_AdvectionRiemann(plist, mesh0));
  auto global_op = op_flux->global_operator();

  // Attach volumetric advection operator to the flux operator.
  // We modify the existing parameter list.
  plist.set<std::string>("matrix type", "advection")
       .set<bool>("gradient operator on test function", true);
  plist.sublist("schema domain").set<std::string>("base", "cell");
  plist.sublist("schema range") = plist.sublist("schema domain");

  Teuchos::RCP<PDE_Abstract> op_adv = Teuchos::rcp(new PDE_Abstract(plist, global_op));

  // approximate accumulation term using the reaction operator
  plist.sublist("schema")
      .set<std::string>("base", "cell")
      .set<Teuchos::Array<std::string> >("location", std::vector<std::string>({"cell"}))
      .set<Teuchos::Array<std::string> >("type", std::vector<std::string>({"scalar"}))
      .set<Teuchos::Array<int> >("number", std::vector<int>({nk}));

  Teuchos::RCP<PDE_Reaction> op_reac = Teuchos::rcp(new PDE_Reaction(plist, mesh0));

  // explicit time integration
  CompositeVector p1aux(*p1);
  remap.Setup(op_adv.ptr(), op_flux.ptr(), op_reac.ptr());
  Explicit_TI::RK<CompositeVector> rk(remap, rk_method, p1aux);

  remap.ChangeVariables(0.0, *p1, p1aux, true);

  int nstep(0), nstep_dbg(0);
  double gcl, gcl_err(0.0);
  double t(0.0), tend(1.0);
  while(t < tend - dt/2) {
    remap.L2Norm(t, p1aux);
    rk.TimeStep(t, dt, p1aux, *p1);

    *p1aux.ViewComponent("cell") = *p1->ViewComponent("cell");

    t += dt;
    nstep++;
  }

  remap.ChangeVariables(1.0, *p1, p2, false);

  // calculate error in the new basis
  double l20_err(0.0), inf0_err(0.0);
  double pl2_err(0.0), pinf_err(0.0), area(0.0);
  double ql2_err(0.0), qinf_err(0.0), mass1(0.0);

  Entity_ID_List nodes;
  std::vector<int> dirs;
  AmanziGeometry::Point v0(dim), v1(dim), tau(dim);

  CompositeVectorSpace cvs3;
  cvs3.SetMesh(mesh1)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, 1);

  CompositeVector p2err(cvs3);
  Epetra_MultiVector& p2c_err = *p2err.ViewComponent("cell");
  p2err.PutScalar(0.0);

  CompositeVector q2(p2);
  Epetra_MultiVector& q2c = *q2.ViewComponent("cell");
  q2c = p2c;

  for (int c = 0; c < ncells_owned; ++c) {
    const AmanziGeometry::Point& xc0 = mesh0->cell_centroid(c);
    double area_c(mesh1->cell_volume(c));

    WhetStone::DenseVector data(nk);
    for (int i = 0; i < nk; ++i) {
      data(i) = p2c[i][c];
    }
    auto poly = dg.cell_basis(c).CalculatePolynomial(mesh0, c, order_p, data);
    numi.ChangeBasisRegularizedToNatural(c, poly);

    // const AmanziGeometry::Point& xg = cell_geometric_center(*mesh1, c);
    const AmanziGeometry::Point& xg = mesh1->cell_centroid(c);
    double err = poly.Value(xc0) - ana.SolutionExact(xg, 1.0);

    inf0_err = std::max(inf0_err, fabs(err));
    l20_err += err * err * area_c;
    p2c_err[0][c] = fabs(err);

    if (nk > 1) {
      mesh0->cell_get_nodes(c, &nodes);
      int nnodes = nodes.size();  
      for (int i = 0; i < nnodes; ++i) {
        mesh0->node_get_coordinates(nodes[i], &v0);
        mesh1->node_get_coordinates(nodes[i], &v1);

        double tmp = poly.Value(v0);
        tmp -= ana.SolutionExact(v1, 1.0);
        pinf_err = std::max(pinf_err, fabs(tmp));
        pl2_err += tmp * tmp * area_c / nnodes;

        p2c_err[0][c] = std::max(p2c_err[0][c], fabs(tmp));
      }
    }

    area += area_c;

    double mass1_c;
    auto& jac = remap.jac();
    int quad_order = (*jac)[c][0].order() + poly.order();

    if (maps_name == "PEM") {
      AmanziMesh::Entity_ID_List faces, nodes;
      mesh0->cell_get_faces(c, &faces);
      int nfaces = faces.size();

      std::vector<AmanziGeometry::Point> xy(3);
      // xy[0] = cell_geometric_center(*mesh0, c);
      xy[0] = mesh0->cell_centroid(c);

      mass1_c = 0.0;
      for (int n = 0; n < nfaces; ++n) {
        int f = faces[n];
        mesh0->face_get_nodes(f, &nodes);
        mesh0->node_get_coordinates(nodes[0], &(xy[1]));
        mesh0->node_get_coordinates(nodes[1], &(xy[2]));

        std::vector<const WhetStone::WhetStoneFunction*> polys(2);
        polys[0] = &(*jac)[c][n];
        polys[1] = &poly;
        mass1_c += numi.IntegrateFunctionsTriangle(xy, polys, quad_order);
      }
    } else {
      WhetStone::Polynomial tmp((*jac)[c][0]);
      tmp.ChangeOrigin(mesh0->cell_centroid(c));
      poly *= tmp;

      mass1_c = numi.IntegratePolynomialCell(c, poly);
    }
    mass1 += mass1_c;

    // optional projection on the space of polynomials 
    if (order_p > 0 && order_p < 3 && dim == 2) {
      poly = dg.cell_basis(c).CalculatePolynomial(mesh0, c, order_p, data);
      numi.ChangeBasisRegularizedToNatural(c, poly);

      maps->ProjectPolynomial(c, poly);
      poly.ChangeOrigin(mesh1->cell_centroid(c));

      if (order_p == 1) {
        poly(0, 0) = mass1_c / mesh1->cell_volume(c);
      }

      // error in the projected solution
      mesh0->cell_get_nodes(c, &nodes);
      int nnodes = nodes.size();  
      for (int i = 0; i < nnodes; ++i) {
        mesh1->node_get_coordinates(nodes[i], &v1);

        double tmp = poly.Value(v1);
        tmp -= ana.SolutionExact(v1, 1.0);
        qinf_err = std::max(qinf_err, fabs(tmp));
        ql2_err += tmp * tmp * area_c / nnodes;
      }

      // save the projected function
      poly.GetPolynomialCoefficients(data);
      for (int i = 0; i < nk; ++i) {
        q2c[i][c] = data(i);
      }
    }
  }

  // parallel collective operations
  double err_in[5] = {pl2_err, area, mass1, ql2_err, l20_err};
  double err_out[5];
  mesh1->get_comm()->SumAll(err_in, err_out, 5);

  double err_tmp = pinf_err;
  mesh1->get_comm()->MaxAll(&err_tmp, &pinf_err, 1);

  err_tmp = qinf_err;
  mesh1->get_comm()->MaxAll(&err_tmp, &qinf_err, 1);

  err_tmp = inf0_err;
  mesh1->get_comm()->MaxAll(&err_tmp, &inf0_err, 1);

  err_tmp = gcl_err;
  mesh1->get_comm()->MaxAll(&err_tmp, &gcl_err, 1);

  // error tests
  pl2_err = std::pow(err_out[0], 0.5);
  ql2_err = std::pow(err_out[3], 0.5);
  l20_err = std::pow(err_out[4], 0.5);
  CHECK(pl2_err < 0.12 / (order_p + 1));

  if (MyPID == 0) {
    printf("nx=%3d  L2=%12.8g %12.8g %12.8g  Inf=%12.8g %12.8g %12.8g dMass=%10.4g  dArea=%10.6g\n", 
        nx, l20_err, pl2_err, ql2_err, 
            inf0_err, pinf_err, qinf_err, err_out[2] - mass0, 1.0 - err_out[1]);
  }

  // visualization
  if (MyPID == 0) {
    GMV::open_data_file(*mesh1, (std::string)"operators.gmv");
    GMV::start_data();
    GMV::write_cell_data(p2c, 0, "remaped");
    GMV::write_cell_data(q2c, 0, "remaped-prj");
    if (order_p > 0) {
      GMV::write_cell_data(q2c, 1, "gradx-prj");
      GMV::write_cell_data(q2c, 2, "grady-prj");
    }

    GMV::write_cell_data(p2c_err, 0, "error");
    GMV::close_data_file();
  }
}

/*
TEST(REMAP_DUAL_FEM) {
  RemapTestsDualRK(0,1, Amanzi::Explicit_TI::heun_euler, "dg modal", "FEM", "", 10,10,0, 0.1);
  RemapTestsDualRK(1,2, Amanzi::Explicit_TI::heun_euler, "dg modal", "FEM", "", 10,10,0, 0.1);
}
*/

TEST(REMAP_DUAL_VEM) {
  RemapTestsDualRK(0,1, Amanzi::Explicit_TI::heun_euler, "VEM", "test/median15x16.exo", 0,0,0, 0.05);
  RemapTestsDualRK(1,2, Amanzi::Explicit_TI::heun_euler, "VEM", "test/median15x16.exo", 0,0,0, 0.05);
  // RemapTestsDualRK(2,3, Amanzi::Explicit_TI::heun_euler, "VEM", "test/median15x16.exo", 0,0,0, 0.05);
  // RemapTestsDualRK(0,1, Amanzi::Explicit_TI::heun_euler, "VEM", "", 5,5,5, 0.2);
  // RemapTestsDualRK(1,2, Amanzi::Explicit_TI::heun_euler, "VEM", "", 5,5,5, 0.1);
}

/*
TEST(REMAP2D_DG_QUADRATURE_ERROR) {
  RemapTestsDualRK(1,1, Amanzi::Explicit_TI::tvd_3rd_order, "VEM", "",  16, 16,0, 0.05);
  RemapTestsDualRK(1,1, Amanzi::Explicit_TI::tvd_3rd_order, "VEM", "",  32, 32,0, 0.05 / 2);
  RemapTestsDualRK(1,1, Amanzi::Explicit_TI::tvd_3rd_order, "VEM", "",  64, 64,0, 0.05 / 4);
  RemapTestsDualRK(1,1, Amanzi::Explicit_TI::tvd_3rd_order, "VEM", "", 128,128,0, 0.05 / 8);
  RemapTestsDualRK(1,1, Amanzi::Explicit_TI::tvd_3rd_order, "VEM", "", 256,256,0, 0.05 / 16);
} 

TEST(REMAP2D_DG_QUADRATURE_ERROR) {
  RemapTestsDualRK(1,1, Amanzi::Explicit_TI::tvd_3rd_order, "PEM", "test/median15x16.exo", 16,0,0, 0.05);
  RemapTestsDualRK(1,1, Amanzi::Explicit_TI::tvd_3rd_order, "PEM", "test/median32x33.exo", 32,0,0, 0.05 / 2);
  RemapTestsDualRK(1,1, Amanzi::Explicit_TI::tvd_3rd_order, "PEM", "test/median63x64.exo", 64,0,0, 0.05 / 4);
  RemapTestsDualRK(1,1, Amanzi::Explicit_TI::tvd_3rd_order, "PEM", "test/median127x128.exo", 128,0,0, 0.05 / 8);
  RemapTestsDualRK(1,1, Amanzi::Explicit_TI::tvd_3rd_order, "PEM", "test/median255x256.exo", 256,0,0, 0.05 / 16);
}

TEST(REMAP2D_DG_QUADRATURE_ERROR) {
  RemapTestsDualRK(1,1, Amanzi::Explicit_TI::tvd_3rd_order, "PEM", "test/mesh_poly20x20.exo", 20,0,0, 0.05);
  RemapTestsDualRK(1,1, Amanzi::Explicit_TI::tvd_3rd_order, "PEM", "test/mesh_poly40x40.exo", 40,0,0, 0.05 / 2);
  RemapTestsDualRK(1,1, Amanzi::Explicit_TI::tvd_3rd_order, "PEM", "test/mesh_poly80x80.exo", 80,0,0, 0.05 / 4);
  RemapTestsDualRK(1,1, Amanzi::Explicit_TI::tvd_3rd_order, "PEM", "test/mesh_poly160x160.exo", 160,0,0, 0.05 / 8);
}

TEST(REMAP2D_DG_QUADRATURE_ERROR) {
  RemapTestsDualRK(2,1, Amanzi::Explicit_TI::tvd_3rd_order, "PEM", "test/random10.exo", 10,0,0, 0.05);
  RemapTestsDualRK(2,1, Amanzi::Explicit_TI::tvd_3rd_order, "PEM", "test/random20.exo", 20,0,0, 0.05 / 2);
  RemapTestsDualRK(2,1, Amanzi::Explicit_TI::tvd_3rd_order, "PEM", "test/random40.exo", 40,0,0, 0.05 / 4);
}

TEST(REMAP2D_DG_QUADRATURE_ERROR) {
  RemapTestsDualRK(2,1, Amanzi::Explicit_TI::tvd_3rd_order, "PEM", "test/triangular8.exo", 0,0,0, 0.025);
  RemapTestsDualRK(2,1, Amanzi::Explicit_TI::tvd_3rd_order, "PEM", "test/triangular16.exo", 0,0,0, 0.025 / 2);
  RemapTestsDualRK(2,1, Amanzi::Explicit_TI::tvd_3rd_order, "PEM", "test/triangular32.exo", 0,0,0, 0.025 / 4);
  RemapTestsDualRK(2,1, Amanzi::Explicit_TI::tvd_3rd_order, "PEM", "test/triangular64.exo", 0,0,0, 0.025 / 8);
  RemapTestsDualRK(2,1, Amanzi::Explicit_TI::tvd_3rd_order, "PEM", "test/triangular128.exo", 0,0,0, 0.025 / 16);
}
*/
