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
#include "LinearOperatorPCG.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "MeshMapsFactory.hh"
#include "NumericalIntegration.hh"
#include "OutputXDMF.hh"
#include "Tensor.hh"
#include "WhetStoneDefs.hh"

// Amanzi::Operators
#include "OperatorDefs.hh"
#include "PDE_Abstract.hh"
#include "PDE_AdvectionRiemann.hh"
#include "PDE_Reaction.hh"

#include "AnalyticDG00.hh"
#include "AnalyticDG04.hh"

namespace Amanzi {

class RemapDG : public Explicit_TI::fnBase<CompositeVector> {
 public:
  RemapDG(const Teuchos::RCP<const AmanziMesh::Mesh> mesh0,
          const Teuchos::RCP<const AmanziMesh::Mesh> mesh1,
          Teuchos::ParameterList& plist,
          bool high_order_velocity);
  ~RemapDG() {};

  virtual void FunctionalTimeDerivative(double t, const CompositeVector& u, CompositeVector& f) override;

  void ChangeVariables(double t, const CompositeVector& p1, CompositeVector& p2, bool flag);

  double L2Norm(double t, const CompositeVector& p1);

  // access 
  const std::vector<WhetStone::VectorPolynomial> velf_vec() const { return velf_vec_; }
  const std::vector<WhetStone::VectorPolynomial> uc() const { return uc_; }
  const std::vector<WhetStone::MatrixPolynomial> J() const { return J_; }
  const std::vector<WhetStone::VectorPolynomial> jac() const { return *jac_; }
  const std::shared_ptr<WhetStone::MeshMaps> maps() const { return maps_; }

 private:
  void UpdateGeometricQuantities_(double t);

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh0_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh1_;
  int ncells_owned_, ncells_wghost_, nfaces_owned_, nfaces_wghost_;
  int dim_;

  const Teuchos::ParameterList plist_;
  std::shared_ptr<WhetStone::MeshMaps> maps_;

  // operators
  Teuchos::RCP<Operators::PDE_Abstract> op_adv_;
  Teuchos::RCP<Operators::PDE_AdvectionRiemann> op_flux_;
  Teuchos::RCP<Operators::PDE_Reaction>op_reac_;

  // geometric data
  std::vector<WhetStone::VectorPolynomial> uc_;
  std::vector<WhetStone::MatrixPolynomial> J_;
  Teuchos::RCP<std::vector<WhetStone::VectorPolynomial> > jac_;

  Teuchos::RCP<std::vector<WhetStone::VectorPolynomial> > velc_;
  Teuchos::RCP<std::vector<WhetStone::Polynomial> > velf_;
  std::vector<WhetStone::VectorPolynomial> velf_vec_;

  // statistics
  int nfun_;
  double tprint_, tl2_, l2norm_;
};


/* *****************************************************************
* Constructor
***************************************************************** */
RemapDG::RemapDG(const Teuchos::RCP<const AmanziMesh::Mesh> mesh0,
                 const Teuchos::RCP<const AmanziMesh::Mesh> mesh1,
                 Teuchos::ParameterList& plist,
                 bool high_order_velocity)
  : mesh0_(mesh0), mesh1_(mesh1), plist_(plist)
{
  // mesh data
  ncells_owned_ = mesh0_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  ncells_wghost_ = mesh0_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
  nfaces_owned_ = mesh0_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  nfaces_wghost_ = mesh0_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);
  dim_ = mesh0_->space_dimension();

  int order = plist.sublist("PK operator")
                   .sublist("flux operator").get<int>("method order");

  // create right-hand side operator
  // -- flux
  auto oplist = plist_.sublist("PK operator").sublist("flux operator");
  op_flux_ = Teuchos::rcp(new Operators::PDE_AdvectionRiemann(oplist, mesh0));
  auto global_op = op_flux_->global_operator();

  // -- advection
  oplist = plist_.sublist("PK operator").sublist("advection operator");
  op_adv_ = Teuchos::rcp(new Operators::PDE_Abstract(oplist, global_op));

  // create left-hand side operator
  oplist = plist_.sublist("PK operator").sublist("reaction operator");
  op_reac_ = Teuchos::rcp(new Operators::PDE_Reaction(oplist, mesh0));

  // factory of mesh maps
  auto map_list = plist_.sublist("maps");
  WhetStone::MeshMapsFactory maps_factory;
  maps_ = maps_factory.Create(map_list, mesh0, mesh1);

  // face velocities
  velf_vec_.resize(nfaces_wghost_);
  if (high_order_velocity) {
    AnalyticDG04 ana(mesh0_, order, true);

    for (int f = 0; f < nfaces_wghost_; ++f) {
      const auto& xf = mesh0_->face_centroid(f);
      ana.VelocityTaylor(xf, 0.0, velf_vec_[f]); 
      velf_vec_[f].ChangeOrigin(AmanziGeometry::Point(dim_));
    }
  } else {
    for (int f = 0; f < nfaces_wghost_; ++f) {
      maps_->VelocityFace(f, velf_vec_[f]);
    }
  } 

  // cell-baced velocities and Jacobian matrices
  WhetStone::Entity_ID_List faces;
  uc_.resize(ncells_wghost_);
  J_.resize(ncells_wghost_);

  for (int c = 0; c < ncells_wghost_; ++c) {
    mesh0_->cell_get_faces(c, &faces);

    std::vector<WhetStone::VectorPolynomial> vvf;
    for (int n = 0; n < faces.size(); ++n) {
      vvf.push_back(velf_vec_[faces[n]]);
    }

    maps_->VelocityCell(c, vvf, uc_[c]);
    maps_->Jacobian(uc_[c], J_[c]);
    // maps_->JacobianCell(c, vvf, J_[c]);
  }

  velf_ = Teuchos::rcp(new std::vector<WhetStone::Polynomial>(nfaces_wghost_));
  velc_ = Teuchos::rcp(new std::vector<WhetStone::VectorPolynomial>(ncells_wghost_));
  jac_ = Teuchos::rcp(new std::vector<WhetStone::VectorPolynomial>(ncells_wghost_));

  tprint_ = 0.0;
  tl2_ = 0.0;
  nfun_ = 0;
  l2norm_ = -1.0;
}


/* *****************************************************************
* Main routine: evaluation of functional
***************************************************************** */
void RemapDG::FunctionalTimeDerivative(
    double t, const CompositeVector& u, CompositeVector& f) {
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

  if (fabs(tprint_ - t) < 1e-6 && mesh0_->Comm()->getRank() == 0) {
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


/* *****************************************************************
* L2 norm
***************************************************************** */
double RemapDG::L2Norm(double t, const CompositeVector& p1) {
  if (fabs(tl2_ - t) < 1e-6) {
    CompositeVector p2(p1);

    ChangeVariables(t, p1, p2, false);
    p1.Dot(p2, &l2norm_);
    tl2_ += 0.1;
  }
  return l2norm_;
} 


/* *****************************************************************
* TBW
***************************************************************** */
void RemapDG::ChangeVariables(
    double t, const CompositeVector& p1, CompositeVector& p2, bool flag)
{
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


/* *****************************************************************
* Calculates various geometric quantaties of intermediate meshes.
***************************************************************** */
void RemapDG::UpdateGeometricQuantities_(double t)
{
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

    // experimental code:
    /*
    WhetStone::Entity_ID_List faces;
    mesh0_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    std::vector<WhetStone::VectorPolynomial> vvf(nfaces);
    for (int n = 0; n < nfaces; ++n) {
      auto tmp = velf_vec_[faces[n]];
      tmp.ChangeOrigin(C[0][0].origin());
      vvf[n].Multiply(C, tmp, true);
    }

    maps_->VelocityCell(c, vvf, (*velc_)[c]);
    (*velc_)[c] *= -1.0;
    */

    // determinant of Jacobian
    maps_->Determinant(t, J_[c], (*jac_)[c]);
  }
}

} // namespace Amanzi


/* *****************************************************************
* Remap of polynomilas in two dimensions. Explicit time scheme.
* Dual formulation places gradient and jumps on a test function.
***************************************************************** */
void RemapTestsDualRK(const Amanzi::Explicit_TI::method_t& rk_method,
                      std::string map_name, std::string file_name,
                      int nx, int ny, int nz, double dt,
                      int deform = 1) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  int dim = (nz == 0) ? 2 : 3;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm->getRank();

  // read parameter list
  std::string xmlFileName = "test/operator_remap.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  Teuchos::ParameterList plist = xmlreader.getParameters();

  int order = plist.sublist("PK operator")
                   .sublist("flux operator").get<int>("method order");

  int nk = WhetStone::PolynomialSpaceDimension(dim, order);

  // make modifications to the parameter list
  plist.sublist("maps").set<std::string>("map name", map_name);

  // print simulation header
  if (MyPID == 0) {
    const auto& map_list = plist.sublist("maps");
    int vel_order = map_list.get<int>("method order");
    std::string vel_method = map_list.get<std::string>("method");
    std::string vel_projector = map_list.get<std::string>("projector");
    std::string map_name = map_list.get<std::string>("map name");
      
    std::cout << "\nTest: " << dim << "D remap, dual formulation:"
              << " mesh=" << ((ny == 0) ? file_name : "square")
              << " deform=" << deform << std::endl;

    std::cout << "      discretization: order=" << order 
              << ", map=" << map_name << std::endl;

    std::cout << "      map details: order=" << vel_order 
              << ", projector=" << vel_projector 
              << ", method=\"" << vel_method << "\"" << std::endl;
  }

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

  // deform the second mesh using the specified algorithm
  AmanziGeometry::Point xv(dim), yv(dim), xref(dim), uv(dim);
  Entity_ID_List nodeids;
  AmanziGeometry::Point_List new_positions, final_positions;

  for (int v = 0; v < nnodes_wghost; ++v) {
    mesh1->node_get_coordinates(v, &xv);
    yv = xv;

    if (deform == 1) {
      double ds(0.0001);
      for (int i = 0; i < 10000; ++i) {
        if (dim == 2) {
          uv[0] = 0.2 * std::sin(M_PI * xv[0]) * std::cos(M_PI * xv[1]);
          uv[1] =-0.2 * std::cos(M_PI * xv[0]) * std::sin(M_PI * xv[1]);
        } else {
          uv[0] = 0.2 * std::sin(M_PI * xv[0]) * std::cos(M_PI * xv[1]) * std::cos(M_PI * xv[2]);
          uv[1] =-0.1 * std::cos(M_PI * xv[0]) * std::sin(M_PI * xv[1]) * std::cos(M_PI * xv[2]);
          uv[2] =-0.1 * std::cos(M_PI * xv[0]) * std::cos(M_PI * xv[1]) * std::sin(M_PI * xv[2]);
        }
        xv += uv * ds;
      }
    }
    else if (deform == 2) {
      xv[0] = yv[0] * yv[1] + (1.0 - yv[1]) * std::pow(yv[0], 0.8);
      xv[1] = yv[1] * yv[0] + (1.0 - yv[0]) * std::pow(yv[1], 0.8);
    }
    else if (deform == 3) {
      for (int i = 0; i < 2; ++i) {
        xv[i] = yv[i];
        if (xv[i] > 0.01 && xv[i] < 0.99) 
          xv[i] += sin(3 * yv[1 - i]) * sin(4 * yv[i]) / nx / 8;
      }
    }
    else if (deform == 4) {
      xv[0] = yv[0] + yv[0] * yv[1] * (1.0 - yv[0]) / 2;
      xv[1] = yv[1] + yv[0] * yv[1] * (1.0 - yv[1]) / 2;
    }
    else if (deform == 5) {
      xv[0] = yv[0] + yv[0] * (1.0 - yv[0]) / 2;
      xv[1] = yv[1] + yv[1] * (1.0 - yv[1]) / 2;
    }

    nodeids.push_back(v);
    new_positions.push_back(xv);
  }
  mesh1->deform(nodeids, new_positions, false, &final_positions);

  // basic remap algorithm
  bool high_order_velocity(false);
  RemapDG remap(mesh0, mesh1, plist, high_order_velocity);

  // create and initialize cell-based field 
  CompositeVectorSpace cvs1;
  cvs1.SetMesh(mesh0)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, nk);
  Teuchos::RCP<CompositeVector> p1 = Teuchos::rcp(new CompositeVector(cvs1));
  Epetra_MultiVector& p1c = *p1->ViewComponent("cell", true);

  // we need dg to use correct scaling of basis functions
  std::string basis = plist.sublist("PK operator")
                           .sublist("flux operator").get<std::string>("dg basis");
  WhetStone::DG_Modal dg(order, mesh0, basis);
  AnalyticDG04 ana(mesh0, order, true);
  ana.InitialGuess(dg, p1c, 1.0);

  // initial mass
  double mass0(0.0);
  WhetStone::NumericalIntegration numi(mesh0);

  for (int c = 0; c < ncells_owned; c++) {
    WhetStone::DenseVector data(nk);
    for (int i = 0; i < nk; ++i) {
      data(i) = p1c[i][c];
    }
    auto poly = dg.cell_basis(c).CalculatePolynomial(mesh0, c, order, data);
    mass0 += numi.IntegratePolynomialCell(c, poly);
  }
  double mass_tmp(mass0);
  mesh0->Comm()->SumAll(&mass_tmp, &mass0, 1);

  // allocate memory for global variables
  CompositeVectorSpace cvs2;
  cvs2.SetMesh(mesh1)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, nk);
  CompositeVector p2(cvs2);
  Epetra_MultiVector& p2c = *p2.ViewComponent("cell");

  // explicit time integration
  CompositeVector p1aux(*p1);
  Explicit_TI::RK<CompositeVector> rk(remap, rk_method, p1aux);

  remap.ChangeVariables(0.0, *p1, p1aux, true);

  int nstep(0), nstep_dbg(0);
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
    auto poly = dg.cell_basis(c).CalculatePolynomial(mesh0, c, order, data);

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
    int quad_order = jac[c][0].order() + poly.order();

    if (map_name == "PEM") {
      AmanziMesh::Entity_ID_List faces, nodes;
      mesh0->cell_get_faces(c, &faces);
      int nfaces = faces.size();

      std::vector<AmanziGeometry::Point> xy(3);
      xy[0] = mesh0->cell_centroid(c);

      mass1_c = 0.0;
      for (int n = 0; n < nfaces; ++n) {
        int f = faces[n];
        mesh0->face_get_nodes(f, &nodes);
        mesh0->node_get_coordinates(nodes[0], &(xy[1]));
        mesh0->node_get_coordinates(nodes[1], &(xy[2]));

        std::vector<const WhetStone::WhetStoneFunction*> polys(2);
        polys[0] = &jac[c][n];
        polys[1] = &poly;
        mass1_c += numi.IntegrateFunctionsSimplex(xy, polys, quad_order);
      }
    } else {
      WhetStone::Polynomial tmp(jac[c][0]);
      tmp.ChangeOrigin(mesh0->cell_centroid(c));
      poly *= tmp;

      mass1_c = numi.IntegratePolynomialCell(c, poly);
    }
    mass1 += mass1_c;

    // optional projection on the space of polynomials 
    if (order > 0 && order < 3 && dim == 2) {
      poly = dg.cell_basis(c).CalculatePolynomial(mesh0, c, order, data);
      remap.maps()->ProjectPolynomial(c, poly);
      poly.ChangeOrigin(mesh1->cell_centroid(c));

      if (order == 1) {
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
      for (int i = 0; i < nk; ++i) {
        q2c[i][c] = poly(i);
      }
    }
  }

  // error in GCL
  double gcl_err(0.0);
  auto& jac = remap.jac();

  for (int c = 0; c < ncells_owned; ++c) {
    double vol1 = numi.IntegratePolynomialCell(c, jac[c][0]);
    double vol2 = mesh1->cell_volume(c);
    gcl_err = std::fabs(vol1 - vol2);
  }

  // parallel collective operations
  double err_in[5] = {pl2_err, area, mass1, ql2_err, l20_err};
  double err_out[5];
  mesh1->Comm()->SumAll(err_in, err_out, 5);

  double err_tmp = pinf_err;
  mesh1->Comm()->MaxAll(&err_tmp, &pinf_err, 1);

  err_tmp = qinf_err;
  mesh1->Comm()->MaxAll(&err_tmp, &qinf_err, 1);

  err_tmp = inf0_err;
  mesh1->Comm()->MaxAll(&err_tmp, &inf0_err, 1);

  err_tmp = gcl_err;
  mesh1->Comm()->MaxAll(&err_tmp, &gcl_err, 1);

  // error tests
  pl2_err = std::pow(err_out[0], 0.5);
  ql2_err = std::pow(err_out[3], 0.5);
  l20_err = std::pow(err_out[4], 0.5);
  CHECK(pl2_err < 0.12 / (order + 1));

  if (MyPID == 0) {
    printf("nx=%3d  L2=%12.8g %12.8g %12.8g  Inf=%12.8g %12.8g %12.8g dMass=%10.4g  dArea=%10.6g\n", 
        nx, l20_err, pl2_err, ql2_err, 
            inf0_err, pinf_err, qinf_err, err_out[2] - mass0, 1.0 - err_out[1]);
    printf("GCL:    L1=%12.8g\n", gcl_err);
  }

  // initialize I/O
  Teuchos::ParameterList iolist;
  iolist.get<std::string>("file name base", "plot");
  OutputXDMF io(iolist, mesh1, true, false);

  io.InitializeCycle(t, nstep);
  io.WriteVector(*p2c(0), "remapped");
  io.WriteVector(*q2c(0), "remapped-prj");
  io.WriteVector(*p2c_err(0), "error");
  io.FinalizeCycle();
}

TEST(REMAP_DUAL_2D) {
  double dT(0.1);
  auto rk_method = Amanzi::Explicit_TI::heun_euler;
  RemapTestsDualRK(rk_method, "FEM", "", 10,10,0, dT);

  RemapTestsDualRK(rk_method, "VEM", "test/median15x16.exo", 0,0,0, dT/2);
  // RemapTestsDualRK(rk_method, "VEM", "", 5,5,5, dT);

  /*
  double dT(0.02);
  auto rk_method = Amanzi::Explicit_TI::tvd_3rd_order;
  std::string maps = "VEM";
  int deform = 5;
  RemapTestsDualRK(rk_method, maps, "",  16, 16,0, dT,    deform);
  RemapTestsDualRK(rk_method, maps, "",  32, 32,0, dT/2,  deform);
  RemapTestsDualRK(rk_method, maps, "",  64, 64,0, dT/4,  deform);
  RemapTestsDualRK(rk_method, maps, "", 128,128,0, dT/8,  deform);
  RemapTestsDualRK(rk_method, maps, "", 256,256,0, dT/16, deform);
  */

  /*
  double dT(0.02);
  auto rk_method = Amanzi::Explicit_TI::tvd_3rd_order;
  std::string maps = "VEM";
  int deform = 5;
  RemapTestsDualRK(rk_method, maps, "test/median15x16.exo",    16,0,0, dT,   deform);
  RemapTestsDualRK(rk_method, maps, "test/median32x33.exo",    32,0,0, dT/2, deform);
  RemapTestsDualRK(rk_method, maps, "test/median63x64.exo",    64,0,0, dT/4, deform);
  RemapTestsDualRK(rk_method, maps, "test/median127x128.exo", 128,0,0, dT/8, deform);
  RemapTestsDualRK(rk_method, maps, "test/median255x256.exo", 256,0,0, dT/16,deform);
  */

  /*
  double dT(0.05);
  auto rk_method = Amanzi::Explicit_TI::tvd_3rd_order;
  std::string maps = "PEM";
  RemapTestsDualRK(rk_method, maps, "test/mesh_poly20x20.exo",    20,0,0, dT);
  RemapTestsDualRK(rk_method, maps, "test/mesh_poly40x40.exo",    40,0,0, dT/2);
  RemapTestsDualRK(rk_method, maps, "test/mesh_poly80x80.exo",    80,0,0, dT/4);
  RemapTestsDualRK(rk_method, maps, "test/mesh_poly160x160.exo", 160,0,0, dT/8);
  */

  /*
  double dT(0.05);
  auto rk_method = Amanzi::Explicit_TI::tvd_3rd_order;
  std::string maps = "PEM";
  RemapTestsDualRK(rk_method, maps, "test/random10.exo", 10,0,0, dT);
  RemapTestsDualRK(rk_method, maps, "test/random20.exo", 20,0,0, dT/2);
  RemapTestsDualRK(rk_method, maps, "test/random40.exo", 40,0,0, dT/4);
  */

  /*
  double dT(0.025);
  auto rk_method = rk_method;
  std::string maps = "PEM";
  RemapTestsDualRK(rk_method, maps, "test/triangular8.exo",  0,0,0, dT);
  RemapTestsDualRK(rk_method, maps, "test/triangular16.exo", 0,0,0, dT/2);
  RemapTestsDualRK(rk_method, maps, "test/triangular32.exo", 0,0,0, dT/4);
  RemapTestsDualRK(rk_method, maps, "test/triangular64.exo", 0,0,0, dT/8);
  RemapTestsDualRK(rk_method, maps, "test/triangular128.exo",0,0,0, dT/16);
  */
}

