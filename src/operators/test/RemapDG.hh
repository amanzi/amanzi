/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Base class for remap methods.
*/

#ifndef AMANZI_REMAP_DG_HH_
#define AMANZI_REMAP_DG_HH_

#include <cstdlib>
#include <iostream>
#include <vector>

// TPLs
#include "Epetra_MultiVector.h"
#include "Teuchos_RCP.hpp"

// Amanzi
#include "Explicit_TI_RK.hh"
#include "Mesh.hh"
#include "MeshMapsFactory.hh"
#include "VectorPolynomial.hh"
#include "WhetStoneDefs.hh"

// Amanzi::Operators
#include "OperatorDefs.hh"
#include "PDE_Abstract.hh"
#include "PDE_AdvectionRiemann.hh"
#include "PDE_Reaction.hh"

namespace Amanzi {

template<class AnalyticDG>
class RemapDG : public Explicit_TI::fnBase<CompositeVector> {
 public:
  RemapDG(const Teuchos::RCP<const AmanziMesh::Mesh> mesh0,
          const Teuchos::RCP<AmanziMesh::Mesh> mesh1,
          Teuchos::ParameterList& plist) 
    : mesh0_(mesh0),
      mesh1_(mesh1),
      plist_(plist),
      dim_(mesh0->space_dimension()),
      high_order_velocity_(false),
      dt_output_(0.1) {};
  ~RemapDG() {};

  // main members required by high-level interface
  virtual void FunctionalTimeDerivative(double t, const CompositeVector& u, CompositeVector& f) override;

  // create basic structures
  virtual void Init() { 
    InitPrimary();
    InitSecondary();
  }
  void InitPrimary();
  void InitSecondary();

  // geometric tools
  // -- various geometric quantities
  virtual void Jacobian(int c, double t, const WhetStone::MatrixPolynomial& J,
                        WhetStone::MatrixPolynomial& Jt);
  virtual void UpdateGeometricQuantities(double t);

  // -- mesh deformation
  virtual void DeformMesh(int deform, double t);
  AmanziGeometry::Point DeformNode(int deform, double t, const AmanziGeometry::Point& yv,
                                              const AmanziGeometry::Point& rv = AmanziGeometry::Point(3));

  // output
  virtual double global_time(double t) { return t; }
  void set_dt_output(double dt) { dt_output_ = dt; }

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh0_;
  Teuchos::RCP<AmanziMesh::Mesh> mesh1_;
  int ncells_owned_, ncells_wghost_, nfaces_owned_, nfaces_wghost_;
  int dim_;

  const Teuchos::ParameterList plist_;
  std::shared_ptr<WhetStone::MeshMaps> maps_;
  bool high_order_velocity_;

  // operators
  int order_;
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
  double tprint_, dt_output_, l2norm_;
};


/* *****************************************************************
* Initialization of remap: operarot and face velocity.
***************************************************************** */
template<class AnalyticDG>
void RemapDG<AnalyticDG>::InitPrimary()
{
  // mesh data
  ncells_owned_ = mesh0_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  ncells_wghost_ = mesh0_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
  nfaces_owned_ = mesh0_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  nfaces_wghost_ = mesh0_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  order_ = plist_.sublist("PK operator")
                 .sublist("flux operator").template get<int>("method order");

  // create right-hand side operator
  // -- flux
  auto oplist = plist_.sublist("PK operator").sublist("flux operator");
  op_flux_ = Teuchos::rcp(new Operators::PDE_AdvectionRiemann(oplist, mesh0_));
  auto global_op = op_flux_->global_operator();

  // -- advection
  oplist = plist_.sublist("PK operator").sublist("advection operator");
  op_adv_ = Teuchos::rcp(new Operators::PDE_Abstract(oplist, global_op));

  // create left-hand side operator
  oplist = plist_.sublist("PK operator").sublist("reaction operator");
  op_reac_ = Teuchos::rcp(new Operators::PDE_Reaction(oplist, mesh0_));

  // factory of mesh maps
  auto map_list = plist_.sublist("maps");
  WhetStone::MeshMapsFactory maps_factory;
  maps_ = maps_factory.Create(map_list, mesh0_, mesh1_);

  // boundary data
  auto bc = Teuchos::rcp(new Operators::BCs(mesh0_, AmanziMesh::FACE, Operators::DOF_Type::VECTOR));
  std::vector<int>& bc_model = bc->bc_model();
  std::vector<std::vector<double> >& bc_value = bc->bc_value_vector(1);

  WhetStone::Polynomial coefs;

  for (int f = 0; f < nfaces_wghost_; f++) {
    const AmanziGeometry::Point& xf = mesh0_->face_centroid(f);
    if (fabs(xf[0]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6 ||
        fabs(xf[1]) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6 ||
        fabs(xf[dim_ - 1]) < 1e-6 || fabs(xf[dim_ - 1] - 1.0) < 1e-6) {
      bc_model[f] = Operators::OPERATOR_BC_REMOVE;
    }
  }
  op_flux_->SetBCs(bc, bc);

  // face velocities fixed for the pseudo time interval [0, 1]
  velf_vec_.resize(nfaces_wghost_);
  if (high_order_velocity_) {
    int vel_order = plist_.sublist("maps").template get<int>("method order");
    AnalyticDG ana(mesh0_, vel_order, true);

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

  // memory allocation
  velf_ = Teuchos::rcp(new std::vector<WhetStone::Polynomial>(nfaces_wghost_));
  velc_ = Teuchos::rcp(new std::vector<WhetStone::VectorPolynomial>(ncells_wghost_));
  jac_ = Teuchos::rcp(new std::vector<WhetStone::VectorPolynomial>(ncells_wghost_));

  // miscallateous
  tprint_ = 0.0;
  nfun_ = 0;
  l2norm_ = -1.0;
}


/* *****************************************************************
* Initialization of remap: Jacobian and cell velcolity.
***************************************************************** */
template<class AnalyticDG>
void RemapDG<AnalyticDG>::InitSecondary()
{
  WhetStone::Entity_ID_List faces;
  J_.resize(ncells_wghost_);
  uc_.resize(ncells_wghost_);

  for (int c = 0; c < ncells_wghost_; ++c) {
    mesh0_->cell_get_faces(c, &faces);

    std::vector<WhetStone::VectorPolynomial> vvf;
    for (int n = 0; n < faces.size(); ++n) {
      vvf.push_back(velf_vec_[faces[n]]);
    }

    maps_->VelocityCell(c, vvf, uc_[c]);
    maps_->Jacobian(uc_[c], J_[c]);
  }
}


/* *****************************************************************
* Main routine: evaluation of functional
***************************************************************** */
template<class AnalyticDG>
void RemapDG<AnalyticDG>::FunctionalTimeDerivative(
    double t, const CompositeVector& u, CompositeVector& f) {
  UpdateGeometricQuantities(t);

  // -- populate operators
  op_adv_->SetupPolyVector(velc_);
  op_adv_->UpdateMatrices();

  op_flux_->Setup(velc_, velf_);
  op_flux_->UpdateMatrices(velf_.ptr());
  op_flux_->ApplyBCs(true, true, true);

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

  if (fabs(tprint_ - t) < 1e-6 && mesh0_->get_comm()->MyPID() == 0) {
    printf("t=%8.5f  L2=%9.5g  nfnc=%3d  umax: ", global_time(t), l2norm_, nfun_);
    for (int i = 0; i < std::min(nk, 3); ++i) printf("%9.5g ", xmax[i]);
    printf("\n");
    tprint_ += dt_output_;
  } 

  // -- calculate right-hand_side
  op_flux_->global_operator()->Apply(x, f);
}


/* *****************************************************************
* Calculates various geometric quantaties on intermediate meshes.
***************************************************************** */
template<class AnalyticDG>
void RemapDG<AnalyticDG>::Jacobian(
    int c, double t, const WhetStone::MatrixPolynomial& J, WhetStone::MatrixPolynomial& Jt)
{
  int nJ = J.NumRows();
  Jt = J * t;

  for (int i = 0; i < nJ; ++i) {
    Jt(i, i % dim_)(0) += 1.0;
  }
}


/* *****************************************************************
* Calculates various geometric quantaties on intermediate meshes.
***************************************************************** */
template<class AnalyticDG>
void RemapDG<AnalyticDG>::UpdateGeometricQuantities(double t)
{
  WhetStone::VectorPolynomial tmp, fmap, cn;

  for (int f = 0; f < nfaces_wghost_; ++f) {
    // cn = j J^{-t} N dA
    if (high_order_velocity_) {
      maps_->VelocityFace(f, tmp);
      fmap = t * tmp;
      maps_->NansonFormula(f, fmap, cn);
      (*velf_)[f] = tmp * cn;
    } else {
      fmap = t * velf_vec_[f];
      maps_->NansonFormula(f, fmap, cn);
      (*velf_)[f] = velf_vec_[f] * cn;
    }
  }

  WhetStone::MatrixPolynomial Jt, C;
  for (int c = 0; c < ncells_wghost_; ++c) {
    Jacobian(c, t, J_[c], Jt);
    maps_->Determinant(Jt, (*jac_)[c]);
    maps_->Cofactors(Jt, C);
    
    // cell-based pseudo velocity -C^t u 
    int nC = C.NumRows();
    (*velc_)[c].resize(nC);

    int kC = nC / dim_;
      for (int n = 0; n < kC; ++n) {
      int m = n * dim_;
      for (int i = 0; i < dim_; ++i) {
        (*velc_)[c][m + i].Reshape(dim_, 0, true);
        (*velc_)[c][m + i].set_origin(uc_[c][0].origin());

        for (int k = 0; k < dim_; ++k) {
          (*velc_)[c][m + i] -= C(m + k, i) * uc_[c][m + k];
        }
      }
    }
  }
}


/* *****************************************************************
* Deform mesh1
***************************************************************** */
template<class AnalyticDG>
void RemapDG<AnalyticDG>::DeformMesh(int deform, double t)
{
  // create distributed random vector
  CompositeVectorSpace cvs;
  cvs.SetMesh(mesh0_)->SetGhosted(true)->AddComponent("node", AmanziMesh::NODE, dim_);
  CompositeVector random(cvs);

  int gid = mesh0_->node_map(false).MaxAllGID();
  double scale = 0.2 * std::pow(gid, -1.0 / dim_);
  Epetra_MultiVector& random_n = *random.ViewComponent("node", true);

  random_n.Random();
  random_n.Scale(scale);
  random.ScatterMasterToGhosted();

  // relocate mesh nodes
  AmanziGeometry::Point xv(dim_), yv(dim_), uv(dim_), rv(dim_);
  AmanziMesh::Entity_ID_List nodeids;
  AmanziGeometry::Point_List new_positions, final_positions;

  int nnodes = mesh0_->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::ALL);

  for (int v = 0; v < nnodes; ++v) {
    for (int i = 0; i < dim_; ++i) rv[i] = random_n[i][v];
    mesh0_->node_get_coordinates(v, &xv);
    nodeids.push_back(v);
    new_positions.push_back(DeformNode(deform, t, xv, rv));
  }
  mesh1_->deform(nodeids, new_positions, false, &final_positions);
}


/* *****************************************************************
* Deformation functional
***************************************************************** */
template<class AnalyticDG>
AmanziGeometry::Point RemapDG<AnalyticDG>::DeformNode(
  int deform, double t, const AmanziGeometry::Point& yv, const AmanziGeometry::Point& rv)
{
  AmanziGeometry::Point uv(dim_), xv(yv);

  if (deform == 1) {
    double ds(0.0001);
    int n = t / ds;
    for (int i = 0; i < n; ++i) {
      if (dim_ == 2) {
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
    if (fabs(yv[0]) > 1e-6 && fabs(1.0 - yv[0]) > 1e-6 &&
        fabs(yv[1]) > 1e-6 && fabs(1.0 - yv[1]) > 1e-6) {
      xv[0] += rv[0];
      xv[1] += rv[1];
    }
  }
  else if (deform == 4) {
    xv[0] += t * yv[0] * yv[1] * (1.0 - yv[0]) / 2;
    xv[1] += t * yv[0] * yv[1] * (1.0 - yv[1]) / 2;
  }
  else if (deform == 5) {
    xv[0] += t * yv[0] * (1.0 - yv[0]) / 2;
    xv[1] += t * yv[1] * (1.0 - yv[1]) / 2;
  }

  return xv;
}

} // namespace Amanzi

#endif
