/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Base class for remap methods.
*/

#ifndef AMANZI_OPERATOR_REMAP_DG_TESTS_HH_
#define AMANZI_OPERATOR_REMAP_DG_TESTS_HH_

#include "RemapDG.hh"

namespace Amanzi {

template<class AnalyticDG>
class RemapDG_Tests : public RemapDG {
 public:
  RemapDG_Tests(const Teuchos::RCP<const AmanziMesh::Mesh> mesh0,
                const Teuchos::RCP<AmanziMesh::Mesh> mesh1,
                Teuchos::ParameterList& plist) 
    : RemapDG(mesh0, mesh1, plist),
      nfun_(0),
      tprint_(0.0),
      l2norm_(1.0),
      dt_output_(0.1) {};
  ~RemapDG_Tests() {};

  // main function  
  virtual void FunctionalTimeDerivative(double t, const CompositeVector& u, CompositeVector& f) override;

  // mesh deformation
  virtual void DeformMesh(int deform, double t);
  AmanziGeometry::Point DeformNode(int deform, double t, const AmanziGeometry::Point& yv,
                                   const AmanziGeometry::Point& rv = AmanziGeometry::Point(3));

  // output
  virtual double global_time(double t) { return t; }
  void set_dt_output(double dt) { dt_output_ = dt; }

 protected:
  // statistics
  int nfun_;
  double tprint_, dt_output_, l2norm_;
};


/* *****************************************************************
* Main routine: evaluation of functional
***************************************************************** */
template<class AnalyticDG>
void RemapDG_Tests<AnalyticDG>::FunctionalTimeDerivative(
    double t, const CompositeVector& u, CompositeVector& f)
{
  RemapDG::FunctionalTimeDerivative(t, u, f);

  // statistics
  nfun_++;
  Epetra_MultiVector& xc = *field_->ViewComponent("cell");
  int nk = xc.NumVectors();
  double xmax[nk], xmin[nk];
  xc.MaxValue(xmax);
  xc.MinValue(xmin);

  double tglob = global_time(t);
  if (fabs(tprint_ - tglob) < 1e-6 && mesh0_->get_comm()->MyPID() == 0) {
    printf("t=%8.5f  L2=%9.5g  nfnc=%5d  sharp=%5.1f%%  umax: ", tglob, l2norm_, nfun_, sharp_);
    for (int i = 0; i < std::min(nk, 4); ++i) printf("%9.5g ", xmax[i]);
    printf("\n");
    tprint_ += dt_output_;
    sharp_ = 0.0;
  } 
}


/* *****************************************************************
* Deform mesh1
***************************************************************** */
template<class AnalyticDG>
void RemapDG_Tests<AnalyticDG>::DeformMesh(int deform, double t)
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
AmanziGeometry::Point RemapDG_Tests<AnalyticDG>::DeformNode(
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
  else if (deform == 6) {
    double phi = t * 2 * M_PI;
    double cs(std::cos(phi)), sn(std::sin(phi));
    xv[0] = cs * yv[0] - sn * yv[1];
    xv[1] = sn * yv[0] + cs * yv[1];
  }

  return xv;
}

} // namespace Amanzi

#endif
