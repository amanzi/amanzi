/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Mini classes implement mathematical models for special physics, such 
  as serial 1D dual porosity models. 
*/

#ifndef AMANZI_MINI_DIFFUSION_1D_HH_
#define AMANZI_MINI_DIFFUSION_1D_HH_

#include <memory>

#include <DenseVector.hh>

#include <Mini_Operator1D.hh>

namespace Amanzi {
namespace Operators {

class Mini_Diffusion1D : public Mini_Operator1D {
 public:
  Mini_Diffusion1D() : Kconst_(1.0) {};
  ~Mini_Diffusion1D() {};

  // set up operator
  void Setup(double K) { Kconst_ = K; }
  void Setup(const std::shared_ptr<WhetStone::DenseVector> k,
             const std::shared_ptr<WhetStone::DenseVector> dkdp) {
    k_ = k;
    dkdp_ = dkdp;
  }
  void Setup(const std::shared_ptr<WhetStone::DenseVector> K,
             const std::shared_ptr<WhetStone::DenseVector> k,
             const std::shared_ptr<WhetStone::DenseVector> dkdp) {
    K_ = K;
    k_ = k;
    dkdp_ = dkdp;
  }

  // generate linearized operators
  // -- build phisical model
  void UpdateMatrices();
  // -- build Jacobian
  void UpdateJacobian(const WhetStone::DenseVector& p,
                      double bcl, int type_l, double bcr, int type_r);

  // modify matrix due to boundary conditions
  void ApplyBCs(double bcl, int type_l, double bcr, int type_r);

  // postprocessing
  // -- flux calculation uses potential p to calculate flux u
  // void UpdateFlux(const WhetStone::DenseVector& p, WhetStone::DenseVector& u);

  // access
  WhetStone::DenseVector& k() { return *k_; }
  WhetStone::DenseVector& dkdp() { return *dkdp_; }

 private:
  std::shared_ptr<WhetStone::DenseVector> K_, k_, dkdp_;
  double Kconst_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif


