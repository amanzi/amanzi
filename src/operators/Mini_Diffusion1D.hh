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

namespace Amanzi {
namespace Operators {

class Mini_Diffusion1D {
 public:
  Mini_Diffusion1D(std::shared_ptr<const WhetStone::DenseVector> mesh);
  ~Mini_Diffusion1D() {};

  // set up operator
  void Setup(const std::shared_ptr<WhetStone::DenseVector> K,
             const std::shared_ptr<WhetStone::DenseVector> k,
             const std::shared_ptr<WhetStone::DenseVector> dkdp) {
    K_ = K;
    k_ = k;
    dkdp_ = dkdp;
  }

  // generate linearized operator
  void UpdateMatrices();

  // modify matrix due to boundary conditions
  void SetBCs();
  void ApplyBCs(bool primary, bool eliminate, bool essential_eqn);

  // postprocessing
  // -- flux calculation uses potential p to calculate flux u
  void UpdateFlux(const WhetStone::DenseVector& p, WhetStone::DenseVector& u);

 protected:
  std::shared_ptr<const WhetStone::DenseVector> mesh_;
  std::shared_ptr<WhetStone::DenseVector> K_, k_, dkdp_;

 private:
  WhetStone::DenseVector diag_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif


