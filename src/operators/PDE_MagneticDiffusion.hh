/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Magnetic diffusion exesices Amanzi's capability to discretize 
  curl operators.
*/

#ifndef AMANZI_OPERATOR_MAGNETIC_DIFFUSION_HH_
#define AMANZI_OPERATOR_MAGNETIC_DIFFUSION_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

// Amanzi
#include "CompositeVector.hh"
#include "DenseMatrix.hh"
#include "exceptions.hh"
#include "Tensor.hh"

// Amanzi::Operators
#include "BCs.hh"
#include "PDE_Electromagnetics.hh"
#include "Operator.hh"
#include "OperatorDefs.hh"

namespace Amanzi {
namespace Operators {

class PDE_MagneticDiffusion : public PDE_Electromagnetics {
 public:
  PDE_MagneticDiffusion(const Teuchos::RCP<Operator>& global_op)
    : PDE_Electromagnetics(global_op)
  {};

  PDE_MagneticDiffusion(Teuchos::ParameterList& plist,
                          const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
    : PDE_Electromagnetics(plist, mesh)
  {
    pde_type_ = PDE_MAGNETIC_DIFFUSION;
    InitMagneticDiffusion_(plist);
  }

  // main virtual members
  // -- create a linearized operator
  using PDE_HelperDiscretization::UpdateMatrices;
  virtual void UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& u,
                              const Teuchos::Ptr<const CompositeVector>& p) override;

  // -- before solving the problem
  virtual void ModifyMatrices(CompositeVector& E, CompositeVector& B, double dt) override;

  // -- after solving the problem
  virtual void ModifyFields(CompositeVector& E, CompositeVector& B, double dt) override;

  // physical quantities
  // -- energies
  virtual double CalculateOhmicHeating(const CompositeVector& E);
  double CalculateMagneticEnergy(const CompositeVector& B);

  // -- divergence
  double CalculateDivergence(int c, const CompositeVector& B);

 private:
  void InitMagneticDiffusion_(Teuchos::ParameterList& plist);

 protected:
  std::vector<WhetStone::DenseMatrix> mass_op_;
  std::vector<WhetStone::DenseMatrix> curl_op_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif


