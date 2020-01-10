//! Diffusion generates local Ops and global Operators for an elliptic operator.

/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
          Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_OPERATOR_PDE_DIFFUSION_HH_
#define AMANZI_OPERATOR_PDE_DIFFUSION_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "exceptions.hh"
#include "Tensor.hh"
#include "Point.hh"
#include "CompositeVector.hh"
#include "DenseMatrix.hh"

#include "Operator.hh"
#include "OperatorDefs.hh"
#include "PDE_HelperDiscretization.hh"

#include <Kokkos_Core.hpp>
#include <Kokkos_Vector.hpp>

/*!
Example:

.. code-block:: xml

    <ParameterList name="OPERATOR_NAME">
      <Parameter name="discretization primary" type="string" value="mfd: optimized for monotonicity"/>
      <Parameter name="discretization secondary" type="string" value="mfd: two-point flux approximation"/>
      <Parameter name="schema" type="Array(string)" value="{face, cell}"/>
      <Parameter name="preconditioner schema" type="Array(string)" value="{face}"/>
      <Parameter name="gravity" type="bool" value="true"/>
      <Parameter name="gravity term discretization" type="string" value="hydraulic head"/>
      <Parameter name="nonlinear coefficient" type="string" value="upwind: face"/>
      <Parameter name="Newton correction" type="string" value="true Jacobian"/>

      <ParameterList name="consistent faces">
        <ParameterList name="linear solver">
          ...
        </ParameterList>
        <ParameterList name="preconditioner">
          ...
        </ParameterList>
      </ParameterList>
    </ParameterList>
*/


/*
  Ghost elemets of composite vectors k_ and dkdp_ are NOT used to be up to date.
  They always have to be updated from master elements before accessing their values.
*/

namespace Amanzi {
namespace Operators {

class PDE_Diffusion : public PDE_HelperDiscretization {
 public:
  PDE_Diffusion(Teuchos::ParameterList& plist,
                const Teuchos::RCP<Operator>& global_op) :
      PDE_HelperDiscretization(global_op),
      plist_(plist)
  {};

  PDE_Diffusion(Teuchos::ParameterList& plist,
                const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
      PDE_HelperDiscretization(mesh),
      plist_(plist)
  {};

  virtual ~PDE_Diffusion() = default;

  // virtual constructor, creates the space/schema and operator.
  virtual void Init() = 0;
  
  // Setters and Setup
  //
  // Note that these default setters can be overridden to do actual work.
  virtual void SetTensorCoefficient(const Kokkos::vector<WhetStone::Tensor>& K) {
    K_ = K;
  }
  virtual void SetScalarCoefficient(const Teuchos::RCP<const CompositeVector>& k,
          const Teuchos::RCP<const CompositeVector>& dkdp) {
    k_ = k;
    dkdp_ = dkdp;
  }
  // Note that gravity and density can be ignored in non-gravity-affected
  // diffusion.
  virtual void SetGravity(const AmanziGeometry::Point& g) {
    g_ = g;
  }
  virtual void SetDensity(double rho) {
    is_scalar_ = true;
    rho_ = rho;
  }
  virtual void SetDensity(const Teuchos::RCP<const CompositeVector>& rho) {
    is_scalar_ = false;
    if (rho->HasComponent("cell")) {
      rho_cv_ = rho;
    }
  }

  // NOTE: this is very inefficient! Deprecate this!
  // 
  // double GetDensity(const int c) {
  //   return is_scalar_ ? rho_ :
  //   if (is_scalar_) {
  //     return rho_ ;
  //   }else{
  //     return (*rho_cv_->ViewComponent("cell", true))[0][c];
  //   }
  // }
  
  // Lumped Setters for lazy development
  void Setup(const Kokkos::vector<WhetStone::Tensor>& K,
             const Teuchos::RCP<const CompositeVector>& k,
             const Teuchos::RCP<const CompositeVector>& dkdp) {
    SetTensorCoefficient(K);
    SetScalarCoefficient(k, dkdp);
  }
  void Setup(const Kokkos::vector<WhetStone::Tensor> & K,
             const Teuchos::RCP<const CompositeVector>& k,
             const Teuchos::RCP<const CompositeVector>& dkdp,
             const double rho,
             const AmanziGeometry::Point& g) {
    SetTensorCoefficient(K);
    SetScalarCoefficient(k, dkdp);
    SetDensity(rho);
    SetGravity(g);
  }
  void Setup(const Kokkos::vector<WhetStone::Tensor>& K,
             const Teuchos::RCP<const CompositeVector>& k,
             const Teuchos::RCP<const CompositeVector>& dkdp,
             const Teuchos::RCP<const CompositeVector>& rho,
             const AmanziGeometry::Point& g) {
    SetTensorCoefficient(K);
    SetScalarCoefficient(k, dkdp);
    SetDensity(rho);
    SetGravity(g);
  }

  // generate linearized operator
  // -- generate local matrices. We can use parameter to define coefficeints
  //    or/and perform on-a-fly linearization. 
  virtual void UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& u,
          const Teuchos::Ptr<const CompositeVector>& p) = 0;
  virtual void UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& u) {
    UpdateMatrices(u, Teuchos::null);
  }
  virtual void UpdateMatrices() {
    UpdateMatrices(Teuchos::null, Teuchos::null);
  }
  
  // -- populate additional Jacobian local matrices
  virtual void UpdateMatricesNewtonCorrection(
          const Teuchos::Ptr<const CompositeVector>& flux,
          const Teuchos::Ptr<const CompositeVector>& u,
          double scalar_factor = 1.0) = 0;
  
  virtual void UpdateMatricesNewtonCorrection(
          const Teuchos::Ptr<const CompositeVector>& flux,
          const Teuchos::Ptr<const CompositeVector>& u,
          const Teuchos::Ptr<const CompositeVector>& factor) = 0;

  // postprocessing
  // -- flux calculation uses potential p to calculate flux u
  //
  // NOTE: this probably should live in Diffusion?  This is a very
  // Diffusion-looking interface that is not valid for other operators?
  virtual void UpdateFlux(const Teuchos::Ptr<const CompositeVector>& p,
                          const Teuchos::Ptr<CompositeVector>& u) = 0;
  
  // -- additional interface on non-manifolds
  virtual void UpdateFluxNonManifold(const Teuchos::Ptr<const CompositeVector>& u,
          const Teuchos::Ptr<CompositeVector>& flux) {
    Errors::Message msg("Diffusion: This diffusion implementation does not support non-manifolds.");
    Exceptions::amanzi_throw(msg);
  }

  // -- matrix modifications
  virtual void ModifyMatrices(const CompositeVector& u) {
    Errors::Message msg("Diffusion: This diffusion implementation does not support ModifyMatrices.");
    Exceptions::amanzi_throw(msg);
  }
  virtual void ScaleMassMatrices(double s) {
    Errors::Message msg("Diffusion: This diffusion implementation does not support ScaleMassMatrices.");
    Exceptions::amanzi_throw(msg);
  }
    

  // -- working with consistent faces -- may not be implemented
  virtual int UpdateConsistentFaces(CompositeVector& u) {
    Errors::Message msg("Diffusion: This diffusion implementation does not support working with consistent faces.");
    Exceptions::amanzi_throw(msg);
    return 1;
  }
  
  // interface to solvers for treating nonlinear BCs.
  virtual double ComputeTransmissibility(int f) const {
    Errors::Message msg("Diffusion: This diffusion implementation does not support local calculations.");
    Exceptions::amanzi_throw(msg);
    return 1;
  }
    
  virtual double ComputeGravityFlux(int f) const {
    Errors::Message msg("Diffusion: This diffusion implementation does not support local calculations.");
    Exceptions::amanzi_throw(msg);
    return 1;
  }
    

  // access
  int schema_prec_dofs() { return global_op_schema_; }
  int schema_dofs() { return local_op_schema_; }

  Teuchos::RCP<const Op> jacobian_op() const { return jac_op_; }
  Teuchos::RCP<Op> jacobian_op() { return jac_op_; }
  void set_jacobian_op(const Teuchos::RCP<Op>& op);
  int schema_jacobian() { return jac_op_schema_; }

  int little_k() const { return little_k_; }
  CompositeVectorSpace little_k_space() const {
    CompositeVectorSpace out;
    out.SetMesh(mesh_);
    out.SetGhosted();
    if (little_k_ == OPERATOR_LITTLE_K_NONE) {
      return out;
    }
    if (little_k_ != OPERATOR_LITTLE_K_UPWIND) {
      out.AddComponent("cell", AmanziMesh::CELL, 1);
    }
    if (little_k_ != OPERATOR_LITTLE_K_STANDARD) {
      out.AddComponent("face", AmanziMesh::FACE, 1);
    }
    if (little_k_ == OPERATOR_LITTLE_K_DIVK_TWIN || 
        little_k_ == OPERATOR_LITTLE_K_DIVK_TWIN_GRAD) {
      out.AddComponent("twin", AmanziMesh::FACE, 1);
    }
    if (little_k_ == OPERATOR_LITTLE_K_DIVK_TWIN_GRAD) {
      out.AddComponent("grad", AmanziMesh::CELL, mesh_->space_dimension());
    }
    return out;          
  }

 protected:
  Teuchos::ParameterList plist_;
  Kokkos::vector<WhetStone::Tensor> K_; 
  //Teuchos::RCP<const std::vector<WhetStone::Tensor> > K_;
  bool K_symmetric_;

  // nonlinear coefficient and its representation
  Teuchos::RCP<const CompositeVector> k_, dkdp_;
  int little_k_;

  // gravity
  bool is_scalar_;
  double rho_;
  Teuchos::RCP<const CompositeVector> rho_cv_;
  AmanziGeometry::Point g_;
  
  // additional operators
  Teuchos::RCP<Op> jac_op_;
  int global_op_schema_, local_op_schema_, jac_op_schema_;
};


/* ******************************************************************
* Initialization of the operator, scalar coefficient.
****************************************************************** */
inline
void PDE_Diffusion::set_jacobian_op(const Teuchos::RCP<Op>& op)
{
  if (global_operator().get()) {
    if (local_op().get()) {
      auto index = std::find(global_operator()->begin(), global_operator()->end(), jac_op_)
                   - global_operator()->begin();
      if (index != global_operator()->size())
        global_operator()->OpPushBack(op);
      else
        global_operator()->OpReplace(op, index);
    } else {
      global_operator()->OpPushBack(op);
    }
  }
  jac_op_ = op;
}

}  // namespace Operators
}  // namespace Amanzi

#endif


