/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy (dasvyat@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)

  OperatorDiffusionFV implements the OperatorDiffusion interface using
  finite volumes, i.e. the two point flux approximation.


  NOTE on the mesh requirements.
  ------------------------------
  It needs a limited set of the mesh interface, and therefore can be
  defined on things less "mesh-like" and more topological.  To
  facilitate that, the needed mesh interface is:

    - space_dimension()
    - num_entities(CELL,FACE,NODE)
    - face_get_cells()
    - cell_get_faces_and_dirs()
    - cell_map()
    - face_area()
    - face_normal()
    - face_centroid()
    - cell_centroid()
   
    NOTE: actually, cell-to-cell distance, face-to-cell distance, not
    necessarily centroid locations are necessary, but this is not in
    the current mesh interface.
*/

#ifndef AMANZI_OPERATOR_DIFFUSION_FV_HH_
#define AMANZI_OPERATOR_DIFFUSION_FV_HH_

#include <strings.h>

// TPLs
#include "Ifpack.h" 
#include "Teuchos_RCP.hpp"

// Amanzi
#include "CompositeVector.hh"
#include "DenseMatrix.hh"
#include "Preconditioner.hh"
#include "OperatorDiffusion.hh"

namespace Amanzi {
namespace Operators {

class BCs;

class OperatorDiffusionFV : public virtual OperatorDiffusion {
 public:
  OperatorDiffusionFV(Teuchos::ParameterList& plist,
                      const Teuchos::RCP<Operator>& global_op) :
      OperatorDiffusion(global_op),
      transmissibility_initialized_(false)
  {
    operator_type_ = OPERATOR_DIFFUSION_FV;
    InitDiffusion_(plist);
  }

  OperatorDiffusionFV(Teuchos::ParameterList& plist,
                      const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
      OperatorDiffusion(mesh),
      transmissibility_initialized_(false)
  {
    operator_type_ = OPERATOR_DIFFUSION_FV;
    InitDiffusion_(plist);
  }

  OperatorDiffusionFV(Teuchos::ParameterList& plist,
                      const Teuchos::RCP<AmanziMesh::Mesh>& mesh) :
      OperatorDiffusion(mesh),
      transmissibility_initialized_(false)
  {
    operator_type_ = OPERATOR_DIFFUSION_FV;
    InitDiffusion_(plist);
  }

  // main virtual members
  // -- setup
  using OperatorDiffusion::Setup;
  virtual void SetTensorCoefficient(const Teuchos::RCP<std::vector<WhetStone::Tensor> >& K);
  virtual void SetScalarCoefficient(const Teuchos::RCP<const CompositeVector>& k,
                                    const Teuchos::RCP<const CompositeVector>& dkdp);

  // -- create an operator
  virtual void UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& flux,
                              const Teuchos::Ptr<const CompositeVector>& u);
  virtual void UpdateMatricesNewtonCorrection(
      const Teuchos::Ptr<const CompositeVector>& flux,
      const Teuchos::Ptr<const CompositeVector>& u,
      double scalar_limiter=1.);
  virtual void UpdateFlux(const CompositeVector& u, CompositeVector& flux);

  // -- modify an operator
  virtual void ApplyBCs(bool primary, bool eliminate);
  virtual void ModifyMatrices(const CompositeVector& u) {};
  virtual void ScaleMassMatrices(double s) {};

  // Developments
  // -- interface to solvers for treating nonlinear BCs.
  virtual double ComputeTransmissibility(int f) const;
  virtual double ComputeGravityFlux(int f) const { return 0.0; }

  // access
  const CompositeVector& transmissibility() { return *transmissibility_; }

 protected:
  void ComputeTransmissibility_();

  void AnalyticJacobian_(const CompositeVector& solution);

  virtual void ComputeJacobianLocal_(
      int mcells, int f, int bc_model, double bc_value,
      double *pres, double *dkdp_cell, WhetStone::DenseMatrix& Jpp);

  virtual void InitDiffusion_(Teuchos::ParameterList& plist);
  
 protected:
  Teuchos::RCP<CompositeVector> transmissibility_;
  bool transmissibility_initialized_;

  int newton_correction_;
  bool exclude_primary_terms_;
};

}  // namespace Operators
}  // namespace Amanzi


#endif
