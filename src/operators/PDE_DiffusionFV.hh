/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy (dasvyat@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)

  DiffusionFV implements the Diffusion interface using
  finite volumes, i.e. the two point flux approximation.


  NOTE on the mesh requirements.
  ------------------------------
  It needs a limited set of the mesh interface, and therefore can be
  defined on things less "mesh-like" and more topological.  To
  facilitate that, the needed mesh interface is:

    - space_dimension()
    - getNumEntities(CELL,FACE,NODE)
    - face_get_cells()
    - cell_get_faces_and_dirs()
    - cell_map()
    - getFaceArea()
    - face_normal()
    - getFaceCentroid()
    - getCellCentroid()
   
    NOTE: actually, cell-to-cell distance, face-to-cell distance, not
    necessarily centroid locations are necessary, but this is not in
    the current mesh interface.
*/

#ifndef AMANZI_OPERATOR_PDE_DIFFUSION_FV_HH_
#define AMANZI_OPERATOR_PDE_DIFFUSION_FV_HH_

#include <strings.h>

// TPLs
#include "Teuchos_RCP.hpp"

// Amanzi
#include "CompositeVector.hh"
#include "DenseMatrix.hh"
#include "Preconditioner.hh"

// Opertors
#include "PDE_Diffusion.hh"

namespace Amanzi {
namespace Operators {

class BCs;

class PDE_DiffusionFV : public virtual PDE_Diffusion {
 public:
  PDE_DiffusionFV(Teuchos::ParameterList& plist,
                  const Teuchos::RCP<Operator>& global_op) :
      PDE_Diffusion(plist, global_op),
      transmissibility_initialized_(false)
  {}

  PDE_DiffusionFV(Teuchos::ParameterList& plist,
                  const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
      PDE_Diffusion(plist, mesh),
      transmissibility_initialized_(false)
  {}

  virtual void Init() override;
  
  // main virtual members
  // -- setup
  virtual void SetTensorCoefficient(const Teuchos::RCP<const TensorVector>& K) override;
  virtual void SetScalarCoefficient(const Teuchos::RCP<const CompositeVector>& k,
                                    const Teuchos::RCP<const CompositeVector>& dkdp) override;

  // -- create an operator
  virtual void UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& flux,
                              const Teuchos::Ptr<const CompositeVector>& u) override;
  virtual void UpdateMatricesNewtonCorrection(
      const Teuchos::Ptr<const CompositeVector>& flux,
      const Teuchos::Ptr<const CompositeVector>& u,
      double scalar_factor = 1.0) override;

  virtual void UpdateMatricesNewtonCorrection(
      const Teuchos::Ptr<const CompositeVector>& flux,
      const Teuchos::Ptr<const CompositeVector>& u,
      const Teuchos::Ptr<const CompositeVector>& factor) override;
  
  virtual void UpdateFlux(const Teuchos::Ptr<const CompositeVector>& u,
                          const Teuchos::Ptr<CompositeVector>& flux) override;

  // -- modify an operator
  virtual void ApplyBCs(bool primary, bool eliminate, bool essential_eqn) override;
  virtual void ApplyBCsJacobian() override;
  
  virtual void ModifyMatrices(const CompositeVector& u) override {};
  virtual void ScaleMassMatrices(double s) override {};

  // access
  const CompositeVector& transmissibility() { return *transmissibility_; }

  virtual CompositeVectorSpace scalar_coefficient_derivative_space() const override
  {
    CompositeVectorSpace out;
    out.SetMesh(mesh_);
    out.SetGhosted();
    if (little_k_type_ != OPERATOR_LITTLE_K_NONE) {
      out.AddComponent("cell", AmanziMesh::CELL, 1);
    }
    return out;
  }

  
 protected:
  cMultiVectorView_type_<DefaultDevice,double>
  ScalarCoefficientFaces(bool scatter) const {
    if (k_ != Teuchos::null) {
      if (scatter) k_->scatterMasterToGhosted("face");
      if (k_->hasComponent("face")) {
        return k_->viewComponent<DefaultDevice>("face", true);
      }
    }
    MultiVectorView_type_<DefaultDevice,double> k_face("k_face", nfaces_wghost, 1);
    Kokkos::deep_copy(k_face, 1.0);
    return k_face;
  }

public: 
  // This function need to be public for Kokkos Lambda
  void ComputeTransmissibility_();
  virtual void AnalyticJacobian_(const CompositeVector& solution);

  // virtual void ComputeJacobianLocal_(
  //     int mcells, int f, int face_dir_0to1, int bc_model, double bc_value,
  //     double *pres, double *dkdp_cell, WhetStone::DenseMatrix& Jpp);
  
 protected:
  Teuchos::RCP<CompositeVector> transmissibility_;
  bool transmissibility_initialized_;
};

}  // namespace Operators
}  // namespace Amanzi


#endif
