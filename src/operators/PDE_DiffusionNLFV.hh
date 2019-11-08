/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy (dasvyat@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)

  PDE_DiffusionNLFV implements the PDE_Diffusion interface 
  using nonlinear finite volumes.
*/

#ifndef AMANZI_OPERATOR_PDE_DIFFUSION_NLFV_HH_
#define AMANZI_OPERATOR_PDE_DIFFUSION_NLFV_HH_

#include <string>
#include <vector>

// TPLs
#include "Epetra_IntVector.h"
#include "Ifpack.h" 
#include "Teuchos_RCP.hpp"

// Amanzi
#include "CompositeVector.hh"
#include "DenseMatrix.hh"
#include "Preconditioner.hh"

#include "PDE_Diffusion.hh"

namespace Amanzi {
namespace Operators {

class BCs;

class PDE_DiffusionNLFV : public PDE_Diffusion {
 public:
  PDE_DiffusionNLFV(Teuchos::ParameterList& plist,
                    const Teuchos::RCP<Operator>& global_op) :
      PDE_Diffusion(plist, global_op),
      stencil_initialized_(false)
  {}

  PDE_DiffusionNLFV(Teuchos::ParameterList& plist,
                    const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
      PDE_Diffusion(plist, mesh),
      stencil_initialized_(false)
  {}

  virtual void Init() override;
  
  // main virtual members
  // -- setup
  using PDE_Diffusion::Setup;
  virtual void SetScalarCoefficient(const Teuchos::RCP<const CompositeVector>& k,
                                    const Teuchos::RCP<const CompositeVector>& dkdp) override;

  // -- create an operator
  virtual void UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& flux,
                              const Teuchos::Ptr<const CompositeVector>& u) override;
  virtual void UpdateMatricesNewtonCorrection(
      const Teuchos::Ptr<const CompositeVector>& flux,
      const Teuchos::Ptr<const CompositeVector>& u, double scalar_factor) override;
  
  virtual void UpdateMatricesNewtonCorrection(
      const Teuchos::Ptr<const CompositeVector>& flux,
      const Teuchos::Ptr<const CompositeVector>& u,
      const Teuchos::Ptr<const CompositeVector>& factor) override;
  
  // -- after solving the problem: postrocessing
  virtual void UpdateFlux(const Teuchos::Ptr<const CompositeVector>& u,
                          const Teuchos::Ptr<CompositeVector>& flux) override;

  // -- modify an operator
  virtual void ApplyBCs(bool primary, bool eliminate, bool essential_eqn) override;

 protected:
  // -- solution can be modified on boundary faces. This reflects specifics
  //    of nonlinear FV schemes, see implementation in the derived classes.
  virtual double MapBoundaryValue_(int f, double u) { return u; }

 protected:
  void InitStencils_();
  void OneSidedFluxCorrections_(int i0, const CompositeVector& u, CompositeVector& sideflux);
  void OneSidedWeightFluxes_(int i0, const CompositeVector& u, CompositeVector& sideflux);
  int OrderCellsByGlobalId_(const AmanziMesh::Entity_ID_List& cells, int& c1, int& c2);
  
 protected:
  int dim_;
  int newton_correction_;

  bool stencil_initialized_;
  Teuchos::RCP<CompositeVector> stencil_data_;
  std::vector<Teuchos::RCP<Epetra_IntVector> > stencil_faces_;
  std::vector<Teuchos::RCP<Epetra_IntVector> > stencil_cells_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif
