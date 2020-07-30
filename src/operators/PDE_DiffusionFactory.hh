// PDE_DiffusionFactory constructs objects which implement the interface for a PDE_Diffusion.

/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  This documentation is for the entire Diffusion concept, which is maintained 
  here because the input spec for Diffusion objects is defined/used here.
*/

#ifndef AMANZI_OPERATOR_PDE_DIFFUSION_FACTORY_HH_
#define AMANZI_OPERATOR_PDE_DIFFUSION_FACTORY_HH_

#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "PDE_Diffusion.hh"
#include "PDE_DiffusionWithGravity.hh"

/*!

``PDE_Diffusion`` forms local ``Op`` s and global ``Operator`` s for elliptic equations:

.. math::
  \nabla \cdot k \nabla u

with a variety of discretizations. Note also, for reasons that are one part historical 
and potentially not that valid, this also supports and implementation with an advective 
source, i.e.:

.. math::
  \nabla \cdot K k (\nabla u + b g z)

for gravitational terms in Richards equations.

The input spec for a diffusion operator consists of:

* `"discretization primary`" ``[string]`` See below for supported options.

 - `"fv: default`" the standard two-point flux finite volume discretization
 - `"nlfv: default`" the nonlinear finite volume method of ???
 - MFD methods, including:
  - `"mfd: default`"
  - `"mfd: monotone for hex`"
  - `"mfd: optimized for monotonicity`"
  - `"mfd: two-point flux approximation`"
  - `"mfd: optimized for sparsity`"
  - `"mfd: support operator`"

 Note that the most commonly used are `"fv: default`" for simple test
 problems (this method is not particularly accurate for distorted
 meshes), `"mfd: optimized for sparsity`" for most real problems on
 unstructured meshes, and `"mfd: optimized for monotonicity`" for
 orthogonal meshes with diagonal tensor/scalar coefficients.

* `"gravity`" ``[bool]`` **false** specifies if the gravitational flow term is included

* `"Newton correction`" ``[string]`` specifies a model for non-physical terms 
  that must be added to the matrix. These terms represent Jacobian and are needed 
  for the preconditioner. Available options are `"true Jacobian`" and `"approximate Jacobian`".
  The FV scheme accepts only the first options. The other schemes accept only the second option.

* `"scaled constraint equation`" ``[bool]`` **false** rescales flux continuity equations
  on mesh faces.  These equations are formed without the nonlinear
  coefficient. This option allows us to treat the case of zero nonlinear
  coefficient, which otherwise generates zero rows in the operator, which is
  then singular.  At moment this feature does not work with non-zero gravity
  term.

* `"constraint equation scaling cutoff`" ``[double]`` specifies the cutoff value for
  applying rescaling strategy described above.
*/

namespace Amanzi {
namespace Operators {

class BCs;

class PDE_DiffusionFactory {
 public:
  PDE_DiffusionFactory() {};
  PDE_DiffusionFactory(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);
  PDE_DiffusionFactory(Teuchos::ParameterList& oplist,
                       const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);

  // Lazy setup requires minimum input or none.
  // Multiple operator can be generated using a common setup
  void SetConstantTensorCoefficient(const WhetStone::Tensor& K);
  void SetVariableTensorCoefficient(const Teuchos::RCP<const std::vector<WhetStone::Tensor> >& K) { K_ = K; }

  void SetConstantScalarCoefficient(double k);
  void SetVariableScalarCoefficient(const Teuchos::RCP<const CompositeVector>& k,
                                    const Teuchos::RCP<const CompositeVector>& dkdu = Teuchos::null) {
    k_ = k;
    dkdu_ = dkdu;
  }

  void SetConstantGravitationalTerm(const AmanziGeometry::Point& g, double b = 1.0);
  void SetVariableGravitationalTerm(const AmanziGeometry::Point& g,
                                    const Teuchos::RCP<const CompositeVector>& b,
                                    const Teuchos::RCP<const CompositeVector>& dbdu = Teuchos::null) {
    g_ = g;
    b_ = b;
    dbdu_ = dbdu;
  }

  Teuchos::RCP<PDE_Diffusion> Create();

  // Backward compatibility
  // -- Diffusion-type PDEs with optional gravity.
  //    Decision is made based on data in the parameter list.
  Teuchos::RCP<PDE_Diffusion>
  Create(Teuchos::ParameterList& oplist,
         const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
         const Teuchos::RCP<BCs>& bc,
         double rho,
         const AmanziGeometry::Point& g);

  Teuchos::RCP<PDE_Diffusion>
  Create(Teuchos::ParameterList& oplist,
         const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
         const Teuchos::RCP<BCs>& bc,
         const Teuchos::RCP<const CompositeVector>& rho,
         const AmanziGeometry::Point& g);

  // -- Diffusion operators without gravity.
  Teuchos::RCP<PDE_Diffusion>
  Create(Teuchos::ParameterList& oplist,
         const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
         const Teuchos::RCP<BCs>& bc);

  Teuchos::RCP<PDE_Diffusion>
  Create(Teuchos::ParameterList& oplist,
         const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);
  
  Teuchos::RCP<PDE_Diffusion>
  Create(Teuchos::ParameterList& oplist,
         const Teuchos::RCP<Operator>& global_op);

  // Diffusion operators with gravity.
  Teuchos::RCP<PDE_DiffusionWithGravity>
  CreateWithGravity(Teuchos::ParameterList& oplist,
                    const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                    const Teuchos::RCP<BCs>& bc);
                    
  Teuchos::RCP<PDE_DiffusionWithGravity>
  CreateWithGravity(Teuchos::ParameterList& oplist,
                    const Teuchos::RCP<Operator>& global_op,
                    const Teuchos::RCP<BCs>& bc);

  Teuchos::RCP<PDE_DiffusionWithGravity>
  CreateWithGravity(Teuchos::ParameterList& oplist,
                    const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);
                    
  Teuchos::RCP<PDE_DiffusionWithGravity>
  CreateWithGravity(Teuchos::ParameterList& oplist,
                    const Teuchos::RCP<Operator>& global_op);

 private:
  Teuchos::ParameterList oplist_;
  const Teuchos::RCP<const AmanziMesh::Mesh> mesh_;

  // diffusivity
  WhetStone::Tensor const_K_;
  Teuchos::RCP<const std::vector<WhetStone::Tensor> > K_;

  double const_k_;
  Teuchos::RCP<const CompositeVector> k_, dkdu_;

  // gravity
  bool gravity_;
  AmanziGeometry::Point g_;
  double const_b_;
  Teuchos::RCP<const CompositeVector> b_, dbdu_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif
