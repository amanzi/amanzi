/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

// PDE_DiffusionFactory constructs objects which implement the interface for a
// PDE_Diffusion.


#ifndef AMANZI_OPERATOR_PDE_DIFFUSION_FACTORY_HH_
#define AMANZI_OPERATOR_PDE_DIFFUSION_FACTORY_HH_

#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "PDE_Diffusion.hh"
#include "PDE_DiffusionWithGravity.hh"

/*!

``PDE_Diffusion`` forms local ``Op`` s and global ``Operator`` s for elliptic
equations:

.. math::
  \nabla \cdot k \nabla u

with a variety of discretizations. Note also, for reasons that are one part
historical and potentially not that valid, this also supports and implementation
with an advective source, i.e.:

.. math::
  \nabla \cdot k (\nabla u + \hat{z})

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

* `"gravity`" ``[bool]`` **false** specifies if the gravitational flow term is
included

* `"Newton correction`" ``[string]`` specifies a model for non-physical terms
  that must be added to the matrix. These terms represent Jacobian and are
needed for the preconditioner. Available options are `"true Jacobian`" and
`"approximate Jacobian`". The FV scheme accepts only the first options. The
other schemes accept only the second option.

* `"scaled constraint equation`" ``[bool]`` **false** rescales flux continuity
equations on mesh faces.  These equations are formed without the nonlinear
  coefficient. This option allows us to treat the case of zero nonlinear
  coefficient, which otherwise generates zero rows in the operator, which is
  then singular.  At moment this feature does not work with non-zero gravity
  term.

* `"constraint equation scaling cutoff`" ``[double]`` specifies the cutoff value
for applying rescaling strategy described above.
*/

namespace Amanzi {
namespace Operators {

class BCs;

struct PDE_DiffusionFactory {
  // Diffusion-type PDEs with optional gravity.
  // Decision is made based on data in the parameter list.
  Teuchos::RCP<PDE_Diffusion>
  Create(Teuchos::ParameterList& oplist,
         const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
         const Teuchos::RCP<BCs>& bc, double rho,
         const AmanziGeometry::Point& g);

  Teuchos::RCP<PDE_Diffusion>
  Create(Teuchos::ParameterList& oplist,
         const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
         const Teuchos::RCP<BCs>& bc,
         const Teuchos::RCP<const CompositeVector>& rho,
         const AmanziGeometry::Point& g);

  // Diffusion operators without gravity.
  Teuchos::RCP<PDE_Diffusion>
  Create(Teuchos::ParameterList& oplist,
         const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
         const Teuchos::RCP<BCs>& bc);

  Teuchos::RCP<PDE_Diffusion>
  Create(Teuchos::ParameterList& oplist,
         const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);

  Teuchos::RCP<PDE_Diffusion> Create(Teuchos::ParameterList& oplist,
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
  inline void SetCellSchema_(Teuchos::ParameterList& oplist);
  inline void SetCellFaceSchema_(Teuchos::ParameterList& oplist);
};

} // namespace Operators
} // namespace Amanzi

#endif
