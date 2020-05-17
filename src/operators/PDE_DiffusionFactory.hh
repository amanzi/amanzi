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
#include "PDE_DiffusionFV.hh"
#include "PDE_DiffusionFVwithGravity.hh"
// #include "PDE_DiffusionFracturedMatrix.hh"
#include "PDE_DiffusionMFD.hh"
// #include "PDE_DiffusionMFDwithGravity.hh"
// #include "PDE_DiffusionNLFV.hh"
// #include "PDE_DiffusionNLFVwithBndFaces.hh"
// #include "PDE_DiffusionNLFVwithBndFacesGravity.hh"
// #include "PDE_DiffusionNLFVwithGravity.hh"

/*!

``PDE_Diffusion`` forms local ``Op`` s and global ``Operator`` s for elliptic equations:

.. math::
  \nabla \cdot k \nabla u

with a variety of discretizations. Note also, for reasons that are one part historical 
and potentially not that valid, this also supports and implementation with an advective 
source, i.e.:

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

struct PDE_DiffusionFactory {
  // A variety of users with a variety of natural interfaces.

  // Diffusion-type PDEs with optional gravity.
  // Decision is made based on data in the parameter list.
  Teuchos::RCP<PDE_Diffusion>
  Create(Teuchos::ParameterList& oplist,
         const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
         const Teuchos::RCP<BCs>& bc,
         double rho,
         const AmanziGeometry::Point& g)
  {
    if (oplist.get<bool>("gravity", false)) {
      auto op = CreateWithGravity_(oplist, mesh, bc);
      op->SetDensity(rho);
      op->SetGravity(g);
      return op;
      
    } else {
      return CreateWithoutGravity_(oplist, mesh, bc);
    }
  }

  Teuchos::RCP<PDE_Diffusion>
  Create(Teuchos::ParameterList& oplist,
         const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
         const Teuchos::RCP<BCs>& bc,
         const Teuchos::RCP<const CompositeVector>& rho,
         const AmanziGeometry::Point& g)
  {
    if (oplist.get<bool>("gravity", false)) {
      auto op = CreateWithGravity_(oplist, mesh, bc);
      op->SetDensity(rho);
      op->SetGravity(g);
      return op;
    } else {
      return CreateWithoutGravity_(oplist, mesh, bc);
    }
  }

  // Diffusion operators without gravity.
  Teuchos::RCP<PDE_Diffusion>
  Create(Teuchos::ParameterList& oplist,
         const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
         const Teuchos::RCP<BCs>& bc=Teuchos::null)
  {
    return CreateWithoutGravity_(oplist, mesh, bc);
  }
  
  Teuchos::RCP<PDE_Diffusion>
  Create(Teuchos::ParameterList& oplist,
         const Teuchos::RCP<Operator>& global_op,
         const Teuchos::RCP<BCs>& bc=Teuchos::null)
  {
    return CreateWithoutGravity_(oplist, global_op, bc);
  }
         
  // Diffusion operators with gravity.
  Teuchos::RCP<PDE_Diffusion>
  CreateWithGravity(Teuchos::ParameterList& oplist,
                    const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                    const Teuchos::RCP<BCs>& bc=Teuchos::null)
  {
    return CreateWithGravity_(oplist, mesh, bc);
  }
                    
  Teuchos::RCP<PDE_Diffusion>
  CreateWithGravity(Teuchos::ParameterList& oplist,
                    const Teuchos::RCP<Operator>& global_op,
                    const Teuchos::RCP<BCs>& bc=Teuchos::null)
  {
    return CreateWithGravity_(oplist, global_op, bc);
  }

 private:

  template<class Second_ptr_type>
  Teuchos::RCP<PDE_Diffusion>
  CreateWithoutGravity_(Teuchos::ParameterList& oplist,
          const Second_ptr_type& mesh_or_global_op,
          const Teuchos::RCP<BCs>& bc);

  template<class Second_ptr_type>
  Teuchos::RCP<PDE_Diffusion>
  CreateWithGravity_(Teuchos::ParameterList& oplist,
                     const Second_ptr_type& mesh_or_global_op,
                     const Teuchos::RCP<BCs>& bc);
};  


//
// Implementations
//
template<class Second_ptr_type>
Teuchos::RCP<PDE_Diffusion>
PDE_DiffusionFactory::CreateWithoutGravity_(
    Teuchos::ParameterList& oplist,
    const Second_ptr_type& mesh_or_global_op, 
    const Teuchos::RCP<BCs>& bc)
{
  std::string name = oplist.get<std::string>("discretization primary");
  bool fractured_matrix = oplist.isParameter("fracture");

  Teuchos::RCP<PDE_Diffusion> op;
  if (name == "fv: default") {
    op = Teuchos::rcp(new PDE_DiffusionFV(oplist, mesh_or_global_op));
  // } else if (name == "nlfv: default") {
  //   op = Teuchos::rcp(new PDE_DiffusionNLFV(oplist, mesh_or_global_op));
  // } else if (name == "nlfv: bnd_faces") {
  //   op = Teuchos::rcp(new PDE_DiffusionNLFVwithBndFaces(oplist,
  //           mesh_or_global_op));
  // } else if (fractured_matrix) {
  //   op = Teuchos::rcp(new PDE_DiffusionFracturedMatrix(oplist,
  //           mesh_or_global_op));
  } else {
    op = Teuchos::rcp(new PDE_DiffusionMFD(oplist, mesh_or_global_op));
  }
  op->Init();
  if (bc != Teuchos::null) op->SetBCs(bc, bc);
  return op;
}

template<class Second_ptr_type>
Teuchos::RCP<PDE_Diffusion>
PDE_DiffusionFactory::CreateWithGravity_(
    Teuchos::ParameterList& oplist,
    const Second_ptr_type& mesh_or_global_op, 
    const Teuchos::RCP<BCs>& bc)
{
  std::string name = oplist.get<std::string>("discretization primary");
  bool fractured_matrix = oplist.isParameter("fracture");

  Teuchos::RCP<PDE_Diffusion> op_g;
  if (name == "fv: default") {
    op_g = Teuchos::rcp(new PDE_DiffusionFVwithGravity(oplist, mesh_or_global_op));
  // } else if (name == "nlfv: default") {
  //   op_g = Teuchos::rcp(new PDE_DiffusionNLFVwithGravity(oplist, mesh_or_global_op));
  // } else if (name == "nlfv: bnd_faces") {
  //   op_g = Teuchos::rcp(new PDE_DiffusionNLFVwithBndFacesGravity(oplist, mesh_or_global_op));
  // } else if (fractured_matrix) {
  //   op_g = Teuchos::rcp(new PDE_DiffusionFracturedMatrix(oplist, mesh_or_global_op));
  // } else {
  //   op_g = Teuchos::rcp(new PDE_DiffusionMFDwithGravity(oplist, mesh_or_global_op));
  }
  op_g->Init();
  if (bc != Teuchos::null) op_g->SetBCs(bc, bc);
  return op_g;
}





}  // namespace Operators
}  // namespace Amanzi

#endif
