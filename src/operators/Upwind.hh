/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

*/

#ifndef AMANZI_UPWIND_HH_
#define AMANZI_UPWIND_HH_

#include <string>
#include <vector>

#include "Epetra_IntVector.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "CompositeVector.hh"
#include "CompositeVectorSpace.hh"
#include "Mesh.hh"
#include "VerboseObject.hh"

#include "OperatorDefs.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* The base class for all upwind methods.
* Currently only one use case is implemented. The input field
* contains (up to) three components:
*
* "cell" - cell-centered values of the field.
* "boundary_face" - field values evaluated on the boundary either
           from internal cell or as a function of Dirichlet data.
* "grad" - (optional) eatimate of the gradient of the input field.
*          It is used in the second-order upwind schemes.
*
* The output field contains one or two components depending on a
* numerical scheme:
*
*  "face" - upwinded value of the component "cell" on mesh faces.
*  "twin" - second upwinded value on mesh faces. It may be useful
*           when the input field has physical (not numerical)
*           discontinuous, e.g. different permeability curves.
*
* Amanzi combines the input and output field in one variable.
****************************************************************** */

class Upwind {
 public:
  Upwind(){};
  Upwind(Teuchos::RCP<const AmanziMesh::Mesh> mesh) : mesh_(mesh), face_comp_("face"){};
  virtual ~Upwind(){};

  // main methods
  // -- initialization of control parameters
  virtual void Init(Teuchos::ParameterList& plist) = 0;

  // -- upwind of a given cell-centered field on mesh faces
  // -- not all input parameters are use by some algorithms
  virtual void Compute(const CompositeVector& flux,
                       const CompositeVector& solution,
                       const std::vector<int>& bc_model,
                       CompositeVector& field) = 0;

  // -- returns combined map for the original and upwinded fields.
  // -- Currently, composite vector cannot be extended on a fly.
  virtual Teuchos::RCP<CompositeVectorSpace> Map()
  {
    Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
    cvs->SetMesh(mesh_)
      ->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1)
      ->AddComponent("face", AmanziMesh::Entity_kind::FACE, 1);
    return cvs;
  }

  // modifiers
  void set_face_comp(const std::string& name) { face_comp_ = name; }

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  std::string face_comp_; // component where to write the upwinded field.
};

} // namespace Operators
} // namespace Amanzi

#endif
