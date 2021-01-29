/*
  Operators 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
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
* "dirichlet_faces" - boundary_face values of the field evaluated on
*          the boundary (either as a function of the internal cell
*          or as a function of the Dirichlet data)
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

template<class Model>
class Upwind {
 public:
  Upwind() {};
  Upwind(Teuchos::RCP<const AmanziMesh::Mesh> mesh,
         Teuchos::RCP<const Model> model) :
      mesh_(mesh),
      model_(model),
      face_comp_("face") {};
  ~Upwind() {};

  // main methods
  // -- initialization of control parameters
  virtual void Init(Teuchos::ParameterList& plist) = 0;

  // -- upwind of a given cell-centered field on mesh faces
  // -- not all input parameters are use by some algorithms
  virtual void Compute(const CompositeVector& flux, const CompositeVector& solution,
                       const std::vector<int>& bc_model, CompositeVector& field) = 0;

  // -- returns combined map for the original and upwinded fields.
  // -- Currently, composite vector cannot be extended on a fly. 
  virtual Teuchos::RCP<CompositeVectorSpace> Map() {
    Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
    cvs->SetMesh(mesh_)->SetGhosted(true)
        ->AddComponent("cell", AmanziMesh::CELL, 1)
        ->AddComponent("dirichlet_faces", AmanziMesh::BOUNDARY_FACE, 1)
        ->AddComponent("face", AmanziMesh::FACE, 1);
    return cvs;
  }

  void PopulateDirichletFaces(const std::vector<int>& bc_model,
                              CompositeVector& field);

  // modifiers
  void set_face_comp(const std::string& name) { face_comp_ = name; }

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<const Model> model_;

  // component name where to write the upwinded field.
  std::string face_comp_;
};


/* ******************************************************************
* Support function: copy data from cells to dirichlet faces
****************************************************************** */
template<class Model>
void Upwind<Model>::PopulateDirichletFaces(const std::vector<int>& bc_model,
                                           CompositeVector& field)
{
  Epetra_MultiVector& field_df = *field.ViewComponent("dirichlet_faces", true);
  Epetra_MultiVector& field_c = *field.ViewComponent("cell", true);

  const Epetra_Map& ext_face_map = mesh_->exterior_face_map(true);
  const Epetra_Map& face_map = mesh_->face_map(true);

  for (int f = 0; f != face_map.NumMyElements(); ++f) {
    if (bc_model[f] == Operators::OPERATOR_BC_DIRICHLET) {
      int bf = ext_face_map.LID(face_map.GID(f));

      AmanziMesh::Entity_ID_List cells;
      mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);

      field_df[0][bf] = field_c[0][cells[0]];
    }
  }
}

}  // namespace Operators
}  // namespace Amanzi

#endif

