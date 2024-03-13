/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

  A face-based field is defined by volume averaging of cell-centered data.
*/

#ifndef AMANZI_ARITHMETIC_AVERAGE_HH_
#define AMANZI_ARITHMETIC_AVERAGE_HH_

#include <string>
#include <vector>

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

// Amanzi
#include "CompositeVector.hh"
#include "Mesh.hh"
#include "Point.hh"

// Operators
#include "Upwind.hh"

namespace Amanzi {
namespace Operators {

class UpwindArithmeticAverage : public Upwind {
 public:
  UpwindArithmeticAverage(Teuchos::RCP<const AmanziMesh::Mesh> mesh) : Upwind(mesh){};
  ~UpwindArithmeticAverage(){};

  // main methods
  void Init(Teuchos::ParameterList& plist);

  void
  Compute(const CompositeVector& flux, const std::vector<int>& bc_model, CompositeVector& field);

 private:
  int method_, order_;
  double tolerance_;
};


/* ******************************************************************
* Public init method. It is not yet used.
****************************************************************** */
inline void
UpwindArithmeticAverage::Init(Teuchos::ParameterList& plist)
{
  method_ = OPERATOR_UPWIND_ARITHMETIC_AVERAGE;
  tolerance_ = plist.get<double>("tolerance", OPERATOR_UPWIND_RELATIVE_TOLERANCE);
  order_ = plist.get<int>("polynomial order", 1);
}


/* ******************************************************************
* Upwind field is placed in component "face".
* Upwinded field must be calculated on all faces of the owned cells.
****************************************************************** */
inline void
UpwindArithmeticAverage::Compute(const CompositeVector& flux,
                                 const std::vector<int>& bc_model,
                                 CompositeVector& field)
{
  AMANZI_ASSERT(field.HasComponent("cell"));
  AMANZI_ASSERT(field.HasComponent(face_comp_));

  field.ScatterMasterToGhosted("cell");

  Epetra_MultiVector& fld_cell = *field.ViewComponent("cell", true);
  Epetra_MultiVector& upw_face = *field.ViewComponent(face_comp_, true);

  int nfaces_wghost =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::ALL);

  int c1, c2;
  double kc1, kc2;
  for (int f = 0; f < nfaces_wghost; ++f) {
    auto cells = mesh_->getFaceCells(f);
    int ncells = cells.size();

    c1 = cells[0];
    kc1 = fld_cell[0][c1];

    if (ncells == 2) {
      c2 = cells[1];
      kc2 = fld_cell[0][c2];

      double v1 = mesh_->getCellVolume(c1);
      double v2 = mesh_->getCellVolume(c2);

      double tmp = v2 / (v1 + v2);
      upw_face[0][f] = kc1 * tmp + kc2 * (1.0 - tmp);

    } else {
      upw_face[0][f] = kc1;
    }
  }
}

} // namespace Operators
} // namespace Amanzi

#endif
