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

#ifndef AMANZI_LIMITER_CELL_DG_HH_
#define AMANZI_LIMITER_CELL_DG_HH_

#include "LimiterCell.hh"
#include "DG_Modal.hh"

namespace Amanzi {
namespace Operators {

class LimiterCellDG : public LimiterCell {
 public:
  LimiterCellDG(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh);
  ~LimiterCellDG() {};

  void ApplyLimiterDG(const AmanziMesh::Entity_ID_View& ids,
                      Teuchos::RCP<const Epetra_MultiVector> field,
                      const WhetStone::DG_Modal& dg,
                      const std::vector<int>& bc_model,
                      const std::vector<double>& bc_value);

  void ApplyLimiterDG(Teuchos::RCP<const Epetra_MultiVector> field,
                      const WhetStone::DG_Modal& dg,
                      const std::vector<int>& bc_model,
                      const std::vector<double>& bc_value)
  {
    AmanziMesh::Entity_ID_View ids("ids", ncells_owned_);
    for (int c = 0; c < ncells_owned_; ++c) ids[c] = c;
    ApplyLimiterDG(ids, field, dg, bc_model, bc_value);
  }

 private:
  void LimiterScalarDG_(const WhetStone::DG_Modal& dg,
                        const AmanziMesh::Entity_ID_View& ids,
                        const std::vector<int>& bc_model,
                        const std::vector<double>& bc_value,
                        double (*)(double));

  void LimiterHierarchicalDG_(const WhetStone::DG_Modal& dg,
                              const AmanziMesh::Entity_ID_View& ids,
                              const std::vector<int>& bc_model,
                              const std::vector<double>& bc_value,
                              double (*)(double));
};

} // namespace Operators
} // namespace Amanzi

#endif
