/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_LIMITER_CELL_HH_
#define AMANZI_LIMITER_CELL_HH_

#include <vector>

#include "Epetra_IntVector.h"
#include "Epetra_MultiVector.h"
#include "Teuchos_RCP.hpp"

#include "CompositeVector.hh"
#include "DenseMatrix.hh"
#include "DenseVector.hh"
#include "Mesh.hh"
#include "Point.hh"

#include "Reconstruction.hh"


namespace Amanzi {
namespace Operators {

class LimiterCell {  
 public:
  LimiterCell(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh);
  ~LimiterCell() {};

  // limited gradient
  // -- identify inflow boundaries (optional)
  void Init(Teuchos::ParameterList& plist,
            Teuchos::RCP<const Epetra_MultiVector> flux = Teuchos::null);

  // -- limit gradient using boundary data
  void ApplyLimiter(Teuchos::RCP<const Epetra_MultiVector> field, int component,
                    const Teuchos::RCP<CompositeVector>& gradient,
                    const std::vector<int>& bc_model, const std::vector<double>& bc_value);

  // -- apply external limiter 
  void ApplyLimiter(Teuchos::RCP<Epetra_MultiVector> limiter);

  // -- apply limiter in specified cells
  void ApplyLimiter(AmanziMesh::Entity_ID_List& ids,
                    Teuchos::RCP<const Epetra_MultiVector> field, int component,
                    std::vector<AmanziGeometry::Point>& gradient);

  // bounds
  void BoundsFaceToCells();
  void BoundsCellToClosestCells();
  void BoundsCellToAllCells();
  void BoundsCellOnBoundary(const std::vector<int>& bc_model, const std::vector<double>& bc_value);

  // calculate value of a linear function at the given point p
  double getValue(int c, const AmanziGeometry::Point& p);
  double getValue(const AmanziGeometry::Point& gradient, int c, const AmanziGeometry::Point& p);

  // access
  Teuchos::RCP<CompositeVector> gradient() { return gradient_; }
  Teuchos::RCP<Epetra_Vector> limiter() { return limiter_; }
  Teuchos::RCP<CompositeVector> bounds() { return bounds_; }

  // modifiers
  void set_gradient(const Teuchos::RCP<CompositeVector>& gradient) { gradient_ = gradient; }
  void set_bounds(const Teuchos::RCP<CompositeVector>& bounds) { bounds_ = bounds; }
 
 private:
  // internal limiters and supporting routines
  void LimiterBarthJespersenFace_(
      const std::vector<int>& bc_model, const std::vector<double>& bc_value,
      Teuchos::RCP<Epetra_Vector> limiter);

  void LimiterBarthJespersenCell_(
      const std::vector<int>& bc_model, const std::vector<double>& bc_value,
      Teuchos::RCP<Epetra_Vector> limiter);

  void LimiterTensorial_(
      const std::vector<int>& bc_model, const std::vector<double>& bc_value);

  void LimiterKuzmin_(
      const std::vector<int>& bc_model, const std::vector<double>& bc_value);

  void LimiterKuzminSet_(AmanziMesh::Entity_ID_List& ids,
                         std::vector<AmanziGeometry::Point>& gradient);

  void LimiterKuzminCell_(int cell,
                          AmanziGeometry::Point& gradient_c,
                          const std::vector<double>& field_node_min_c,
                          const std::vector<double>& field_node_max_c);

  void CalculateDescentDirection_(std::vector<AmanziGeometry::Point>& normals,
                                  AmanziGeometry::Point& normal_new,
                                  double& L22normal_new, 
                                  AmanziGeometry::Point& direction);

  void ApplyDirectionalLimiter_(AmanziGeometry::Point& normal, 
                                AmanziGeometry::Point& p,
                                AmanziGeometry::Point& direction, 
                                AmanziGeometry::Point& gradient);

  void IdentifyUpwindCells_();

  void LimiterExtensionTransportTensorial_();
  void LimiterExtensionTransportBarthJespersen_(Teuchos::RCP<Epetra_Vector> limiter);
  void LimiterExtensionTransportKuzmin_(
      const std::vector<double>& field_local_min, const std::vector<double>& field_local_max);

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  int dim;
  int ncells_owned, nfaces_owned, nnodes_owned;
  int ncells_wghost, nfaces_wghost, nnodes_wghost;
  int cell_max_nodes;

  Teuchos::RCP<const Epetra_MultiVector> field_;
  Teuchos::RCP<CompositeVector> gradient_, bounds_;
  Teuchos::RCP<Epetra_Vector> limiter_;
  int component_;

  Teuchos::RCP<const Epetra_MultiVector> flux_;  // for limiters
  std::vector<std::vector<int> > upwind_cells_;  // fracture friendly 
  std::vector<std::vector<int> > downwind_cells_;

  double bc_scaling_;
  int limiter_id_, stencil_id_;
  bool limiter_correction_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif
