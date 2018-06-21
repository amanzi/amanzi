/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_RECONSTRUCTION_CELL_HH_
#define AMANZI_RECONSTRUCTION_CELL_HH_

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

class ReconstructionCell : public Reconstruction {  
 public:
  ReconstructionCell() {};
  ReconstructionCell(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh) : Reconstruction(mesh) {};
  ~ReconstructionCell() {};

  // save pointer to the already distributed field.
  void Init(Teuchos::RCP<const Epetra_MultiVector> field, Teuchos::ParameterList& plist, int component = 0);

  // global reconstruction
  // -- compute gradient and keep it internally
  void Compute();

  // -- identify inflow boundaries (optional)
  void InitLimiter(Teuchos::RCP<const Epetra_MultiVector> flux);
  // -- limit gradient using boundary data
  virtual void ApplyLimiter(const std::vector<int>& bc_model,
                            const std::vector<double>& bc_value) override;
  // -- apply external limiter 
  virtual void ApplyLimiter(Teuchos::RCP<Epetra_MultiVector> limiter) override;

  // -- calculate value of a reconstructed function at given point p
  virtual double getValue(int c, const AmanziGeometry::Point& p) override;
  double getValue(AmanziGeometry::Point& gradient, int c, const AmanziGeometry::Point& p);

  // local reconstruction
  // -- compute gradient in specified cells
  void ComputeGradient(const AmanziMesh::Entity_ID_List& ids,
                       std::vector<AmanziGeometry::Point>& gradient);

  // -- apply limiter in specified cells
  void ApplyLimiter(AmanziMesh::Entity_ID_List& ids,
                    std::vector<AmanziGeometry::Point>& gradient);

  // access
  Teuchos::RCP<CompositeVector> gradient() { return gradient_; }
  Teuchos::RCP<Epetra_Vector> limiter() { return limiter_; }
 
 private:
  void PopulateLeastSquareSystem_(AmanziGeometry::Point& centroid,
                                  double field_value,
                                  WhetStone::DenseMatrix& matrix,
                                  WhetStone::DenseVector& rhs);

  // internal limiters and supporting routines
  void LimiterBarthJespersen_(
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

  void LimiterExtensionTransportTensorial_(
      const std::vector<double>& field_local_min, const std::vector<double>& field_local_max);

  void LimiterExtensionTransportKuzmin_(
      const std::vector<double>& field_local_min, const std::vector<double>& field_local_max);

  void LimiterExtensionTransportBarthJespersen_(
      const std::vector<double>& field_local_min, const std::vector<double>& field_local_max,
      Teuchos::RCP<Epetra_Vector> limiter);

  // On intersecting manifolds, we extract neighboors living in the same manifold
  // using a smoothness criterion.
  void CellFaceAdjCellsNonManifold_(AmanziMesh::Entity_ID c,
                                    AmanziMesh::Parallel_type ptype,
                                    std::vector<AmanziMesh::Entity_ID>& cells) const;
 private:
  int dim;
  int ncells_owned, nfaces_owned, nnodes_owned;
  int ncells_wghost, nfaces_wghost, nnodes_wghost;
  int cell_max_nodes;

  Teuchos::RCP<CompositeVector> gradient_;
  Teuchos::RCP<Epetra_Vector> limiter_;

 private: 
  Teuchos::RCP<const Epetra_MultiVector> flux_;  // for limiters
  std::vector<std::vector<int> > upwind_cells_;  // fracture friendly 
  std::vector<std::vector<int> > downwind_cells_;

  double bc_scaling_;
  int limiter_id_, poly_order_;
  bool limiter_correction_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif
