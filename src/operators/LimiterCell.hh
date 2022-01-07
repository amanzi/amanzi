/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Implementation of different limiters uses a few common rules:
  1. Dirichlet boundary data are used to update limiter bounds.
  2. Limiters are modified optionally so the the stable time step
     of first-order scheme is reduce not more than twice. This
     step requires to specify a face-based flux field.
  3. At the moment, we require both the input field and boundary
     conditions to be defined at ghost positions.
*/

#ifndef AMANZI_LIMITER_CELL_HH_
#define AMANZI_LIMITER_CELL_HH_

#include <functional>
#include <vector>

#include "Epetra_IntVector.h"
#include "Epetra_MultiVector.h"
#include "Teuchos_RCP.hpp"

#include "CompositeVector.hh"
#include "DenseMatrix.hh"
#include "DenseVector.hh"
#include "DG_Modal.hh"
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
                    const std::vector<int>& bc_model, const std::vector<double>& bc_value) {
    AmanziMesh::Entity_ID_List ids(ncells_owned_);
    for (int c = 0; c < ncells_owned_; ++c) ids[c] = c;
    ApplyLimiter(ids, field, component, gradient, bc_model, bc_value); 
  }

  void ApplyLimiter(Teuchos::RCP<const Epetra_MultiVector> field, int component,
                    const Teuchos::RCP<CompositeVector>& gradient) {
    std::vector<int> bc_model;
    std::vector<double> bc_value;
    AmanziMesh::Entity_ID_List ids(ncells_owned_);
    for (int c = 0; c < ncells_owned_; ++c) ids[c] = c;
    ApplyLimiter(ids, field, component, gradient, bc_model, bc_value); 
  }

  // -- apply limiter in specified cells
  // -- apply limiter in specified cells
  void ApplyLimiter(const AmanziMesh::Entity_ID_List& ids,
                    Teuchos::RCP<const Epetra_MultiVector> field, int component,
                    const Teuchos::RCP<CompositeVector>& gradient,
                    const std::vector<int>& bc_model, const std::vector<double>& bc_value);

  // -- apply external limiter 
  void ApplyLimiter(Teuchos::RCP<Epetra_MultiVector> limiter);

  // limited dG solution 
  // -- apply limiter in spcified cells
  void ApplyLimiter(Teuchos::RCP<const Epetra_MultiVector> field, const WhetStone::DG_Modal& dg,
                    const std::vector<int>& bc_model, const std::vector<double>& bc_value) {
    AmanziMesh::Entity_ID_List ids(ncells_owned_);
    for (int c = 0; c < ncells_owned_; ++c) ids[c] = c;
    ApplyLimiter(ids, field, dg, bc_model, bc_value); 
  }

  void ApplyLimiter(const AmanziMesh::Entity_ID_List& ids,
                    Teuchos::RCP<const Epetra_MultiVector> field, const WhetStone::DG_Modal& dg,
                    const std::vector<int>& bc_model, const std::vector<double>& bc_value);

  // bounds for FV fields: if reset=true they are recalculated 
  Teuchos::RCP<CompositeVector> BoundsForCells(
      const Epetra_MultiVector& field, 
      const std::vector<int>& bc_model, const std::vector<double>& bc_value, int stencil);
  Teuchos::RCP<CompositeVector> BoundsForFaces(
      const Epetra_MultiVector& field,
      const std::vector<int>& bc_model, const std::vector<double>& bc_value, int stencil);
  Teuchos::RCP<CompositeVector> BoundsForEdges(
      const Epetra_MultiVector& field,
      const std::vector<int>& bc_model, const std::vector<double>& bc_value, int stencil);
  Teuchos::RCP<CompositeVector> BoundsForNodes(
      const Epetra_MultiVector& field,
      const std::vector<int>& bc_model, const std::vector<double>& bc_value, int stencil);

  // calculate value of a linear function at the given point p
  void getBounds(int c, int f, int stencil, double* umin, double* umax);
  double getValue(int c, const AmanziGeometry::Point& p);
  double getValue(const AmanziGeometry::Point& gradient, int c, const AmanziGeometry::Point& p);

  // access
  Teuchos::RCP<CompositeVector> gradient() { return gradient_; }
  Teuchos::RCP<Epetra_Vector> limiter() { return limiter_; }
  Teuchos::RCP<CompositeVector> bounds() { return bounds_; }
  int type() { return type_; }

  // modifiers
  void set_gradient(const Teuchos::RCP<CompositeVector>& gradient) { gradient_ = gradient; }
  void set_bounds(const Teuchos::RCP<CompositeVector>& bounds) { bounds_ = bounds; }
 
 private:
  // scalar limiters
  void LimiterScalar_(
      const AmanziMesh::Entity_ID_List& ids,
      const std::vector<int>& bc_model, const std::vector<double>& bc_value,
      Teuchos::RCP<Epetra_Vector> limiter, double (*)(double));

  void LimiterScalarDG_(
      const WhetStone::DG_Modal& dg, const AmanziMesh::Entity_ID_List& ids,
      const std::vector<int>& bc_model, const std::vector<double>& bc_value, double (*)(double));

  void LimiterTensorial_(
      const AmanziMesh::Entity_ID_List& ids,
      const std::vector<int>& bc_model, const std::vector<double>& bc_value);

  void LimiterKuzmin_(
      const AmanziMesh::Entity_ID_List& ids,
      const std::vector<int>& bc_model, const std::vector<double>& bc_value);

  // hierarchical limiters
  void LimiterHierarchicalDG_(
      const WhetStone::DG_Modal& dg, const AmanziMesh::Entity_ID_List& ids,
      const std::vector<int>& bc_model, const std::vector<double>& bc_value, double (*)(double));

  // supprting routines for limiters
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

  void LimiterExtensionTransportScalar_(Teuchos::RCP<Epetra_Vector> limiter);
  void LimiterExtensionTransportTensorial_();
  void LimiterExtensionTransportKuzmin_(
      const std::vector<double>& field_local_min, const std::vector<double>& field_local_max);

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  int dim;
  int ncells_owned_, nfaces_owned_, nedges_owned_, nnodes_owned_;
  int ncells_wghost_, nfaces_wghost_, nedges_wghost_, nnodes_wghost_;

  Teuchos::RCP<const Epetra_MultiVector> field_;
  Teuchos::RCP<CompositeVector> gradient_, bounds_;
  Teuchos::RCP<Epetra_Vector> limiter_;
  int component_;

  Teuchos::RCP<const Epetra_MultiVector> flux_;  // for limiters
  std::vector<std::vector<int> > upwind_cells_;  // fracture friendly 
  std::vector<std::vector<int> > downwind_cells_;

  int type_, stencil_id_;
  AmanziMesh::Entity_kind location_;
  bool limiter_correction_, external_bounds_;
  double cfl_;

  int limiter_points_;  // number of Gauss points on faces where limiting occurs
};

}  // namespace Operators
}  // namespace Amanzi

#endif
