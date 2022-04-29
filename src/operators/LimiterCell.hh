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

#include "Epetra_MultiVector.h"
#include "Teuchos_RCP.hpp"

#include "BCs.hh"
#include "CompositeVector.hh"
#include "Mesh.hh"
#include "Point.hh"

#include "Reconstruction.hh"

namespace Amanzi {
namespace Operators {

class LimiterCell {  
 public:
  LimiterCell(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh);
  ~LimiterCell() {};

  // limiting reconstruction data (gradient or full solution)
  // -- identify inflow boundaries (optional)
  void Init(Teuchos::ParameterList& plist,
            Teuchos::RCP<const Epetra_MultiVector> flux = Teuchos::null);

  // -- limit reconstructed data (typically gradeient) using neighboors
  //    and boundary data. Actual work is done by the 4th function.
  void ApplyLimiter(Teuchos::RCP<const Epetra_MultiVector> field, int component,
                    const Teuchos::RCP<Reconstruction>& lifting,
                    const Teuchos::RCP<const BCs>& bc);

  void ApplyLimiter(Teuchos::RCP<const Epetra_MultiVector> field, int component,
                    const Teuchos::RCP<Reconstruction>& lifting);

  void ApplyLimiter(Teuchos::RCP<const Epetra_MultiVector> field, int component,
                    const Teuchos::RCP<Reconstruction>& lifting,
                    const std::vector<int>& bc_model, const std::vector<double>& bc_value);

  void ApplyLimiter(const AmanziMesh::Entity_ID_List& ids,
                    Teuchos::RCP<const Epetra_MultiVector> field, int component,
                    const Teuchos::RCP<Reconstruction>& lifting,
                    const std::vector<int>& bc_model, const std::vector<double>& bc_value);

  // -- apply external limiter 
  void ApplyLimiter(Teuchos::RCP<Epetra_MultiVector> limiter);

  // bounds for cell-centered fields 
  Teuchos::RCP<CompositeVector> BoundsForCells(
      const Epetra_MultiVector& field, 
      const std::vector<int>& bc_model, const std::vector<double>& bc_value, int stencil) const;
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
  Teuchos::RCP<Epetra_Vector> limiter() { return limiter_; }
  Teuchos::RCP<CompositeVector> get_bounds() { return bounds_; }
  Teuchos::RCP<const CompositeVector> get_bounds() const { return bounds_; }

  int get_type() { return type_; }
  bool get_external_bounds() const { return external_bounds_; }

  // modifiers 
  void set_bounds(const Teuchos::RCP<CompositeVector>& bounds) { bounds_ = bounds; }
  void set_controls(const Teuchos::RCP<std::vector<std::vector<AmanziGeometry::Point> > >& controls) { controls_ = controls; }
 
 private:
  // scalar limiters
  void LimiterScalar_(
      const AmanziMesh::Entity_ID_List& ids,
      const std::vector<int>& bc_model, const std::vector<double>& bc_value,
      Teuchos::RCP<Epetra_Vector> limiter, double (*)(double));

  void LimiterTensorial_(
      const AmanziMesh::Entity_ID_List& ids,
      const std::vector<int>& bc_model, const std::vector<double>& bc_value);

  void LimiterKuzmin_(
      const AmanziMesh::Entity_ID_List& ids,
      const std::vector<int>& bc_model, const std::vector<double>& bc_value);

  // supporting routines for limiters
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

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  int dim;
  int ncells_owned_, nfaces_owned_, nedges_owned_, nnodes_owned_;
  int ncells_wghost_, nfaces_wghost_, nedges_wghost_, nnodes_wghost_;

  Teuchos::RCP<Reconstruction> lifting_;
  Teuchos::RCP<Epetra_MultiVector> data_;
  Teuchos::RCP<Epetra_Vector> limiter_;

  Teuchos::RCP<const Epetra_MultiVector> field_;
  Teuchos::RCP<CompositeVector> bounds_;
  int component_;

  Teuchos::RCP<const Epetra_MultiVector> flux_;  // for limiters
  std::vector<std::vector<int> > upwind_cells_;  // fracture friendly 
  std::vector<std::vector<int> > downwind_cells_;

  int type_, stencil_id_, location_;
  bool limiter_correction_, external_bounds_;

  int limiter_points_;  // number of Gauss points on faces where limiting occurs
  bool external_controls_;
  Teuchos::RCP<std::vector<std::vector<AmanziGeometry::Point> > > controls_;

  double cfl_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif
