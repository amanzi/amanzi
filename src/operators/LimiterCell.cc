/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

  Implementation of different limiters uses a few common rules:
  1. Dirichlet boundary data are used to update limiter bounds.
  2. Limiters are modified optionally so the the stable timestep
     of first-order scheme is reduce not more than twice. This
     step requires to specify a face-based flux field.
  3. At the moment, we require the input field and boundary data
     to have valid values in ghost positions. Exception, is the
     limiter for DG fields.
*/

#include <algorithm>
#include <vector>

// TPLs
#include "Epetra_Vector.h"

// Amanzi
#include "DenseMatrix.hh"
#include "DenseVector.hh"
#include "exceptions.hh"
#include "Quadrature1D.hh"
#include "WhetStoneMeshUtils.hh"

// Amanzi::Operators
#include "OperatorDefs.hh"
#include "OperatorUtils.hh"
#include "LimiterCell.hh"

namespace Amanzi {
namespace Operators {

/* *******************************************************************
* Work in progress: initialization like in Transport PK.
******************************************************************* */
LimiterCell::LimiterCell(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh)
  : mesh_(mesh),
    lifting_(Teuchos::null),
    data_(Teuchos::null),
    component_(0),
    flux_(Teuchos::null),
    external_bounds_(false),
    external_controls_(false),
    cfl_(1.0)
{
  dim = mesh_->getSpaceDimension();

  ncells_owned_ =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  nfaces_owned_ =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
  nnodes_owned_ =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_kind::OWNED);

  ncells_wghost_ =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::ALL);
  nfaces_wghost_ =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::ALL);
  nnodes_wghost_ =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_kind::ALL);

  if (mesh_->hasEdges()) {
    nedges_owned_ =
      mesh_->getNumEntities(AmanziMesh::Entity_kind::EDGE, AmanziMesh::Parallel_kind::OWNED);
    nedges_wghost_ =
      mesh_->getNumEntities(AmanziMesh::Entity_kind::EDGE, AmanziMesh::Parallel_kind::ALL);
  }
}


/* *******************************************************************
* Work in progress: initialization like in Transport PK.
******************************************************************* */
void
LimiterCell::Init(Teuchos::ParameterList& plist, const Teuchos::RCP<const Epetra_MultiVector>& flux)
{
  if (flux != Teuchos::null) SetFlux(flux);

  // process parameters for limiters
  std::string stencil;
  std::string name = plist.get<std::string>("limiter", "Barth-Jespersen");
  type_ = 0;
  if (name == "Barth-Jespersen") {
    type_ = OPERATOR_LIMITER_BARTH_JESPERSEN;
    stencil = plist.get<std::string>("limiter stencil", "face to cells");
  } else if (name == "Barth-Jespersen dg") {
    type_ = OPERATOR_LIMITER_BARTH_JESPERSEN_DG;
    stencil = plist.get<std::string>("limiter stencil", "face to cells");
  } else if (name == "Barth-Jespersen dg hierarchical") {
    type_ = OPERATOR_LIMITER_BARTH_JESPERSEN_DG_HIERARCHICAL;
    stencil = plist.get<std::string>("limiter stencil", "face to cells");
  } else if (name == "Michalak-Gooch") {
    type_ = OPERATOR_LIMITER_MICHALAK_GOOCH;
    stencil = plist.get<std::string>("limiter stencil", "face to cells");
  } else if (name == "Michalak-Gooch dg") {
    type_ = OPERATOR_LIMITER_MICHALAK_GOOCH_DG;
    stencil = plist.get<std::string>("limiter stencil", "face to cells");
  } else if (name == "tensorial") {
    type_ = OPERATOR_LIMITER_TENSORIAL;
    stencil = plist.get<std::string>("limiter stencil", "face to cells");
  } else if (name == "Kuzmin") {
    type_ = OPERATOR_LIMITER_KUZMIN;
    stencil = plist.get<std::string>("limiter stencil", "node to cells");
  }

  if (stencil == "node to cells") {
    stencil_id_ = OPERATOR_LIMITER_STENCIL_N2C;
    location_ = AmanziMesh::Entity_kind::NODE;
  } else if (stencil == "face to cells") {
    stencil_id_ = OPERATOR_LIMITER_STENCIL_F2C;
    location_ = AmanziMesh::Entity_kind::FACE;
  } else if (stencil == "edge to cells") {
    stencil_id_ = OPERATOR_LIMITER_STENCIL_E2C;
    location_ = AmanziMesh::Entity_kind::EDGE;
  } else if (stencil == "cell to closest cells") {
    stencil_id_ = OPERATOR_LIMITER_STENCIL_C2C_CLOSEST;
    location_ = AmanziMesh::Entity_kind::FACE;
  } else if (stencil == "cell to all cells") {
    stencil_id_ = OPERATOR_LIMITER_STENCIL_C2C_ALL;
    location_ = AmanziMesh::Entity_kind::FACE;
  } else {
    Errors::Message msg("unknown stencil for limiter");
    Exceptions::amanzi_throw(msg);
  }

  if (plist.isParameter("limiter location")) {
    name = plist.get<std::string>("limiter location");
    int tmp = AmanziMesh::createEntityKind(name);
    AMANZI_ASSERT((location_ == tmp) || (location_ == AmanziMesh::Entity_kind::FACE &&
                                         tmp == AmanziMesh::Entity_kind::NODE));
    location_ = tmp;
  }

  cfl_ = plist.get<double>("limiter cfl", 1.0);
  external_bounds_ = plist.get<bool>("use external bounds", false);
  external_controls_ = plist.get<bool>("use external controls", false);
  limiter_points_ = plist.get<int>("limiter points", 1);
  limiter_correction_ = plist.get<bool>("limiter extension for transport", false);
}


void
LimiterCell::SetFlux(const Teuchos::RCP<const Epetra_MultiVector>& flux)
{
  flux_ = flux;
  IdentifyUpwindCells_();
}


/* ******************************************************************
* Apply an internal limiter using BC object
****************************************************************** */
void
LimiterCell::ApplyLimiter(const Teuchos::RCP<const Epetra_MultiVector>& field,
                          int component,
                          const Teuchos::RCP<Reconstruction>& lifting,
                          const Teuchos::RCP<const BCs>& bc)
{
  const auto& bc_model = bc->bc_model();
  const auto& bc_value = bc->bc_value();

  AmanziMesh::Entity_ID_View ids("ids", ncells_owned_);
  for (int c = 0; c < ncells_owned_; ++c) ids[c] = c;

  ApplyLimiter(ids, field, component, lifting, bc_model, bc_value);
}


void
LimiterCell::ApplyLimiter(const Teuchos::RCP<const Epetra_MultiVector>& field,
                          int component,
                          const Teuchos::RCP<Reconstruction>& lifting,
                          const std::vector<int>& bc_model,
                          const std::vector<double>& bc_value)
{
  AmanziMesh::Entity_ID_View ids("ids", ncells_owned_);
  for (int c = 0; c < ncells_owned_; ++c) ids[c] = c;
  ApplyLimiter(ids, field, component, lifting, bc_model, bc_value);
}


void
LimiterCell::ApplyLimiter(const Teuchos::RCP<const Epetra_MultiVector>& field,
                          int component,
                          const Teuchos::RCP<Reconstruction>& lifting)
{
  std::vector<int> bc_model;
  std::vector<double> bc_value;
  AmanziMesh::Entity_ID_View ids("ids", ncells_owned_);
  for (int c = 0; c < ncells_owned_; ++c) ids[c] = c;
  ApplyLimiter(ids, field, component, lifting, bc_model, bc_value);
}


/* ******************************************************************
* Apply an internal limiter.
****************************************************************** */
void
LimiterCell::ApplyLimiter(const AmanziMesh::cEntity_ID_View& ids,
                          const Teuchos::RCP<const Epetra_MultiVector>& field,
                          int component,
                          const Teuchos::RCP<Reconstruction>& lifting,
                          const std::vector<int>& bc_model,
                          const std::vector<double>& bc_value)
{
  field_ = field;
  lifting_ = lifting;
  component_ = component;

  if (external_bounds_ && bounds_ == Teuchos::null) {
    Errors::Message msg("External bounds for limiters are requested but not provided");
    Exceptions::amanzi_throw(msg);
  }

  if (external_controls_ && controls_ == Teuchos::null) {
    Errors::Message msg("External control points for limiters are requested but not provided");
    Exceptions::amanzi_throw(msg);
  }

  // low level access is need for frequntly used functions
  data_ = lifting_->data()->ViewComponent("cell");

  if (limiter_ == Teuchos::null)
    limiter_ = Teuchos::rcp(new Epetra_Vector(mesh_->getMap(AmanziMesh::Entity_kind::CELL, false)));
  limiter_->PutScalar(1.0);

  if (type_ == OPERATOR_LIMITER_BARTH_JESPERSEN) {
    LimiterScalar_(ids, bc_model, bc_value, limiter_, [](double x) { return std::min(1.0, x); });
    ApplyLimiter(limiter_);
  } else if (type_ == OPERATOR_LIMITER_MICHALAK_GOOCH) {
    LimiterScalar_(ids, bc_model, bc_value, limiter_, [](double x) {
      return (x < 1.5) ? x - 4 * x * x * x / 27 : 1.0;
    });
    ApplyLimiter(limiter_);
  } else if (type_ == OPERATOR_LIMITER_TENSORIAL) {
    LimiterTensorial_(ids, bc_model, bc_value);
  } else if (type_ == OPERATOR_LIMITER_KUZMIN) {
    LimiterKuzmin_(ids, bc_model, bc_value);
  } else {
    Errors::Message msg("Unknown limiter");
    Exceptions::amanzi_throw(msg);
  }

  // apply safety factor
  if (cfl_ != 1.0) {
    limiter_->Scale(cfl_);
    data_->Scale(cfl_);
  }

  lifting->data()->ScatterMasterToGhosted("cell");
}


/* ******************************************************************
* The limiter must be between 0 and 1
****************************************************************** */
void
LimiterCell::ApplyLimiter(Teuchos::RCP<Epetra_MultiVector> limiter)
{
  int npoly = data_->NumVectors();

  for (int c = 0; c < ncells_owned_; c++) {
    for (int i = 0; i < npoly; i++) (*data_)[i][c] *= (*limiter)[0][c];
  }
}


/* *******************************************************************
* Tensorial limiter limits the gradient directly, to avoid
* calculation of a 3x3 matrix.
******************************************************************* */
void
LimiterCell::LimiterTensorial_(const AmanziMesh::cEntity_ID_View& ids,
                               const std::vector<int>& bc_model,
                               const std::vector<double>& bc_value)
{
  AMANZI_ASSERT(data_->NumVectors() == dim);
  auto& grad = *lifting_->data()->ViewComponent("cell");

  double u1, u1f, umin, umax, tol;
  double tol_base = sqrt(OPERATOR_LIMITER_TOLERANCE);

  AmanziGeometry::Point gradient_c1(dim);

  // Step 1: limit gradient to a feasiable set excluding Dirichlet boundary
  if (!external_bounds_) {
    if (stencil_id_ == OPERATOR_LIMITER_STENCIL_F2C)
      BoundsForFaces(*field_, bc_model, bc_value, stencil_id_, bounds_);
    else BoundsForCells(*field_, bc_model, bc_value, stencil_id_, bounds_);
  }

  for (int n = 0; n < ids.size(); ++n) {
    int c = ids[n];
    const auto& faces = mesh_->getCellFaces(c);
    int nfaces = faces.size();

    const AmanziGeometry::Point& xc = mesh_->getCellCentroid(c);
    for (int i = 0; i < dim; i++) gradient_c1[i] = grad[i][c];
    (*limiter_)[c] = norm(gradient_c1);

    u1 = (*field_)[component_][c];
    tol = tol_base * std::max(OPERATOR_LIMITER_FIELD_TOLERANCE, fabs(u1));

    // loop 1: project of planes
    // loop 2: verify and compute BJ factor
    double bj_fac(1.0);
    for (int loop = 0; loop < 2; loop++) {
      for (int i = 0; i < nfaces; ++i) {
        int f = faces[i];
        const AmanziGeometry::Point& xf = mesh_->getFaceCentroid(f);

        getBounds(c, f, stencil_id_, &umin, &umax);

        u1f = getValue(gradient_c1, c, xf);

        if (loop == 0) {
          if (u1f < umin) {
            ProjectOnPlane_(umin - u1, xf, xc, gradient_c1);
          } else if (u1f > umax) {
            ProjectOnPlane_(umax - u1, xf, xc, gradient_c1);
          }
        } else {
          if (u1f < umin - tol) {
            bj_fac = std::min(bj_fac, (umin - u1) / (u1f - u1));
          } else if (u1f > umax + tol) {
            bj_fac = std::min(bj_fac, (umax - u1) / (u1f - u1));
          }
        }
      }
    }

    if (bj_fac < 1.0) gradient_c1 *= bj_fac;
    double grad_norm = norm(gradient_c1);
    if (grad_norm < OPERATOR_LIMITER_TOLERANCE) gradient_c1.set(0.0);

    for (int i = 0; i < dim; i++) grad[i][c] = gradient_c1[i];
  }

  // Step 3: enforce a priori timestep estimate (division of dT by 2).
  if (limiter_correction_) {
    BoundsForCells(*field_, bc_model, bc_value, OPERATOR_LIMITER_STENCIL_C2C_CLOSEST, bounds_);
    LimiterExtensionTransportTensorial_();
  }

  // approximate estimate of scalar limiter (mainly for statistics)
  for (int c = 0; c < ncells_owned_; c++) {
    double grad_norm0 = (*limiter_)[c];
    if (grad_norm0 == 0.0) {
      (*limiter_)[c] = 1.0;
    } else {
      for (int i = 0; i < dim; i++) gradient_c1[i] = (*data_)[i][c];
      (*limiter_)[c] = norm(gradient_c1) / grad_norm0;
    }
  }
}


/* *******************************************************************
* Extension of the tensorial limiter. Routine changes gradient to
* satisfy an a prioty estimate of the stable timestep. That estimate
* assumes that the weigthed flux is smaller than the first-order flux.
******************************************************************* */
void
LimiterCell::LimiterExtensionTransportTensorial_()
{
  double u1f, u1;

  auto& grad_c = *lifting_->data()->ViewComponent("cell");
  auto& bounds_c = *bounds_->ViewComponent("cell", true);

  for (int c = 0; c < ncells_owned_; c++) {
    const auto& faces = mesh_->getCellFaces(c);
    int nfaces = faces.size();

    double a, b;
    double flux, outflux = 0.0, outflux_weigted = 0.0;
    for (int i = 0; i < nfaces; i++) {
      int f = faces[i];
      int c1 = (upwind_cells_[f].size() > 0) ? upwind_cells_[f][0] : -1;

      if (c == c1) {
        const AmanziGeometry::Point& xcf = mesh_->getFaceCentroid(f);
        u1f = getValue(c, xcf);
        u1 = (*field_)[component_][c];

        a = u1f - u1;
        if (fabs(a) > OPERATOR_LIMITER_TOLERANCE * (fabs(u1f) + fabs(u1))) {
          if (a > 0) b = u1 - bounds_c[0][c];
          else b = u1 - bounds_c[1][c];

          flux = fabs((*flux_)[0][f]);
          outflux += flux;
          if (b) {
            outflux_weigted += flux * a / b;
          } else {
            for (int k = 0; k < dim; k++) grad_c[k][c] = 0.0;
            break;
          }
        }
      }
    }

    if (outflux_weigted > outflux) {
      double psi = outflux / outflux_weigted;
      for (int i = 0; i < dim; i++) grad_c[i][c] *= psi;
    }
  }
}


/* *******************************************************************
* Routine calculates a scalar limiter using the cell-based algorithm.
* First, it limits face-centered value of a reconstracted function
* by min-max of two or more cell-centered values.
* Second, it limits outflux values which gives factor 0.5 in the
* timestep estimate.
******************************************************************* */
void
LimiterCell::LimiterScalar_(const AmanziMesh::cEntity_ID_View& ids,
                            const std::vector<int>& bc_model,
                            const std::vector<double>& bc_value,
                            Teuchos::RCP<Epetra_Vector> limiter,
                            double (*func)(double))
{
  limiter->PutScalar(1.0);

  double u1, u1x, u1_add, umin, umax;
  AmanziGeometry::Point gradient_c(dim);
  AmanziMesh::cEntity_ID_View ents;

  // limiting gradient inside domain
  if (!external_bounds_) {
    if (stencil_id_ == OPERATOR_LIMITER_STENCIL_F2C)
      BoundsForFaces(*field_, bc_model, bc_value, stencil_id_, bounds_);
    else BoundsForCells(*field_, bc_model, bc_value, stencil_id_, bounds_);
  }

  double tol_base = sqrt(OPERATOR_LIMITER_TOLERANCE);
  for (int n = 0; n < ids.size(); ++n) {
    int c = ids[n];
    u1 = (*field_)[component_][c];
    double tol = tol_base * std::max(OPERATOR_LIMITER_FIELD_TOLERANCE, fabs(u1));

    // only two options for control points are supported
    int x, nents(0);
    if (external_controls_) {
      nents = (*controls_)[c].size();
    } else if (location_ == AmanziMesh::Entity_kind::FACE) {
      ents = mesh_->getCellFaces(c);
      nents = ents.size();
    } else if (location_ == AmanziMesh::Entity_kind::NODE) {
      ents = mesh_->getCellNodes(c);
      nents = ents.size();
    } else {
      AMANZI_ASSERT(false);
    }

    for (int i = 0; i < nents; i++) {
      if (external_controls_) {
        x = c;
        u1_add = lifting_->getValueSlope(c, (*controls_)[c][i]);
      } else if (location_ == AmanziMesh::Entity_kind::FACE) {
        x = ents[i];
        const AmanziGeometry::Point& xf = mesh_->getFaceCentroid(x);
        u1_add = lifting_->getValueSlope(c, xf);
      } else if (location_ == AmanziMesh::Entity_kind::NODE) {
        x = ents[i];
        const auto& xv = mesh_->getNodeCoordinate(x);
        u1_add = lifting_->getValueSlope(c, xv);
      }

      u1x = u1 + u1_add;

      getBounds(c, x, stencil_id_, &umin, &umax);

      if (u1x < u1 - tol) {
        (*limiter)[c] = std::min((*limiter)[c], func((umin - u1) / u1_add));
      } else if (u1x > u1 + tol) {
        (*limiter)[c] = std::min((*limiter)[c], func((umax - u1) / u1_add));
      }
    }
  }

  // enforce an a priori timestep estimate (dT / 2).
  if (limiter_correction_) {
    if (stencil_id_ == OPERATOR_LIMITER_STENCIL_F2C)
      BoundsForCells(*field_, bc_model, bc_value, OPERATOR_LIMITER_STENCIL_C2C_CLOSEST, bounds_);
    LimiterExtensionTransportScalar_(limiter);
  }
}


/* *******************************************************************
* Extension of the scalarlimiter. Routine changes gradient to
* satisfy an a prioty estimate of the stable timestep. That estimate
* assumes that the weigthed flux is smaller that the first-order flux.
******************************************************************* */
void
LimiterCell::LimiterExtensionTransportScalar_(Teuchos::RCP<Epetra_Vector> limiter)
{
  AMANZI_ASSERT(upwind_cells_.size() > 0);

  double u1, u1f, tol;
  double tol_base = sqrt(OPERATOR_LIMITER_TOLERANCE);
  AmanziGeometry::Point gradient_c1(dim);

  auto& grad_c = *lifting_->data()->ViewComponent("cell");
  auto& bounds_c = *bounds_->ViewComponent("cell", true);

  for (int c = 0; c < ncells_owned_; c++) {
    const AmanziGeometry::Point& xc = mesh_->getCellCentroid(c);
    const auto& faces = mesh_->getCellFaces(c);
    int nfaces = faces.size();

    double a, b;
    double flux, outflux = 0.0, outflux_weigted = 0.0;
    for (int i = 0; i < nfaces; i++) {
      int f = faces[i];
      int c1 = (upwind_cells_[f].size() > 0) ? upwind_cells_[f][0] : -1;

      if (c == c1) {
        const AmanziGeometry::Point& xcf = mesh_->getFaceCentroid(f);
        u1 = (*field_)[component_][c];
        for (int k = 0; k < dim; k++) gradient_c1[k] = grad_c[k][c1];
        u1f = u1 + gradient_c1 * (xcf - xc);

        a = u1f - u1;
        tol = OPERATOR_LIMITER_TOLERANCE * (fabs(u1f) + fabs(u1));
        if (fabs(a) > tol) {
          if (a > 0) b = u1 - bounds_c[0][c];
          else b = u1 - bounds_c[1][c];

          flux = fabs((*flux_)[0][f]);
          outflux += flux;
          if (fabs(b) > tol * tol_base) {
            outflux_weigted += flux * a / b;
          } else {
            (*limiter)[c] = 0.0;
            break;
          }
        }
      }
    }

    if (outflux_weigted > outflux) {
      double psi = outflux / outflux_weigted;
      (*limiter)[c] *= psi;
    }
  }
}


/* *******************************************************************
* Kuzmin's limiter use all neighbors of a computational cell.
******************************************************************* */
void
LimiterCell::LimiterKuzmin_(const AmanziMesh::cEntity_ID_View& ids,
                            const std::vector<int>& bc_model,
                            const std::vector<double>& bc_value)
{
  AMANZI_ASSERT(data_->NumVectors() == dim);
  auto& grad = *lifting_->data()->ViewComponent("cell");

  // calculate local extrema at nodes
  if (!external_bounds_) BoundsForNodes(*field_, bc_model, bc_value, stencil_id_, bounds_);
  auto& bounds_v = *bounds_->ViewComponent("node", true);

  // limit reconstructed gradients at cell nodes
  AmanziGeometry::Point xp(dim);

  AmanziGeometry::Point gradient_c(dim), p(dim), normal_new(dim), direction(dim);
  std::vector<AmanziGeometry::Point> normals;

  for (int n = 0; n < ids.size(); ++n) {
    int c = ids[n];
    auto nodes = mesh_->getCellNodes(c);
    int nnodes = nodes.size();
    std::vector<double> field_min_cell(nnodes), field_max_cell(nnodes);

    for (int i = 0; i < nnodes; i++) {
      int v = nodes[i];
      field_min_cell[i] = bounds_v[0][v];
      field_max_cell[i] = bounds_v[1][v];
    }

    for (int i = 0; i < dim; i++) gradient_c[i] = grad[i][c];
    (*limiter_)[c] = norm(gradient_c);

    LimiterKuzminCell_(c, gradient_c, field_min_cell, field_max_cell);

    for (int i = 0; i < dim; i++) grad[i][c] = gradient_c[i];
  }

  // enforce an a priori timestep estimate (dT / 2).
  if (limiter_correction_) {
    AmanziMesh::Entity_ID_View cells;
    std::vector<double> field_local_min(ncells_wghost_);
    std::vector<double> field_local_max(ncells_wghost_);

    for (int c = 0; c < ncells_owned_; c++) {
      auto nodes = mesh_->getCellNodes(c);
      field_local_min[c] = field_local_max[c] = (*field_)[component_][c];

      for (int i = 0; i < nodes.size(); i++) {
        int v = nodes[i];
        field_local_min[c] = std::min(field_local_min[c], bounds_v[0][v]);
        field_local_max[c] = std::max(field_local_max[c], bounds_v[1][v]);
      }
    }

    LimiterExtensionTransportKuzmin_(field_local_min, field_local_max);
  }

  // approximate estimate of scalar limiter (mainly for statistics)
  for (int c = 0; c < ncells_owned_; c++) {
    double grad_norm0 = (*limiter_)[c];
    if (grad_norm0 == 0.0) {
      (*limiter_)[c] = 1.0;
    } else {
      for (int i = 0; i < dim; i++) gradient_c[i] = grad[i][c];
      (*limiter_)[c] = norm(gradient_c) / grad_norm0;
    }
  }
}


/* *******************************************************************
* Kuzmin's limiter use all neighbors of the given cell and limit
* gradient in this cell only.
******************************************************************* */
void
LimiterCell::LimiterKuzminCell_(int c,
                                AmanziGeometry::Point& gradient_c,
                                const std::vector<double>& field_node_min_c,
                                const std::vector<double>& field_node_max_c)
{
  double u1p, u1, bj_fac(1.0), tol;
  double tol_base = sqrt(OPERATOR_LIMITER_TOLERANCE);
  AmanziGeometry::Point xp(dim);

  auto nodes = mesh_->getCellNodes(c);
  int nnodes = nodes.size();

  u1 = (*field_)[component_][c];
  tol = tol_base * std::max(OPERATOR_LIMITER_FIELD_TOLERANCE, fabs(u1));

  const AmanziGeometry::Point& xc = mesh_->getCellCentroid(c);

  for (int loop = 0; loop < 2; loop++) {
    for (int i = 0; i < nnodes; i++) {
      int v = nodes[i];
      double umin = field_node_min_c[i];
      double umax = field_node_max_c[i];

      xp = mesh_->getNodeCoordinate(v);
      u1p = getValue(gradient_c, c, xp);

      if (loop == 0) {
        if (u1p < umin) {
          ProjectOnPlane_(umin - u1, xp, xc, gradient_c);
        } else if (u1p > umax) {
          ProjectOnPlane_(umax - u1, xp, xc, gradient_c);
        }
      } else {
        if (u1p < umin - tol) {
          bj_fac = std::min(bj_fac, (umin - u1) / (u1p - u1));
        } else if (u1p > umax + tol) {
          bj_fac = std::min(bj_fac, (umax - u1) / (u1p - u1));
        }
      }
    }

    if (bj_fac < 1.0) gradient_c *= bj_fac;
    double grad_norm = norm(gradient_c);
    if (grad_norm < OPERATOR_LIMITER_TOLERANCE) gradient_c.set(0.0);
  }
}


/* *******************************************************************
* Extension of Kuzmin's limiter. Routine changes gradient to
* satisfy an a prioty estimate of the stable timestep. That estimate
* assumes that the weigthed flux is smaller that the first-order flux.
******************************************************************* */
void
LimiterCell::LimiterExtensionTransportKuzmin_(const std::vector<double>& field_local_min,
                                              const std::vector<double>& field_local_max)
{
  AMANZI_ASSERT(upwind_cells_.size() > 0);

  double u1, up;
  AmanziGeometry::Point xp(dim);

  auto& grad = *lifting_->data()->ViewComponent("cell");

  for (int c = 0; c < ncells_owned_; c++) {
    const auto& faces = mesh_->getCellFaces(c);
    int nfaces = faces.size();

    double a, b;
    double flux, outflux = 0.0, outflux_weigted = 0.0;
    for (int i = 0; i < nfaces; i++) {
      int f = faces[i];
      int c1 = (upwind_cells_[f].size() > 0) ? upwind_cells_[f][0] : -1;

      if (c == c1) {
        auto nodes = mesh_->getFaceNodes(f);
        int nnodes = nodes.size();

        // define dimensionless quadrature weigths
        std::vector<double> weights(nnodes, 0.5);
        if (dim == 3) {
          double area = mesh_->getFaceArea(f);
          WhetStone::PolygonCentroidWeights(*mesh_, nodes, area, weights);
        }

        for (int j = 0; j < nnodes; j++) {
          int v = nodes[j];
          xp = mesh_->getNodeCoordinate(v);
          up = getValue(c, xp);
          u1 = (*field_)[component_][c];

          a = up - u1;
          if (fabs(a) > OPERATOR_LIMITER_TOLERANCE * (fabs(up) + fabs(u1))) {
            if (a > 0) b = u1 - field_local_min[c];
            else b = u1 - field_local_max[c];

            flux = fabs((*flux_)[0][f]) * weights[j];
            outflux += flux;
            if (b) {
              outflux_weigted += flux * a / b;
            } else {
              for (int k = 0; k < dim; k++) grad[k][c] = 0.0;
              break;
            }
          }
        }
      }
    }

    if (outflux_weigted > outflux) {
      double psi = outflux / outflux_weigted;
      for (int i = 0; i < dim; i++) grad[i][c] *= psi;
    }
  }
}


/* *******************************************************************
* Project gradient on a plane defined by du + grad * direction = 0.
******************************************************************* */
void
LimiterCell::ProjectOnPlane_(double du,
                             const AmanziGeometry::Point& xf,
                             const AmanziGeometry::Point& xc,
                             AmanziGeometry::Point& gradient)
{
  AmanziGeometry::Point direction = xf - xc;
  double len = norm(direction);
  direction /= len;

  double a = du / len - gradient * direction;
  gradient += a * direction;
}


/* *******************************************************************
* Identify flux direction based on orientation of the face normal
* and sign of the Darcy velocity.
******************************************************************* */
void
LimiterCell::IdentifyUpwindCells_()
{
  upwind_cells_.clear();
  downwind_cells_.clear();

  upwind_cells_.resize(nfaces_wghost_);
  downwind_cells_.resize(nfaces_wghost_);

  for (int c = 0; c < ncells_wghost_; c++) {
    const auto& [faces, dirs] = mesh_->getCellFacesAndDirections(c);

    for (int i = 0; i < faces.size(); i++) {
      int f = faces[i];
      double tmp = (*flux_)[0][f] * dirs[i];
      if (tmp > 0.0) {
        upwind_cells_[f].push_back(c);
      } else if (tmp < 0.0) {
        downwind_cells_[f].push_back(c);
      } else if (dirs[i] > 0) {
        upwind_cells_[f].push_back(c);
      } else {
        downwind_cells_[f].push_back(c);
      }
    }
  }
}


/* ******************************************************************
* Calculate internal bounds.
****************************************************************** */
void
LimiterCell::BoundsForCells(const Epetra_MultiVector& field,
                            const std::vector<int>& bc_model,
                            const std::vector<double>& bc_value,
                            int stencil,
                            Teuchos::RCP<CompositeVector>& bounds) const
{
  if (bounds == Teuchos::null || !bounds->HasComponent("cell")) {
    CompositeVectorSpace cvs;
    if (bounds != Teuchos::null) cvs = bounds->Map();
    cvs.SetMesh(mesh_)->SetGhosted(true)->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 2);
    bounds = Teuchos::rcp(new CompositeVector(cvs));
  }
  auto& bounds_c = *bounds->ViewComponent("cell", true);

  for (int c = 0; c < ncells_wghost_; ++c) {
    bounds_c[0][c] = OPERATOR_LIMITER_INFINITY;
    bounds_c[1][c] = -OPERATOR_LIMITER_INFINITY;
  }

  if (stencil == OPERATOR_LIMITER_STENCIL_C2C_CLOSEST) {
    for (int c = 0; c < ncells_owned_; c++) {
      double value = field[component_][c];
      bounds_c[0][c] = std::min(bounds_c[0][c], value);
      bounds_c[1][c] = std::max(bounds_c[1][c], value);

      auto cells = AmanziMesh::getCellFaceAdjacentCells(*mesh_, c, AmanziMesh::Parallel_kind::ALL);
      for (int i = 0; i < cells.size(); i++) {
        value = field[component_][cells[i]];
        bounds_c[0][c] = std::min(bounds_c[0][c], value);
        bounds_c[1][c] = std::max(bounds_c[1][c], value);
      }
    }
  } else if (stencil == OPERATOR_LIMITER_STENCIL_C2C_ALL) {
    for (int c = 0; c < ncells_owned_; c++) {
      double value = field[component_][c];
      bounds_c[0][c] = std::min(bounds_c[0][c], value);
      bounds_c[1][c] = std::max(bounds_c[1][c], value);

      auto nodes = mesh_->getCellNodes(c);
      for (int i = 0; i < nodes.size(); i++) {
        int v = nodes[i];
        auto cells = mesh_->getNodeCells(v);

        for (int k = 0; k < cells.size(); ++k) {
          value = field[component_][cells[k]];
          bounds_c[0][c] = std::min(bounds_c[0][c], value);
          bounds_c[1][c] = std::max(bounds_c[1][c], value);
        }
      }
    }
  }

  // add boundary conditions to the bounds
  if (bc_model.size() > 0) {
    for (int f = 0; f < nfaces_wghost_; ++f) {
      if (bc_model[f] == OPERATOR_BC_DIRICHLET) {
        auto cells = mesh_->getFaceCells(f);

        for (int n = 0; n < cells.size(); ++n) {
          int c = cells[n];
          bounds_c[0][c] = std::min(bounds_c[0][c], bc_value[f]);
          bounds_c[1][c] = std::max(bounds_c[1][c], bc_value[f]);
        }
      }
    }
  }
}


/* ******************************************************************
* Calculate internal bounds for the face to closest cells stencil.
****************************************************************** */
void
LimiterCell::BoundsForFaces(const Epetra_MultiVector& field,
                            const std::vector<int>& bc_model,
                            const std::vector<double>& bc_value,
                            int stencil,
                            Teuchos::RCP<CompositeVector>& bounds) const
{
  if (bounds == Teuchos::null || !bounds->HasComponent("face")) {
    CompositeVectorSpace cvs;
    if (bounds != Teuchos::null) cvs = bounds->Map();
    cvs.SetMesh(mesh_)->SetGhosted(true)->AddComponent("face", AmanziMesh::Entity_kind::FACE, 2);
    bounds = Teuchos::rcp(new CompositeVector(cvs));
  }
  auto& bounds_f = *bounds->ViewComponent("face", true);

  for (int f = 0; f < nfaces_wghost_; ++f) {
    bounds_f[0][f] = OPERATOR_LIMITER_INFINITY;
    bounds_f[1][f] = -OPERATOR_LIMITER_INFINITY;
  }

  for (int f = 0; f < nfaces_wghost_; ++f) {
    auto cells = mesh_->getFaceCells(f);

    for (int i = 0; i < cells.size(); ++i) {
      int c = cells[i];
      double value = field[component_][c];
      bounds_f[0][f] = std::min(bounds_f[0][f], value);
      bounds_f[1][f] = std::max(bounds_f[1][f], value);
    }
  }

  // add boundary conditions to the bounds
  if (bc_model.size() > 0) {
    for (int f = 0; f < nfaces_owned_; ++f) {
      if (bc_model[f] == OPERATOR_BC_DIRICHLET) {
        bounds_f[0][f] = std::min(bounds_f[0][f], bc_value[f]);
        bounds_f[1][f] = std::max(bounds_f[1][f], bc_value[f]);
      }
    }
  }
}


/* ******************************************************************
* Calculate internal bounds for the edge to closest cells stencil.
****************************************************************** */
void
LimiterCell::BoundsForEdges(const Epetra_MultiVector& field,
                            const std::vector<int>& bc_model,
                            const std::vector<double>& bc_value,
                            int stencil,
                            Teuchos::RCP<CompositeVector>& bounds) const
{
  if (bounds == Teuchos::null || !bounds->HasComponent("edge")) {
    CompositeVectorSpace cvs;
    if (bounds != Teuchos::null) cvs = bounds->Map();
    cvs.SetMesh(mesh_)->SetGhosted(true)->AddComponent("edge", AmanziMesh::Entity_kind::FACE, 2);
    bounds = Teuchos::rcp(new CompositeVector(cvs));
  }

  auto& bounds_e = *bounds->ViewComponent("edge", true);

  for (int e = 0; e < nedges_wghost_; ++e) {
    bounds_e[0][e] = OPERATOR_LIMITER_INFINITY;
    bounds_e[1][e] = -OPERATOR_LIMITER_INFINITY;
  }

  for (int e = 0; e < nedges_wghost_; ++e) {
    auto cells = mesh_->getEdgeCells(e);

    for (int i = 0; i < cells.size(); ++i) {
      int c = cells[i];
      double value = field[component_][c];
      bounds_e[0][e] = std::min(bounds_e[0][e], value);
      bounds_e[1][e] = std::max(bounds_e[1][e], value);
    }
  }

  // add boundary conditions to the bounds
  if (bc_model.size() > 0) {
    for (int e = 0; e < nedges_owned_; ++e) {
      if (bc_model[e] == OPERATOR_BC_DIRICHLET) {
        bounds_e[0][e] = std::min(bounds_e[0][e], bc_value[e]);
        bounds_e[1][e] = std::max(bounds_e[1][e], bc_value[e]);
      }
    }
  }
}


/* ******************************************************************
* Calculate internal bounds for the node to cells stencil.
****************************************************************** */
void
LimiterCell::BoundsForNodes(const Epetra_MultiVector& field,
                            const std::vector<int>& bc_model,
                            const std::vector<double>& bc_value,
                            int stencil,
                            Teuchos::RCP<CompositeVector>& bounds) const
{
  if (bounds == Teuchos::null || !bounds->HasComponent("node")) {
    CompositeVectorSpace cvs;
    if (bounds != Teuchos::null) cvs = bounds->Map();
    cvs.SetMesh(mesh_)->SetGhosted(true)->AddComponent("node", AmanziMesh::Entity_kind::FACE, 2);
    bounds = Teuchos::rcp(new CompositeVector(cvs));
  }
  auto& bounds_v = *bounds->ViewComponent("node", true);

  for (int v = 0; v < nnodes_wghost_; ++v) {
    bounds_v[0][v] = OPERATOR_LIMITER_INFINITY;
    bounds_v[1][v] = -OPERATOR_LIMITER_INFINITY;
  }

  for (int v = 0; v < nnodes_wghost_; ++v) {
    auto cells = mesh_->getNodeCells(v);

    for (int i = 0; i < cells.size(); ++i) {
      int c = cells[i];
      double value = field[component_][c];
      bounds_v[0][v] = std::min(bounds_v[0][v], value);
      bounds_v[1][v] = std::max(bounds_v[1][v], value);
    }
  }

  // add boundary conditions to the bounds
  if (bc_model.size() > 0) {
    for (int v = 0; v < nnodes_wghost_; ++v) {
      if (bc_model[v] == OPERATOR_BC_DIRICHLET) {
        bounds_v[0][v] = std::min(bounds_v[0][v], bc_value[v]);
        bounds_v[1][v] = std::max(bounds_v[1][v], bc_value[v]);
      }
    }
  }
}


/* ******************************************************************
* Extract face limiters.
****************************************************************** */
void
LimiterCell::getBounds(int c, int f, int stencil, double* umin, double* umax)
{
  if (stencil == OPERATOR_LIMITER_STENCIL_F2C) {
    auto& bounds_f = *bounds_->ViewComponent("face", true);
    *umin = bounds_f[0][f];
    *umax = bounds_f[1][f];
  } else {
    auto& bounds_c = *bounds_->ViewComponent("cell", true);
    *umin = bounds_c[0][c];
    *umax = bounds_c[1][c];
  }
}


/* ******************************************************************
* Calculates reconstructed value at point p.
****************************************************************** */
double
LimiterCell::getValue(int c, const AmanziGeometry::Point& p)
{
  const auto& xc = mesh_->getCellCentroid(c);

  double value = (*field_)[component_][c];
  for (int i = 0; i < dim; i++) value += (*data_)[i][c] * (p[i] - xc[i]);
  return value;
}


/* ******************************************************************
* Calculates reconstructed value at point p using external gradient.
****************************************************************** */
double
LimiterCell::getValue(const AmanziGeometry::Point& gradient, int c, const AmanziGeometry::Point& p)
{
  const auto& xc = mesh_->getCellCentroid(c);
  return (*field_)[component_][c] + gradient * (p - xc);
}

} // namespace Operators
} // namespace Amanzi
