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
  2. Limiters are modified optionally so the the stable time step
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
#include "exceptions.hh"
#include "Quadrature1D.hh"
#include "WhetStoneDefs.hh"
#include "WhetStoneMeshUtils.hh"

// Amanzi::Operators
#include "OperatorDefs.hh"
#include "LimiterCellDG.hh"

namespace Amanzi {
namespace Operators {

/* *******************************************************************
* Work in progress: initialization like in Transport PK.
******************************************************************* */
LimiterCellDG::LimiterCellDG(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh)
  : LimiterCell(mesh){};


/* ******************************************************************
* Apply internal limiter for DG schemes
****************************************************************** */
void
LimiterCellDG::ApplyLimiterDG(const AmanziMesh::Entity_ID_List& ids,
                              Teuchos::RCP<const Epetra_MultiVector> field,
                              const WhetStone::DG_Modal& dg,
                              const std::vector<int>& bc_model,
                              const std::vector<double>& bc_value)
{
  field_ = field;

  if (external_bounds_ && bounds_ == Teuchos::null) {
    Errors::Message msg("External bounds for limiters are requested but not provided");
    Exceptions::amanzi_throw(msg);
  }

  limiter_ = Teuchos::rcp(new Epetra_Vector(mesh_->getMap(AmanziMesh::Entity_kind::CELL, false)));
  limiter_->PutScalar(1.0);

  if (type_ == OPERATOR_LIMITER_BARTH_JESPERSEN_DG) {
    LimiterScalarDG_(dg, ids, bc_model, bc_value, [](double x) { return std::min(1.0, x); });
  } else if (type_ == OPERATOR_LIMITER_MICHALAK_GOOCH_DG) {
    LimiterScalarDG_(dg, ids, bc_model, bc_value, [](double x) {
      return (x < 1.5) ? x - 4 * x * x * x / 27 : 1.0;
    });
  } else if (type_ == OPERATOR_LIMITER_BARTH_JESPERSEN_DG_HIERARCHICAL) {
    LimiterHierarchicalDG_(dg, ids, bc_model, bc_value, [](double x) { return x; });
  } else {
    Errors::Message msg("Unknown limiter");
    Exceptions::amanzi_throw(msg);
  }
}


/* *******************************************************************
* The scalar limiter for modal DG schemes.
******************************************************************* */
void
LimiterCellDG::LimiterScalarDG_(const WhetStone::DG_Modal& dg,
                                const AmanziMesh::Entity_ID_List& ids,
                                const std::vector<int>& bc_model,
                                const std::vector<double>& bc_value,
                                double (*func)(double))
{
  AMANZI_ASSERT(dg.cell_basis(0).id() == WhetStone::TAYLOR_BASIS_NORMALIZED_ORTHO);

  double u1, u1f, umin, umax;
  AmanziMesh::Entity_ID_List nodes;

  int nk = field_->NumVectors();
  WhetStone::DenseVector data(nk);
  AmanziGeometry::Point x1(dim), x2(dim), xm(dim);
  int order = WhetStone::PolynomialSpaceOrder(dim, nk);

  if (!external_bounds_) { bounds_ = BoundsForCells(*field_, bc_model, bc_value, stencil_id_); }

  int ilast = (dim == 2) ? 1 : 0;

  for (int n = 0; n < ids.size(); ++n) {
    int c = ids[n];
    const auto& faces = mesh_->getCellFaces(c);
    int nfaces = faces.size();

    for (int i = 0; i < nk; ++i) data(i) = (*field_)[i][c];
    auto poly = dg.cell_basis(c).CalculatePolynomial(mesh_, c, order, data);

    u1 = (*field_)[0][c];
    double tol = sqrt(OPERATOR_LIMITER_TOLERANCE) * fabs(u1);

    int mq = limiter_points_ - 1;
    double climiter(1.0);

    for (int m = 0; m < nfaces; ++m) {
      int f = faces[m];
      nodes = mesh_->getFaceNodes(f);
      int nnodes = nodes.size();

      getBounds(c, f, stencil_id_, &umin, &umax);

      for (int i = 0; i < nnodes - ilast; ++i) {
        int j = (i + 1) % nnodes;
        x1 = mesh_->getNodeCoordinate(nodes[i]);
        x2 = mesh_->getNodeCoordinate(nodes[j]);

        for (int k = 0; k < limiter_points_; ++k) {
          xm = x1 * WhetStone::q1d_points[mq][k] + x2 * (1.0 - WhetStone::q1d_points[mq][k]);
          u1f = poly.Value(xm);
          double u1_add = u1f - u1;

          if (u1f < u1 - tol) {
            climiter = std::min(climiter, func((umin - u1) / u1_add));
          } else if (u1f > u1 + tol) {
            climiter = std::min(climiter, func((umax - u1) / u1_add));
          }
        }
      }
    }

    (*limiter_)[c] = climiter;
  }
}


/* *******************************************************************
* The hierarchical limiter for modal DG schemes.
******************************************************************* */
void
LimiterCellDG::LimiterHierarchicalDG_(const WhetStone::DG_Modal& dg,
                                      const AmanziMesh::Entity_ID_List& ids,
                                      const std::vector<int>& bc_model,
                                      const std::vector<double>& bc_value,
                                      double (*func)(double))
{
  AMANZI_ASSERT(dim == 2);
  AMANZI_ASSERT(dg.cell_basis(0).id() == WhetStone::TAYLOR_BASIS_NORMALIZED_ORTHO);

  double u1, u1f, umin, umax;
  AmanziMesh::Entity_ID_List nodes;

  int nk = field_->NumVectors();
  WhetStone::DenseVector data(nk);
  AmanziGeometry::Point x1(dim), x2(dim), xm(dim);
  int order = WhetStone::PolynomialSpaceOrder(dim, nk);

  // calculate bounds
  // -- for mean values
  component_ = 0;
  std::vector<Teuchos::RCP<CompositeVector>> bounds(1 + dim);
  bounds[0] = BoundsForCells(*field_, bc_model, bc_value, stencil_id_);

  // -- for gradient
  {
    std::vector<int> bc_model_none(nfaces_wghost_, OPERATOR_BC_NONE);
    Epetra_MultiVector field_tmp(mesh_->getMap(AmanziMesh::Entity_kind::CELL, true), dim);

    for (int c = 0; c < ncells_owned_; ++c) {
      for (int i = 0; i < nk; ++i) data(i) = (*field_)[i][c];
      auto poly = dg.cell_basis(c).CalculatePolynomial(mesh_, c, order, data);
      for (int i = 0; i < dim; ++i) field_tmp[i][c] = poly(i + 1);
    }

    for (int i = 0; i < dim; ++i) {
      component_ = i;
      bounds[i + 1] = BoundsForCells(field_tmp, bc_model_none, bc_value, stencil_id_);
    }
  }

  for (int n = 0; n < ids.size(); ++n) {
    int c = ids[n];
    const auto& faces = mesh_->getCellFaces(c);
    int nfaces = faces.size();

    for (int i = 0; i < nk; ++i) data(i) = (*field_)[i][c];
    auto poly = dg.cell_basis(c).CalculatePolynomial(mesh_, c, order, data);
    auto grad = Gradient(poly);
    // poly.Reshape(dim, 1);

    u1 = (*field_)[0][c];
    double tol = sqrt(OPERATOR_LIMITER_TOLERANCE) * fabs(u1);

    int m = limiter_points_ - 1;
    double climiter(1.0), hlimiter(1.0);

    for (int i = 0; i < nfaces; i++) {
      int f = faces[i];
      nodes = mesh_->getFaceNodes(f);

      x1 = mesh_->getNodeCoordinate(nodes[0]);
      x2 = mesh_->getNodeCoordinate(nodes[1]);

      // limit mean values
      bounds_ = bounds[0];
      getBounds(c, f, stencil_id_, &umin, &umax);

      for (int k = 0; k < limiter_points_; ++k) {
        xm = x1 * WhetStone::q1d_points[m][k] + x2 * (1.0 - WhetStone::q1d_points[m][k]);
        u1f = poly.Value(xm);
        double u1_add = u1f - u1;

        if (u1f < umin - tol) {
          climiter = std::min(climiter, func((umin - u1) / u1_add));
        } else if (u1f > umax + tol) {
          climiter = std::min(climiter, func((umax - u1) / u1_add));
        }
      }

      // limit gradient values
      for (int l = 0; l < dim; ++l) {
        bounds_ = bounds[l + 1];
        getBounds(c, f, stencil_id_, &umin, &umax);

        for (int k = 0; k < limiter_points_; ++k) {
          xm = x1 * WhetStone::q1d_points[m][k] + x2 * (1.0 - WhetStone::q1d_points[m][k]);
          u1f = grad[l].Value(xm);
          double u1_add = u1f - grad[l](0);

          if (u1f < umin) {
            hlimiter = std::min(hlimiter, func((umin - u1) / u1_add));
          } else if (u1f > umax) {
            hlimiter = std::min(hlimiter, func((umax - u1) / u1_add));
          }
        }
      }
    }

    climiter = std::max(climiter, hlimiter);
    (*limiter_)[c] = climiter;
  }
}

} // namespace Operators
} // namespace Amanzi
