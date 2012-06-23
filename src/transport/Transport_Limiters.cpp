/*
This is the transport component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided Reconstruction.cppin the top-level COPYRIGHT file.

Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <algorithm>
#include <vector>

#include "Transport_PK.hpp"

namespace Amanzi {
namespace AmanziTransport {

/* *******************************************************************
 * Tensorial limiter limits the gradient directly, to avoid 
 * calculation of a 3x3 matrix.
 ****************************************************************** */
void Transport_PK::LimiterTensorial(const int component,
                                    Teuchos::RCP<Epetra_Vector> scalar_field,
                                    Teuchos::RCP<Epetra_MultiVector> gradient)
{
  const Epetra_Vector& darcy_flux = TS_nextBIG->ref_darcy_flux();

  double u1, u2, u1f, u2f, umin, umax, L22normal_new;
  AmanziGeometry::Point gradient_c1(dim), gradient_c2(dim);
  AmanziGeometry::Point normal_new(dim), direction(dim), p(dim);

  std::vector<AmanziGeometry::Point> normals;
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> fdirs;

  // Step 1: limit gradient to a feasiable set excluding Dirichlet boundary
  for (int c = 0; c <= cmax_owned; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &fdirs);
    int nfaces = faces.size();

    const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
    for (int i = 0; i < dim; i++) gradient_c1[i] = (*gradient)[i][c];

    normals.clear();  // normals to planes the define the feasiable set
    for (int loop = 0; loop < 2; loop++) {
      for (int i = 0; i < nfaces; i++) {
        int f = faces[i];
        int c1 = (*upwind_cell_)[f];
        int c2 = (*downwind_cell_)[f];

        if (c1 >= 0 && c2 >= 0) {
          u1 = (*scalar_field)[c1];
          u2 = (*scalar_field)[c2];
        } else if (c1 == c) {
          u1 = u2 = (*scalar_field)[c1];
        } else {
          continue;
        }
        umin = std::min(u1, u2);
        umax = std::max(u1, u2);

        const AmanziGeometry::Point& xcf = mesh_->face_centroid(f);
        u1f = lifting.getValue(gradient_c1, c, xcf);

        // check if umin <= u1f <= umax
        if (u1f < umin) {
          normal_new = xcf - xc;
          CalculateDescentDirection(normals, normal_new, L22normal_new, direction);

          p = ((umin - u1) / sqrt(L22normal_new)) * direction;
          ApplyDirectionalLimiter(normal_new, p, direction, gradient_c1);
        } else if (u1f > umax) {
          normal_new = xcf - xc;
          CalculateDescentDirection(normals, normal_new, L22normal_new, direction);

          p = ((umax - u1) / sqrt(L22normal_new)) * direction;
          ApplyDirectionalLimiter(normal_new, p, direction, gradient_c1);
        }
      }
      if (normals.size() == 0) break;  // No limiters were imposed.
    }

    double grad_norm = norm(gradient_c1);
    if (grad_norm < TRANSPORT_LIMITER_TOLERANCE * bc_scaling) gradient_c1.set(0.0);

    for (int i = 0; i < dim; i++) (*gradient)[i][c] = gradient_c1[i];
  }

  // Local extrema are calculated here and updated in Step 2.
  AmanziMesh::Entity_ID_List cells;

  for (int c = 0; c <= cmax_owned; c++) {
    mesh_->cell_get_face_adj_cells(c, AmanziMesh::USED, &cells);
    component_local_min_[c] = component_local_max_[c] = (*scalar_field)[c];

    for (int i = 0; i < cells.size(); i++) {
      double value = (*scalar_field)[cells[i]];
      component_local_min_[c] = std::min(component_local_min_[c], value);
      component_local_max_[c] = std::max(component_local_max_[c], value);
    }
  }

  // Step 2: limit gradient on the Dirichlet boundary
  for (int n = 0; n < bcs.size(); n++) {
    if (component == bcs_tcc_index[n]) {
      for (Amanzi::Iterator bc = bcs[n]->begin(); bc != bcs[n]->end(); ++bc) {
        int f = bc->first;
        int c1 = (*upwind_cell_)[f];
        int c2 = (*downwind_cell_)[f];

        if (c2 >= 0 && c2 <= cmax_owned) {
          u2 = (*scalar_field)[c2];
          u1 = bc->second;
          umin = std::min(u1, u2);
          umax = std::max(u1, u2);

          bc_scaling = std::max(bc_scaling, u1);
          component_local_max_[c2] = std::max(component_local_max_[c2], u1);
          component_local_min_[c2] = std::min(component_local_min_[c2], u1);

          const AmanziGeometry::Point& xc2 = mesh_->cell_centroid(c2);
          const AmanziGeometry::Point& xcf = mesh_->face_centroid(f);
          u2f = lifting.getValue(c2, xcf);
          for (int i = 0; i < dim; i++) gradient_c2[i] = (*gradient)[i][c2];
          direction = xcf - xc2;

          if (u2f < umin) {
            p = ((umin - u2) / L22(direction)) * direction;
            ApplyDirectionalLimiter(direction, p, direction, gradient_c2);

          } else if (u2f > umax) {
            p = ((umax - u2) / L22(direction)) * direction;
            ApplyDirectionalLimiter(direction, p, direction, gradient_c2);
          }

          for (int i = 0; i < dim; i++) (*gradient)[i][c2] = gradient_c2[i];
        }
      }
    }
  }

  // Step 3: enforcing a priori time step estimate (division of dT by 2).
  for (int c = 0; c <= cmax_owned; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &fdirs);
    int nfaces = faces.size();

    double a, b;
    double flux, outflux = 0.0, outflux_weigted = 0.0;
    for (int i = 0; i < nfaces; i++) {
      int f = faces[i];
      int c1 = (*upwind_cell_)[f];

      if (c == c1) {
        const AmanziGeometry::Point& xcf = mesh_->face_centroid(f);
        u1f = lifting.getValue(c, xcf);
        u1 = (*scalar_field)[c];

        a = u1f - u1;
        if (fabs(a) > TRANSPORT_LIMITER_TOLERANCE * (fabs(u1f) + fabs(u1))) {
          if (a > 0) b = u1 - component_local_min_[c];
          else       b = u1 - component_local_max_[c];

          flux = fabs(darcy_flux[f]);
          outflux += flux;
          if (b) {
            outflux_weigted += flux * a / b;
          } else {
            for (int k = 0; k < dim; k++) (*gradient)[k][c] = 0.0;
            break;
          }
        }
      }
    }

    if (outflux_weigted > outflux) {
      double psi = outflux / outflux_weigted;
      for (int i = 0; i < dim; i++) (*gradient)[i][c] *= psi;
    }
  }

  TS_nextBIG->CopyMasterMultiCell2GhostMultiCell(*gradient);
}


/* *******************************************************************
 * Routine calculates BJ limiter adjusted to linear advection. 
 * First, it limits face-centered value of a reconstracted function 
 * by min-max of two cell-centered values.
 * Second, it limits outflux values which gives factor 0.5 in the
 * time step estimate.
 ****************************************************************** */
void Transport_PK::LimiterBarthJespersen(const int component,
                                         Teuchos::RCP<Epetra_Vector> scalar_field,
                                         Teuchos::RCP<Epetra_MultiVector> gradient,
                                         Teuchos::RCP<Epetra_Vector> limiter)
{
  const Epetra_Vector& darcy_flux = TS_nextBIG->ref_darcy_flux();

  for (int c = 0; c <= cmax; c++) (*limiter)[c] = 1.0;

  double u1, u2, u1f, u2f, umin, umax;  // cell and inteface values
  AmanziGeometry::Point gradient_c1(dim), gradient_c2(dim);

  // Step 1: limiting gradient inside domain
  for (int f = 0; f <= fmax_owned; f++) {
    int c1, c2;
    c1 = (*upwind_cell_)[f];
    c2 = (*downwind_cell_)[f];
    if (c1 < 0 || c2 < 0) continue;

    u1 = (*scalar_field)[c1];
    u2 = (*scalar_field)[c2];
    umin = std::min(u1, u2);
    umax = std::max(u1, u2);
    double tol = sqrt(TRANSPORT_LIMITER_TOLERANCE) * (fabs(u1) + fabs(u2));

    const AmanziGeometry::Point& xc1 = mesh_->cell_centroid(c1);
    const AmanziGeometry::Point& xc2 = mesh_->cell_centroid(c2);
    const AmanziGeometry::Point& xcf = mesh_->face_centroid(f);

    for (int i = 0; i < dim; i++) gradient_c1[i] = (*gradient)[i][c1];
    double u1_add = gradient_c1 * (xcf - xc1);
    u1f = u1 + u1_add;

    if (u1f < umin - tol) {
      (*limiter)[c1] = std::min((*limiter)[c1], (umin - u1) / u1_add);
    } else if (u1f > umax + tol) {
      (*limiter)[c1] = std::min((*limiter)[c1], (umax - u1) / u1_add);
    }

    for (int i = 0; i < dim; i++) gradient_c2[i] = (*gradient)[i][c2];
    double u2_add = gradient_c2 * (xcf - xc2);
    u2f = u2 + u2_add;

    if (u2f < umin - tol) {
      (*limiter)[c2] = std::min((*limiter)[c2], (umin - u2) / u2_add);
    } else if (u2f > umax + tol) {
      (*limiter)[c2] = std::min((*limiter)[c2], (umax - u2) / u2_add);
    }
  }

  // Local extrema are calculated here and updated in Step 2.
  AmanziMesh::Entity_ID_List cells;

  for (int c = 0; c <= cmax_owned; c++) {
    mesh_->cell_get_face_adj_cells(c, AmanziMesh::USED, &cells);
    component_local_min_[c] = component_local_max_[c] = (*scalar_field)[c];

    for (int i = 0; i < cells.size(); i++) {
      double value = (*scalar_field)[cells[i]];
      component_local_min_[c] = std::min(component_local_min_[c], value);
      component_local_max_[c] = std::max(component_local_max_[c], value);
    }
  }

  // Step 2: limiting gradient on the Dirichlet boundary
  for (int n = 0; n < bcs.size(); n++) {
    if (component == bcs_tcc_index[n]) {
      for (Amanzi::Iterator bc = bcs[n]->begin(); bc != bcs[n]->end(); ++bc) {
        int f = bc->first;
        int c2 = (*downwind_cell_)[f];

        if (c2 >= 0) {
          u2 = (*scalar_field)[c2];
          u1 = bc->second;
          umin = std::min(u1, u2);
          umax = std::max(u1, u2);
          double tol = sqrt(TRANSPORT_LIMITER_TOLERANCE) * (fabs(u1) + fabs(u2));

          component_local_max_[c2] = std::max(component_local_max_[c2], u1);
          component_local_min_[c2] = std::min(component_local_min_[c2], u1);

          const AmanziGeometry::Point& xc2 = mesh_->cell_centroid(c2);
          const AmanziGeometry::Point& xcf = mesh_->face_centroid(f);

          for (int i = 0; i < dim; i++) gradient_c2[i] = (*gradient)[i][c2];
          double u_add = gradient_c2 * (xcf - xc2);
          u2f = u2 + u_add;

          if (u2f < umin - tol) {
            (*limiter)[c2] = std::min((*limiter)[c2], (umin - u2) / u_add);
          } else if (u2f > umax + tol) {
            (*limiter)[c2] = std::min((*limiter)[c2], (umax - u2) / u_add);
          }
        }
      }
    }
  }

  // Step 3: enforcing a priori time step estimate (division of dT by 2).
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> fdirs;

  for (int c = 0; c <= cmax_owned; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &fdirs);
    int nfaces = faces.size();

    double a, b;
    double flux, outflux = 0.0, outflux_weigted = 0.0;
    for (int i = 0; i < nfaces; i++) {
      int f = faces[i];
      int c1 = (*upwind_cell_)[f];

      if (c == c1) {
        const AmanziGeometry::Point& xcf = mesh_->face_centroid(f);
        u1f = lifting.getValue(c, xcf);
        u1 = (*scalar_field)[c];

        a = u1f - u1;
        if (fabs(a) > TRANSPORT_LIMITER_TOLERANCE * (fabs(u1f) + fabs(u1))) {
          if (a > 0) b = u1 - component_local_min_[c];
          else       b = u1 - component_local_max_[c];

          flux = fabs(darcy_flux[f]);
          outflux += flux;
          if (b) {
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

  TS_nextBIG->CopyMasterCell2GhostCell(*limiter);
}


/* *******************************************************************
 * Kuzmin's limiter use all neighbors of a computational cell.  
 ****************************************************************** */
void Transport_PK::LimiterKuzmin(const int component,
                                 Teuchos::RCP<Epetra_Vector> scalar_field,
                                 Teuchos::RCP<Epetra_MultiVector> gradient)
{
  // Step 1: local extrema are calculated here at nodes and updated later
  std::vector<double> component_node_min(vmax + 1);
  std::vector<double> component_node_max(vmax + 1);

  AmanziMesh::Entity_ID_List nodes;

  component_node_min.assign(vmax + 1,  TRANSPORT_CONCENTRATION_INFINITY);
  component_node_max.assign(vmax + 1, -TRANSPORT_CONCENTRATION_INFINITY);

  for (int c = 0; c <= cmax_owned; c++) {
    mesh_->cell_get_nodes(c, &nodes);

    double value = (*scalar_field)[c];
    for (int i = 0; i < nodes.size(); i++) {
      int v = nodes[i];
      component_node_min[v] = std::min(component_node_min[v], value);
      component_node_max[v] = std::max(component_node_max[v], value);
    }
  }

  // Update min/max at nodes from influx boundary data
  for (int n = 0; n < bcs.size(); n++) {
    if (component == bcs_tcc_index[n]) {
      for (Amanzi::Iterator bc = bcs[n]->begin(); bc != bcs[n]->end(); ++bc) {
        int f = bc->first;
        int c2 = (*downwind_cell_)[f];

        if (c2 >= 0) {
          mesh_->face_get_nodes(f, &nodes);
          double value = (*scalar_field)[c2];

          for (int i = 0; i < nodes.size(); i++) {
            int v = nodes[i];
            component_node_min[v] = std::min(component_node_min[v], value);
            component_node_max[v] = std::max(component_node_max[v], value);
          }
        }
      }
    }
  }

  // Step 2: limit reconstructed gradients at cell nodes
  double umin, umax, up, u1;
  AmanziGeometry::Point xp(dim);

  double L22normal_new;
  AmanziGeometry::Point gradient_c(dim), p(dim), normal_new(dim), direction(dim);

  std::vector<AmanziGeometry::Point> normals;

  for (int c = 0; c <= cmax_owned; c++) {
    mesh_->cell_get_nodes(c, &nodes);
    int nnodes = nodes.size();

    const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
    for (int i = 0; i < dim; i++) gradient_c[i] = (*gradient)[i][c];

    double umin = component_node_min[c];
    double umax = component_node_max[c];

    normals.clear();  // normals to planes the define the feasiable set
    for (int loop = 0; loop < 2; loop++) {
      for (int i = 0; i < nnodes; i++) {
        int v = nodes[i];

        mesh_->node_get_coordinates(v, &xp);
        up = lifting.getValue(gradient_c, c, xp);

        // check if umin <= up <= umax
        if (up < umin) {
          normal_new = xp - xc;
          CalculateDescentDirection(normals, normal_new, L22normal_new, direction);

          p = ((umin - up) / sqrt(L22normal_new)) * direction;
          ApplyDirectionalLimiter(normal_new, p, direction, gradient_c);
        } else if (up > umax) {
          normal_new = xp - xc;
          CalculateDescentDirection(normals, normal_new, L22normal_new, direction);

          p = ((umax - up) / sqrt(L22normal_new)) * direction;
          ApplyDirectionalLimiter(normal_new, p, direction, gradient_c);
        }
      }
      if (normals.size() == 0) break;  // No limiters were imposed.
    }

    double grad_norm = norm(gradient_c);
    if (grad_norm < TRANSPORT_LIMITER_TOLERANCE * bc_scaling) gradient_c.set(0.0);

    for (int i = 0; i < dim; i++) (*gradient)[i][c] = gradient_c[i];
  }

  // Step 3: extrema are calculated for cells.
  AmanziMesh::Entity_ID_List cells;

  for (int c = 0; c <= cmax_owned; c++) {
    mesh_->cell_get_nodes(c, &nodes);
    component_local_min_[c] = component_local_max_[c] = (*scalar_field)[c];

    for (int i = 0; i < nodes.size(); i++) {
      int v = nodes[i];
      component_local_min_[c] = std::min(component_local_min_[c], component_node_min[v]);
      component_local_max_[c] = std::max(component_local_max_[c], component_node_max[v]);
    }
  }

  // Step 4: enforcing a priori time step estimate (division of dT by 2).
  // Experimental version is limited to 2D (lipnikov@lanl.gov).
  const Epetra_Vector& darcy_flux = TS_nextBIG->ref_darcy_flux();
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> fdirs;

  for (int c = 0; c <= cmax_owned; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &fdirs);
    int nfaces = faces.size();

    double a, b;
    double flux, outflux = 0.0, outflux_weigted = 0.0;
    for (int i = 0; i < nfaces; i++) {
      int f = faces[i];
      int c1 = (*upwind_cell_)[f];

      if (c == c1) {
        mesh_->face_get_nodes(f, &nodes);

        for (int j = 0; j < nodes.size(); j++) {
          int v = nodes[j];
          mesh_->node_get_coordinates(v, &xp);
          up = lifting.getValue(c, xp);
          u1 = (*scalar_field)[c];

          a = up - u1;
          if (fabs(a) > TRANSPORT_LIMITER_TOLERANCE * (fabs(up) + fabs(u1))) {
            if (a > 0)
              b = u1 - component_local_min_[c];
            else
              b = u1 - component_local_max_[c];

            flux = fabs(darcy_flux[f]);
            outflux += flux;
            if (b) {
              outflux_weigted += flux * a / b;
            } else {
              for (int k = 0; k < dim; k++) (*gradient)[k][c] = 0.0;
              break;
            }
          }
        }
      }
    }

    if (outflux_weigted > outflux) {
      double psi = outflux / outflux_weigted;
      for (int i = 0; i < dim; i++) (*gradient)[i][c] *= psi;
    }
  }

  TS_nextBIG->CopyMasterMultiCell2GhostMultiCell(*gradient);
}


/* *******************************************************************
 * Descent direction is obtained by orthogonalizing normal direction
 * 'normal_new' to previous normals. A few exceptions are analyzed.  
 ****************************************************************** */
void Transport_PK::CalculateDescentDirection(std::vector<AmanziGeometry::Point>& normals,
                                             AmanziGeometry::Point& normal_new,
                                             double& L22normal_new,
                                             AmanziGeometry::Point& direction)
{
  L22normal_new = L22(normal_new);
  normal_new /= sqrt(L22normal_new);
  direction = normal_new;

  int nnormals = normals.size();
  if (nnormals == dim) {
    normals.clear();
  } else {
    double a;
    for (int n = 0; n < nnormals; n++) {
      a = normals[n] * direction;
      for (int i = 0; i < dim; i++) direction[i] -= a * normals[n][i];
    }

    // verify new direction
    a = L22(direction);
    if (a < L22normal_new * TRANSPORT_LIMITER_TOLERANCE) {
      normals.clear();
      direction = normal_new;
    } else {
      a = sqrt(a);
      for (int i = 0; i < dim; i++) direction[i] /= a;
    }
  }
  normals.push_back(direction);
}


/* *******************************************************************
 * Routine projects gradient on a plane defined by normal and point p.
 ****************************************************************** */
void Transport_PK::ApplyDirectionalLimiter(AmanziGeometry::Point& normal,
                                           AmanziGeometry::Point& p,
                                           AmanziGeometry::Point& direction,
                                           AmanziGeometry::Point& gradient)
{
  double a = ((p - gradient) * normal) / (direction * normal);
  gradient += a * direction;
}

}  // namespace AmanziTransport
}  // namespace Amanzi

