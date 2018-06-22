/*
  Operators 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Implementation of three different limiters.
*/

#include <algorithm>
#include <vector>

// TPLs
#include "Epetra_Vector.h"

// Amanzi
#include "MeshUtils.hh"

// Amanzi::Operators
#include "OperatorDefs.hh"
#include "ReconstructionCell.hh"

namespace Amanzi {
namespace Operators {

/* *******************************************************************
* Work in progress: initialization like in Transport PK.
******************************************************************* */
void ReconstructionCell::InitLimiter(Teuchos::RCP<const Epetra_MultiVector> flux)
{
  flux_ = flux;
  IdentifyUpwindCells_();
}


/* *******************************************************************
* Tensorial limiter limits the gradient directly, to avoid 
* calculation of a 3x3 matrix.
******************************************************************* */
void ReconstructionCell::LimiterTensorial_(
    const std::vector<int>& bc_model, const std::vector<double>& bc_value)
{
  double u1, u2, u1f, u2f, umin, umax, L22normal_new;
  AmanziGeometry::Point gradient_c1(dim), gradient_c2(dim);
  AmanziGeometry::Point normal_new(dim), direction(dim), p(dim);

  Epetra_MultiVector& grad = *gradient_->ViewComponent("cell", false);

  std::vector<AmanziGeometry::Point> normals;
  AmanziMesh::Entity_ID_List faces;
  Teuchos::RCP<Epetra_Vector> limiter = Teuchos::rcp(new Epetra_Vector(mesh_->cell_map(false)));

  // Step 1: limit gradient to a feasiable set excluding Dirichlet boundary
  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
    for (int i = 0; i < dim; i++) gradient_c1[i] = grad[i][c];
    (*limiter_)[c] = norm(gradient_c1); 

    normals.clear();  // normals to planes that define a feasiable set
    for (int loop = 0; loop < 2; loop++) {
      for (int i = 0; i < nfaces; ++i) {
        int f = faces[i];
        // Limiting on inflow boundary is done in Step 2.
        if (upwind_cells_[f].size() == 0) continue;

        // We enforce non-decreasing reconstruction on outflow boundary.
        if (downwind_cells_[f].size() == 0) {
          u1 = u2 = (*field_)[component_][c];
        } else {
          int c1 = upwind_cells_[f][0];
          for (int n = 1; n < upwind_cells_[f].size(); ++n) {
            if (upwind_cells_[f][n] == c) c1 = c;
          }   

          int c2 = downwind_cells_[f][0];
          for (int n = 1; n < downwind_cells_[f].size(); ++n) {
            if (downwind_cells_[f][n] == c) c2 = c;
          }

          u1 = (*field_)[component_][c];
          u2 = (*field_)[component_][c1 + c2 - c];
        }

        umin = std::min(u1, u2);
        umax = std::max(u1, u2);

        const AmanziGeometry::Point& xcf = mesh_->face_centroid(f);
        u1f = getValue(gradient_c1, c, xcf);

        // check if umin <= u1f <= umax
        if (u1f < umin) {
          normal_new = xcf - xc;
          CalculateDescentDirection_(normals, normal_new, L22normal_new, direction);

          // p = ((umin - u1) / sqrt(L22normal_new)) * direction;
          p = ((umin - u1) / sqrt(L22normal_new)) * normal_new;
          ApplyDirectionalLimiter_(normal_new, p, direction, gradient_c1);

        } else if (u1f > umax) {
          normal_new = xcf - xc;
          CalculateDescentDirection_(normals, normal_new, L22normal_new, direction);

          // p = ((umax - u1) / sqrt(L22normal_new)) * direction;
          p = ((umax - u1) / sqrt(L22normal_new)) * normal_new;
          ApplyDirectionalLimiter_(normal_new, p, direction, gradient_c1);
        }
      }
      if (normals.size() == 0) break;  // No limiters were imposed.
    }

    double grad_norm = norm(gradient_c1);
    if (grad_norm < OPERATOR_LIMITER_TOLERANCE * bc_scaling_) gradient_c1.set(0.0);

    for (int i = 0; i < dim; i++) grad[i][c] = gradient_c1[i];
  }

  // Local extrema are calculated here and updated in Step 2.
  AmanziMesh::Entity_ID_List cells;
  std::vector<double> field_local_min(ncells_wghost);
  std::vector<double> field_local_max(ncells_wghost);

  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_face_adj_cells(c, AmanziMesh::Parallel_type::ALL, &cells);
    field_local_min[c] = field_local_max[c] = (*field_)[component_][c];

    for (int i = 0; i < cells.size(); i++) {
      double value = (*field_)[component_][cells[i]];
      field_local_min[c] = std::min(field_local_min[c], value);
      field_local_max[c] = std::max(field_local_max[c], value);
    }
  }

  // Step 2: limit gradient on the Dirichlet boundary
  for (int f = 0; f < nfaces_owned; f++) {
    if (bc_model[f] == OPERATOR_BC_DIRICHLET) {
      int c1 = (upwind_cells_[f].size() > 0) ? upwind_cells_[f][0] : -1;
      int c2 = (downwind_cells_[f].size() > 0) ? downwind_cells_[f][0] : -1;

      if (c2 >= 0 && c2 < ncells_owned) {
        u2 = (*field_)[component_][c2];
        u1 = bc_value[f];
        umin = std::min(u1, u2);
        umax = std::max(u1, u2);

        bc_scaling_ = std::max(bc_scaling_, u1);
        field_local_max[c2] = std::max(field_local_max[c2], u1);
        field_local_min[c2] = std::min(field_local_min[c2], u1);

        const AmanziGeometry::Point& xc2 = mesh_->cell_centroid(c2);
        const AmanziGeometry::Point& xcf = mesh_->face_centroid(f);
        u2f = getValue(c2, xcf);
        for (int k = 0; k < dim; k++) gradient_c2[k] = grad[k][c2];
        direction = xcf - xc2;

        if (u2f < umin) {
          p = ((umin - u2) / L22(direction)) * direction;
          ApplyDirectionalLimiter_(direction, p, direction, gradient_c2);

        } else if (u2f > umax) {
          p = ((umax - u2) / L22(direction)) * direction;
          ApplyDirectionalLimiter_(direction, p, direction, gradient_c2);
        }

        for (int k = 0; k < dim; k++) grad[k][c2] = gradient_c2[k];
      }
    }
  }

  // Step 3: enforce a priori time step estimate (division of dT by 2).
  if (limiter_correction_) {
    LimiterExtensionTransportTensorial_(field_local_min, field_local_max);
  }    

  gradient_->ScatterMasterToGhosted("cell");

  // approximate estimate of scalar limiter (mainly for statistics)
  for (int c = 0; c < ncells_owned; c++) {
    double grad_norm0 = (*limiter_)[c];
    if (grad_norm0 == 0.0) { 
      (*limiter_)[c] = 1.0;
    } else {
      for (int i = 0; i < dim; i++) gradient_c1[i] = grad[i][c];
      (*limiter_)[c] = norm(gradient_c1) / grad_norm0;
    }
  }
}


/* *******************************************************************
* Extension of the tensorial limiter. Routine changes gradient to 
* satisfy an a prioty estimate of the stable time step. That estimate
* assumes that the weigthed flux is smaller than the first-order flux.
******************************************************************* */
void ReconstructionCell::LimiterExtensionTransportTensorial_(
    const std::vector<double>& field_local_min, const std::vector<double>& field_local_max)
{
  double u1f, u1;
  AmanziMesh::Entity_ID_List faces;

  Epetra_MultiVector& grad = *gradient_->ViewComponent("cell", false);

  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    double a, b;
    double flux, outflux = 0.0, outflux_weigted = 0.0;
    for (int i = 0; i < nfaces; i++) {
      int f = faces[i];
      int c1 = (upwind_cells_[f].size() > 0) ? upwind_cells_[f][0] : -1;

      if (c == c1) {
        const AmanziGeometry::Point& xcf = mesh_->face_centroid(f);
        u1f = getValue(c, xcf);
        u1 = (*field_)[component_][c];

        a = u1f - u1;
        if (fabs(a) > OPERATOR_LIMITER_TOLERANCE * (fabs(u1f) + fabs(u1))) {
          if (a > 0) b = u1 - field_local_min[c];
          else       b = u1 - field_local_max[c];

          flux = fabs((*flux_)[0][f]);
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

    if (outflux_weigted > outflux) {
      double psi = outflux / outflux_weigted;
      for (int i = 0; i < dim; i++) grad[i][c] *= psi;
    }
  }
}


/* *******************************************************************
* Routine calculates BJ limiter adjusted to linear advection. 
* First, it limits face-centered value of a reconstracted function 
* by min-max of two cell-centered values.
* Second, it limits outflux values which gives factor 0.5 in the
* time step estimate.
******************************************************************* */
void ReconstructionCell::LimiterBarthJespersen_(
    const std::vector<int>& bc_model, const std::vector<double>& bc_value,
    Teuchos::RCP<Epetra_Vector> limiter)
{
  limiter->PutScalar(1.0);
  Teuchos::RCP<Epetra_MultiVector> grad = gradient_->ViewComponent("cell", false);

  double u1, u2, u1f, u2f, umin, umax;  // cell and inteface values
  AmanziGeometry::Point gradient_c1(dim), gradient_c2(dim);

  // Step 1: limiting gradient inside domain
  for (int f = 0; f < nfaces_owned; f++) {
    int c1 = (upwind_cells_[f].size() > 0) ? upwind_cells_[f][0] : -1;
    int c2 = (downwind_cells_[f].size() > 0) ? downwind_cells_[f][0] : -1;
    if (c1 < 0 || c2 < 0) continue;

    u1 = (*field_)[component_][c1];
    u2 = (*field_)[component_][c2];
    umin = std::min(u1, u2);
    umax = std::max(u1, u2);
    double tol = sqrt(OPERATOR_LIMITER_TOLERANCE) * (fabs(u1) + fabs(u2));

    const AmanziGeometry::Point& xc1 = mesh_->cell_centroid(c1);
    const AmanziGeometry::Point& xc2 = mesh_->cell_centroid(c2);
    const AmanziGeometry::Point& xcf = mesh_->face_centroid(f);

    for (int i = 0; i < dim; i++) gradient_c1[i] = (*grad)[i][c1];
    double u1_add = gradient_c1 * (xcf - xc1);
    u1f = u1 + u1_add;

    if (u1f < umin - tol) {
      (*limiter)[c1] = std::min((*limiter)[c1], (umin - u1) / u1_add);
    } else if (u1f > umax + tol) {
      (*limiter)[c1] = std::min((*limiter)[c1], (umax - u1) / u1_add);
    }

    for (int i = 0; i < dim; i++) gradient_c2[i] = (*grad)[i][c2];
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
  std::vector<double> field_local_min(ncells_wghost);
  std::vector<double> field_local_max(ncells_wghost);

  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_face_adj_cells(c, AmanziMesh::Parallel_type::ALL, &cells);
    field_local_min[c] = field_local_max[c] = (*field_)[component_][c];

    for (int i = 0; i < cells.size(); i++) {
      double value = (*field_)[component_][cells[i]];
      field_local_min[c] = std::min(field_local_min[c], value);
      field_local_max[c] = std::max(field_local_max[c], value);
    }
  }

  // Step 2: limiting gradient on the Dirichlet boundary
  for (int f = 0; f < nfaces_owned; ++f) {
    if (bc_model[f] == OPERATOR_BC_DIRICHLET) {
      if (downwind_cells_[f].size() > 0) {
        int c2 = downwind_cells_[f][0];

        u2 = (*field_)[component_][c2];
        u1 = bc_value[f];
        umin = std::min(u1, u2);
        umax = std::max(u1, u2);
        double tol = sqrt(OPERATOR_LIMITER_TOLERANCE) * (fabs(u1) + fabs(u2));

        field_local_max[c2] = std::max(field_local_max[c2], u1);
        field_local_min[c2] = std::min(field_local_min[c2], u1);

        const AmanziGeometry::Point& xc2 = mesh_->cell_centroid(c2);
        const AmanziGeometry::Point& xcf = mesh_->face_centroid(f);

        for (int k = 0; k < dim; k++) gradient_c2[k] = (*grad)[k][c2];
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

  // Step 3: enforcing a priori time step estimate (division of dT by 2).
  if (limiter_correction_) {
    LimiterExtensionTransportBarthJespersen_(field_local_min, field_local_max, limiter);
  }    
}


/* *******************************************************************
* Extension of Barth-JEspersen's limiter. Routine changes gradient to 
* satisfy an a prioty estimate of the stable time step. That estimate
* assumes that the weigthed flux is smaller that the first-order flux.
******************************************************************* */
void ReconstructionCell::LimiterExtensionTransportBarthJespersen_(
    const std::vector<double>& field_local_min, const std::vector<double>& field_local_max,
    Teuchos::RCP<Epetra_Vector> limiter)
{
  double u1, u1f;
  AmanziGeometry::Point gradient_c1(dim);
  AmanziMesh::Entity_ID_List faces;

  Teuchos::RCP<Epetra_MultiVector> grad = gradient_->ViewComponent("cell", false);

  for (int c = 0; c < ncells_owned; c++) {
    const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    double a, b;
    double flux, outflux = 0.0, outflux_weigted = 0.0;
    for (int i = 0; i < nfaces; i++) {
      int f = faces[i];
      int c1 = (upwind_cells_[f].size() > 0) ? upwind_cells_[f][0] : -1;

      if (c == c1) {
        const AmanziGeometry::Point& xcf = mesh_->face_centroid(f);
        u1 = (*field_)[component_][c];
        for (int i = 0; i < dim; i++) gradient_c1[i] = (*grad)[i][c1];
        u1f = u1 + gradient_c1 * (xcf - xc);

        a = u1f - u1;
        if (fabs(a) > OPERATOR_LIMITER_TOLERANCE * (fabs(u1f) + fabs(u1))) {
          if (a > 0) b = u1 - field_local_min[c];
          else       b = u1 - field_local_max[c];

          flux = fabs((*flux_)[0][f]);
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
}


/* *******************************************************************
* Kuzmin's limiter use all neighbors of a computational cell.  
******************************************************************* */
void ReconstructionCell::LimiterKuzmin_(
    const std::vector<int>& bc_model, const std::vector<double>& bc_value)
{
  Teuchos::RCP<Epetra_MultiVector> grad = gradient_->ViewComponent("cell", false);

  // Step 1: local extrema are calculated here at nodes and updated later
  std::vector<double> field_node_min(nnodes_wghost);
  std::vector<double> field_node_max(nnodes_wghost);

  AmanziMesh::Entity_ID_List nodes;

  field_node_min.assign(nnodes_wghost,  OPERATOR_LIMITER_INFINITY);
  field_node_max.assign(nnodes_wghost, -OPERATOR_LIMITER_INFINITY);

  for (int c = 0; c < ncells_wghost; c++) {
    mesh_->cell_get_nodes(c, &nodes);

    double value = (*field_)[component_][c];
    for (int i = 0; i < nodes.size(); i++) {
      int v = nodes[i];
      field_node_min[v] = std::min(field_node_min[v], value);
      field_node_max[v] = std::max(field_node_max[v], value);
    }
  }

  // Update min/max at nodes from influx boundary data
  for (int f = 0; f < nfaces_owned; ++f) {
    int c1 = (upwind_cells_[f].size() > 0) ? upwind_cells_[f][0] : -1;
    int c2 = (downwind_cells_[f].size() > 0) ? downwind_cells_[f][0] : -1;

    if (c2 >= 0 && c1 < 0) {
      mesh_->face_get_nodes(f, &nodes);
      int nnodes = nodes.size();

      for (int i = 0; i < nnodes; i++) {
        int v = nodes[i];
        if (bc_model[v] == OPERATOR_BC_DIRICHLET) {
          double value = bc_value[v];
          field_node_min[v] = std::min(field_node_min[v], value);
          field_node_max[v] = std::max(field_node_max[v], value);
        }
      }
    }
  }

  // Step 2: limit reconstructed gradients at cell nodes
  double umin, umax, up, u1;
  AmanziGeometry::Point xp(dim);

  double L22normal_new;
  AmanziGeometry::Point gradient_c(dim), p(dim), normal_new(dim), direction(dim);
  std::vector<double> field_min_cell(cell_max_nodes), field_max_cell(cell_max_nodes);

  std::vector<AmanziGeometry::Point> normals;

  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_nodes(c, &nodes);
    int nnodes = nodes.size();
    for (int i = 0; i < nnodes; i++) {
      int v = nodes[i];
      field_min_cell[i] = field_node_min[v];
      field_max_cell[i] = field_node_max[v];
    }

    for (int i = 0; i < dim; i++) gradient_c[i] = (*grad)[i][c];

    LimiterKuzminCell_(c, gradient_c, field_min_cell, field_max_cell);

    for (int i = 0; i < dim; i++) (*grad)[i][c] = gradient_c[i];
  }

  if (limiter_correction_) {
    // Step 3: extrema are calculated for cells.
    AmanziMesh::Entity_ID_List cells;
    std::vector<double> field_local_min(ncells_wghost);
    std::vector<double> field_local_max(ncells_wghost);

    for (int c = 0; c < ncells_owned; c++) {
      mesh_->cell_get_nodes(c, &nodes);
      field_local_min[c] = field_local_max[c] = (*field_)[component_][c];

      for (int i = 0; i < nodes.size(); i++) {
        int v = nodes[i];
        field_local_min[c] = std::min(field_local_min[c], field_node_min[v]);
        field_local_max[c] = std::max(field_local_max[c], field_node_max[v]);
      }
    }

    // Step 4: enforcing a priori time step estimate (division of dT by 2).
    LimiterExtensionTransportKuzmin_(field_local_min, field_local_max);
  }    

  gradient_->ScatterMasterToGhosted("cell");
}



/* *******************************************************************
* Kuzmin's limiter over a subset of cells.  
******************************************************************* */
void ReconstructionCell::LimiterKuzminSet_(AmanziMesh::Entity_ID_List& ids,
                                           std::vector<AmanziGeometry::Point>& gradient)
{
  // Step 1: local extrema are calculated here at nodes and updated later
  std::vector<double> field_node_min(nnodes_wghost);
  std::vector<double> field_node_max(nnodes_wghost);
  std::vector<double> field_min_cell(cell_max_nodes), field_max_cell(cell_max_nodes);

  AmanziMesh::Entity_ID_List nodes;

  field_node_min.assign(nnodes_wghost,  OPERATOR_LIMITER_INFINITY);
  field_node_max.assign(nnodes_wghost, -OPERATOR_LIMITER_INFINITY);
  bc_scaling_ = 1.0;

  for (int c = 0; c < ncells_wghost; ++c) {
    mesh_->cell_get_nodes(c, &nodes);
    double value = (*field_)[component_][c];
    for (int i = 0; i < nodes.size(); ++i) {
      int v = nodes[i];
      field_node_min[v] = std::min(field_node_min[v], value);
      field_node_max[v] = std::max(field_node_max[v], value);
    }
  }

  for (int k = 0; k < ids.size(); ++k){
    int c = ids[k];
    mesh_->cell_get_nodes(c, &nodes);
    int nnodes = nodes.size();
    for (int i = 0; i < nnodes; i++) {
      int v = nodes[i];
      field_min_cell[i] = field_node_min[v];
      field_max_cell[i] = field_node_max[v];
    }
    LimiterKuzminCell_(c, gradient[k], field_min_cell, field_max_cell);
  }
}


/* *******************************************************************
* Kuzmin's limiter use all neighbors of the given cell and limit 
* gradient in this cell only.  
******************************************************************* */
void ReconstructionCell::LimiterKuzminCell_(int cell,
                                            AmanziGeometry::Point& gradient_c,
                                            const std::vector<double>& field_node_min_c,
                                            const std::vector<double>& field_node_max_c)
{
  double umin, umax, up, u1;
  AmanziGeometry::Point xp(dim);

  double L22normal_new;
  AmanziGeometry::Point p(dim), normal_new(dim), direction(dim);
  AmanziMesh::Entity_ID_List nodes;
  std::vector<AmanziGeometry::Point> normals;

  mesh_->cell_get_nodes(cell, &nodes);
  int nnodes = nodes.size();

  u1 = (*field_)[component_][cell];

  const AmanziGeometry::Point& xc = mesh_->cell_centroid(cell);

  normals.clear();  // normals to planes the define the feasiable set
  for (int loop = 0; loop < 2; loop++) {
    for (int i = 0; i < nnodes; i++) {
      int v = nodes[i];
      double umin = field_node_min_c[i];
      double umax = field_node_max_c[i];

      mesh_->node_get_coordinates(v, &xp);
      up = getValue(gradient_c, cell, xp);

      // check if umin <= up <= umax
      if (up < umin) {
        normal_new = xp - xc;
        CalculateDescentDirection_(normals, normal_new, L22normal_new, direction);

        // p = ((umin - u1) / sqrt(L22normal_new)) * direction;
        p = ((umin - u1) / sqrt(L22normal_new)) * normal_new;
        ApplyDirectionalLimiter_(normal_new, p, direction, gradient_c);
      } else if (up > umax) {
        normal_new = xp - xc;
        CalculateDescentDirection_(normals, normal_new, L22normal_new, direction);

        // p = ((umax - u1) / sqrt(L22normal_new)) * direction;
        p = ((umax - u1) / sqrt(L22normal_new)) * normal_new;
        ApplyDirectionalLimiter_(normal_new, p, direction, gradient_c);
      }
    }
    if (normals.size() == 0) break;  // No limiters were imposed.
  
    double grad_norm = norm(gradient_c);
    if (grad_norm < OPERATOR_LIMITER_TOLERANCE * bc_scaling_) gradient_c.set(0.0);
  }
}


/* *******************************************************************
* Extension of Kuzmin's limiter. Routine changes gradient to 
* satisfy an a prioty estimate of the stable time step. That estimate
* assumes that the weigthed flux is smaller that the first-order flux.
******************************************************************* */
void ReconstructionCell::LimiterExtensionTransportKuzmin_(
    const std::vector<double>& field_local_min, const std::vector<double>& field_local_max)
{
  double u1, up;
  AmanziGeometry::Point xp(dim);
  AmanziMesh::Entity_ID_List faces, nodes;

  Teuchos::RCP<Epetra_MultiVector> grad = gradient_->ViewComponent("cell", false);

  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    double a, b;
    double flux, outflux = 0.0, outflux_weigted = 0.0;
    for (int i = 0; i < nfaces; i++) {
      int f = faces[i];
      int c1 = (upwind_cells_[f].size() > 0) ? upwind_cells_[f][0] : -1;

      if (c == c1) {
        mesh_->face_get_nodes(f, &nodes);
        int nnodes = nodes.size();

        // define dimensionless quadrature weigths
        std::vector<double> weights(nnodes, 0.5);
        if (dim == 3) {
          double area = mesh_->face_area(f);
          WhetStone::PolygonCentroidWeights(*mesh_, nodes, area, weights);
        }

        for (int j = 0; j < nnodes; j++) {
          int v = nodes[j];
          mesh_->node_get_coordinates(v, &xp);
          up = getValue(c, xp);
          u1 = (*field_)[component_][c];

          a = up - u1;
          if (fabs(a) > OPERATOR_LIMITER_TOLERANCE * (fabs(up) + fabs(u1))) {
            if (a > 0)
              b = u1 - field_local_min[c];
            else
              b = u1 - field_local_max[c];

            flux = fabs((*flux_)[0][f]) * weights[j];
            outflux += flux;
            if (b) {
              outflux_weigted += flux * a / b;
            } else {
              for (int k = 0; k < dim; k++) (*grad)[k][c] = 0.0;
              break;
            }
          }
        }
      }
    }

    if (outflux_weigted > outflux) {
      double psi = outflux / outflux_weigted;
      for (int i = 0; i < dim; i++) (*grad)[i][c] *= psi;
    }
  }
}


/* *******************************************************************
* Descent direction is obtained by orthogonalizing normal direction
* 'normal_new' to previous normals. A few exceptions are analyzed.  
******************************************************************* */
void ReconstructionCell::CalculateDescentDirection_(
    std::vector<AmanziGeometry::Point>& normals,
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
    if (a < L22normal_new * OPERATOR_LIMITER_TOLERANCE) {
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
******************************************************************* */
void ReconstructionCell::ApplyDirectionalLimiter_(AmanziGeometry::Point& normal,
                                                  AmanziGeometry::Point& p,
                                                  AmanziGeometry::Point& direction,
                                                  AmanziGeometry::Point& gradient)
{
  double a = ((p - gradient) * normal) / (direction * normal);
  gradient += a * direction;
}


/* *******************************************************************
* Identify flux direction based on orientation of the face normal 
* and sign of the Darcy velocity.                               
******************************************************************* */
void ReconstructionCell::IdentifyUpwindCells_()
{
  upwind_cells_.clear();
  downwind_cells_.clear();

  upwind_cells_.resize(nfaces_wghost);
  downwind_cells_.resize(nfaces_wghost);

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  for (int c = 0; c < ncells_wghost; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);

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

}  // namespace Operators
}  // namespace Amanzi

