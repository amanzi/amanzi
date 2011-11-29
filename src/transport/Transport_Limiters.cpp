/*
This is the transport component of the Amanzi code. 
License: BSD
Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "Transport_PK.hpp"

namespace Amanzi {
namespace AmanziTransport {

/* *******************************************************************
 * Tensorial limiter limits the gradient directly, to avoid 
 * calculation of a 3x3 matrix.
 ****************************************************************** */
void Transport_PK::limiterTensorial(const int component,
                                    Teuchos::RCP<Epetra_Vector> scalar_field, 
                                    Teuchos::RCP<Epetra_MultiVector> gradient)
{
  const Epetra_Vector& darcy_flux = TS_nextBIG->ref_darcy_flux();

  double u1, u2, u1f, u2f, umin, umax, L22normal_new;
  AmanziGeometry::Point gradient_c1(dim), gradient_c2(dim);
  AmanziGeometry::Point normal_new(dim), direction(dim), p(dim);

  std::vector<AmanziGeometry::Point> normals;
  AmanziMesh::Entity_ID_List faces;

  for (int c=cmin; c<=cmax_owned; c++) {
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
    for (int i=0; i<dim; i++) gradient_c1[i] = (*gradient)[i][c]; 
 
    // project along direction to the feasible set 
    normals.clear();
    for (int loop=0; loop<2; loop++) {
      for (int i=0; i<nfaces; i++) {
        int f = faces[i];
        int c1 = (*upwind_cell_)[f];
        int c2 = (*downwind_cell_)[f];

        if (c1>=0 && c2>=0) {
          u1 = (*scalar_field)[c1];
          u2 = (*scalar_field)[c2];
        } else if (c1==c) {
          u1 = u2 = (*scalar_field)[c1];
        } else {  // This case is analyze by influx b.c.
          continue;
        }
        umin = std::min(u1, u2);
        umax = std::max(u1, u2);

        const AmanziGeometry::Point& xcf = mesh_->face_centroid(f);
        u1f = lifting.getValue(gradient_c1, c, xcf);

        // check if umin <= u1f <= umax
        if (u1f < umin) {
          normal_new = xcf - xc;
          calculate_descent_direction(normals, normal_new, L22normal_new, direction);
 
          p = ((umin - u1) / sqrt(L22normal_new)) * direction;
          apply_directional_limiter(normal_new, p, direction, gradient_c1);
        } else if (u1f > umax) {
          normal_new = xcf - xc;
          calculate_descent_direction(normals, normal_new, L22normal_new, direction);
 
          p = ((umax - u1) / sqrt(L22normal_new)) * direction;
          apply_directional_limiter(normal_new, p, direction, gradient_c1);
        }
      }
      if (normals.size() == 0) break;  // No limiters were imposed. 
    }

    double grad_norm = norm(gradient_c1);
    if (grad_norm < TRANSPORT_LIMITER_TOLERANCE * bc_scaling) gradient_c1 = 0.0;

    for (int i=0; i<dim; i++) (*gradient)[i][c] = gradient_c1[i];
  }
 
  // limiting gradient on the Dirichlet boundary
  std::vector<double> influx(number_wghost_cells, 0.0);

  for (int n=0; n<bcs.size(); n++) {
    if (component == bcs_tcc_index[n]) {
      for (BoundaryFunction::Iterator bc=bcs[n]->begin(); bc != bcs[n]->end(); ++bc) {
        int f = bc->first;
        int c1 = (*upwind_cell_)[f]; 
        int c2 = (*downwind_cell_)[f]; 

        if (c2 >= 0) {
          u2 = (*scalar_field)[c2];
          u1 = bc->second;
          umin = std::min(u1, u2);
          umax = std::max(u1, u2);

          bc_scaling = std::max(bc_scaling, u1);

          const AmanziGeometry::Point& xc2 = mesh_->cell_centroid(c2);
          const AmanziGeometry::Point& xcf = mesh_->face_centroid(f);
          u2f = lifting.getValue(c2, xcf);
          for (int i=0; i<dim; i++) gradient_c2[i] = (*gradient)[i][c2]; 
          direction = xcf - xc2;
 
          if (u2f < umin) {
            p = ((umin - u2) / L22(direction)) * direction;
            apply_directional_limiter(direction, p, direction, gradient_c2);
          }
          else if (u2f > umax) {
            p = ((umax - u2) / L22(direction)) * direction;
            apply_directional_limiter(direction, p, direction, gradient_c2);
          }

          for (int i=0; i<dim; i++) (*gradient)[i][c2] = gradient_c2[i]; 
          influx[c2] += fabs(u1);
        } 
      }
    }
  }
 
  // cleaning local extrema (experimental lipnikov@lanl.gov)
  //for (int c=cmin; c<=cmax_owned; c++) {
  //  u1 = (*scalar_field)[c];
  //   if (u1 <= lifting.get_field_local_min()[c] || u1 >= lifting.get_field_local_max()[c]) {
  //    for (int i=0; i<dim; i++) (*gradient)[i][c] = 0.0;
  //  } 
  //}

  TS_nextBIG->distribute_cell_multivector(*gradient);
}


/* *******************************************************************
 * Routine calculates BJ limiter adjusted to linear advection. 
 * First, it limits face-centered value of a reconstracted function 
 * by min-max of two cell-centered values.
 * Second, it limits outflux values which gives factor 0.5 in the
 * time step estimate.
 ****************************************************************** */
void Transport_PK::limiterBarthJespersen(const int component,
                                         Teuchos::RCP<Epetra_Vector> scalar_field, 
                                         Teuchos::RCP<Epetra_MultiVector> gradient,
                                         std::vector<double>& field_local_min,
                                         std::vector<double>& field_local_max,
                                         Teuchos::RCP<Epetra_Vector> limiter)
{
  for (int c=cmin; c<=cmax; c++) (*limiter)[c] = 1.0;
 
  std::vector<double> total_influx(number_wghost_cells, 0.0);  // follows calculation of dT
  std::vector<double> total_outflux(number_wghost_cells, 0.0);
  const Epetra_Vector& darcy_flux = TS_nextBIG->ref_darcy_flux();

  double u1, u2, u1f, u2f, umin, umax;  // cell and inteface values
  AmanziGeometry::Point gradient_c1(dim), gradient_c2(dim);

  for (int f=fmin; f<=fmax_owned; f++) {
    int c1, c2;
    c1 = (*upwind_cell_)[f];
    c2 = (*downwind_cell_)[f];
    if (c1 < 0 || c2 < 0) return;  // the desing should be modified (lipnikov@lanl.gov)

    const AmanziGeometry::Point& xc1 = mesh_->cell_centroid(c1);
    const AmanziGeometry::Point& xc2 = mesh_->cell_centroid(c2);
    const AmanziGeometry::Point& xcf = mesh_->face_centroid(f);
 
    for (int i=0; i<dim; i++) gradient_c1[i] = (*gradient)[i][c1]; 
    double u1_add = gradient_c1 * (xcf - xc1);
    u1f = u1 + u1_add;

    if (u1f < umin && u1_add) {
      (*limiter)[c1] = std::min((*limiter)[c1], (umin - u1) / u1_add);
    }
    else if (u1f > umax && u1_add) {
      (*limiter)[c1] = std::min((*limiter)[c1], (umax - u1) / u1_add); 
    }

    for (int i=0; i<dim; i++) gradient_c2[i] = (*gradient)[i][c2]; 
    double u2_add = gradient_c2 * (xcf - xc2);
    u2f = u2 + u2_add;

    if (u2f < umin && u2_add) {
      (*limiter)[c2] = std::min((*limiter)[c2], (umin - u2) / u2_add);
    }
    else if (u2f > umax && u2_add) {
      (*limiter)[c2] = std::min((*limiter)[c2], (umax - u2) / u2_add); 
    }

    total_influx[c2] += fabs(darcy_flux[f]);
    total_outflux[c1] += fabs(darcy_flux[f]); 
  } 
 
  // limiting second term in the flux formula (outflow darcy_flux)
  // has to be moved down, after BC. (lipnikov@lan.gov)
  for (int c=cmin; c<=cmax_owned; c++) {
    if (total_outflux[c]) {
      double psi = std::min(total_influx[c] / total_outflux[c], 1.0);
      (*limiter)[c] = std::min((*limiter)[c], psi);
    }
  }

  // limiting gradient on the Dirichlet boundary
  for (int n=0; n<bcs.size(); n++) {
    if (component == bcs_tcc_index[n]) {
      for (BoundaryFunction::Iterator bc=bcs[n]->begin(); bc != bcs[n]->end(); ++bc) {
        int f = bc->first;
        int c1 = (*upwind_cell_)[f]; 
        int c2 = (*downwind_cell_)[f]; 

        if (c2 >= 0) {
          u2 = (*scalar_field)[c2];
          u1 = bc->second;
          umin = std::min(u1, u2);
          umax = std::max(u1, u2);

          const AmanziGeometry::Point& xc2 = mesh_->cell_centroid(c2);
          const AmanziGeometry::Point& xcf = mesh_->face_centroid(f);

          for (int i=0; i<dim; i++) gradient_c2[i] = (*gradient)[i][c2]; 
          double u_add = gradient_c2 * (xcf - xc2);
          u2f = u2 + u_add;

          if (u2f < umin && u_add) {
            (*limiter)[c2] = std::min((*limiter)[c2], (umin - u2) / u_add);
          }
          else if (u2f > umax && u_add) {
            (*limiter)[c2] = std::min((*limiter)[c2], (umax - u2) / u_add); 
          }
        } 
        else if (c1 >= 0) {
          //(*limiter)[c1] = 0.0;  // ad-hoc testing (lipnikov@lanl.gov)
        }
      }
    }
  }

  TS_nextBIG->distribute_cell_vector(*limiter_);
}
 

/* *******************************************************************
 * Descent direction is obtained by orthogonalizing normal direction
 * 'normal_new' to previous normals. A few exceptions are analyzed.  
 ****************************************************************** */
void Transport_PK::calculate_descent_direction(std::vector<AmanziGeometry::Point>& normals,
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
    for (int n=0; n<nnormals; n++) {
      a = normals[n] * direction;
      for (int i=0; i<dim; i++) direction[i] -= a * normals[n][i];
    }

    // verify new direction
    a = L22(direction);
    if (a < L22normal_new * TRANSPORT_LIMITER_TOLERANCE) {
      normals.clear();
      direction = normal_new;
    } else {
      for (int i=0; i<dim; i++) direction[i] /= sqrt(a);
    }
  }
  normals.push_back(direction); 
}


/* *******************************************************************
 * Routine projects gradient on a plane defined by normal and point p.
 ****************************************************************** */
void Transport_PK::apply_directional_limiter(AmanziGeometry::Point& normal, 
                                             AmanziGeometry::Point& p,
                                             AmanziGeometry::Point& direction, 
                                             AmanziGeometry::Point& gradient)
{
  double a = ((p - gradient) * normal) / (direction * normal); 
  gradient += a * direction;
}
  
}  // namespace AmanziTransport
}  // namespace Amanzi

