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
  std::vector<double> total_influx(number_wghost_cells, 0.0);  // follows calculation of dT
  std::vector<double> total_outflux(number_wghost_cells, 0.0);
  const Epetra_Vector& darcy_flux = TS_nextBIG->ref_darcy_flux();

  double u1, u2, u1f, u2f, umin, umax;
  AmanziGeometry::Point gradient_c1(dim), gradient_c2(dim), direction(dim), p(dim);

  for (int c=cmin; c<=cmax_owned; c++) {
    AmanziMesh::Entity_ID_List faces;
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    const AmanziGeometry::Point& xc1 = mesh_->cell_centroid(c);
    for (int i=0; i<dim; i++) gradient_c1[i] = (*gradient)[i][c]; 

    // find the first violation of feasible set 
    bool flag;
    for (int loop=0; loop<2; loop++) {
      flag = true;
      for (int i=0; i<nfaces; i++) {
        int f = faces[i];
        int c1 = (*upwind_cell_)[f];
        int c2 = (*downwind_cell_)[f];
        if (c1<0 || c2<0) continue; // boundary faces are considered separately

        const AmanziGeometry::Point& xcf = mesh_->face_centroid(f);

        u1 = (*scalar_field)[c1];
        u2 = (*scalar_field)[c2];
        umin = std::min(u1, u2);
        umax = std::max(u1, u2);

        // check if umin <= u1f <= umax
        u1f = lifting.getValue(gradient_c1, c1, xcf);
        direction = xcf - xc1;

        if (u1f < umin) {
          p = ((umin - u1) / L22(direction)) * direction;
          apply_directional_limiter(direction, p, direction, gradient_c1);
          flag = false;
        }
        else if (u1f > umax) {
          p = ((umax - u1) / L22(direction)) * direction;
          apply_directional_limiter(direction, p, direction, gradient_c1);
          flag = false;
        }

        if (loop == 0) { 
          total_influx[c2] += fabs(darcy_flux[f]);
          total_outflux[c1] += fabs(darcy_flux[f]); 
        }
      }
      if (flag) break;  // No limiters were needed. 
    }
    for (int i=0; i<dim; i++) (*gradient)[i][c] = gradient_c1[i];
  }
 
  // limiting gradient on the Dirichlet boundary
  for ( int n=0; n<bcs.size(); n++) {
    for (int k=0; k<bcs[n].faces.size(); k++) {
      int f = bcs[n].faces[k];
      int c1 = (*upwind_cell_)[f]; 
      int c2 = (*downwind_cell_)[f]; 

      if (c2 >= 0) {
        u2 = (*scalar_field)[c2];

        if (bcs[n].type == TRANSPORT_BC_CONSTANT_TCC) {
          u1 = bcs[n].values[component];
          umin = std::min(u1, u2);
          umax = std::max(u1, u2);

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
          total_influx[c2] += fabs(u1);
        } 
      }
    }
  }

  // limiting second term in the flux formula (outflow darcy_flux)
  for (int c=cmin; c<=cmax_owned; c++) {
    if (total_outflux[c]) {
      double psi = std::min(total_influx[c] / total_outflux[c], 1.0);
      for (int i=0; i<dim; i++) (*gradient)[i][c] *= psi; 
    }
  }

#ifdef HAVE_MPI
  const Epetra_BlockMap& source_fmap = (*limiter_).Map();
  const Epetra_BlockMap& target_fmap = (*limiter_).Map();

  Epetra_Import importer(target_fmap, source_fmap);
  (*limiter_).Import(*limiter_, importer, Insert);
#endif
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
  for ( int n=0; n<bcs.size(); n++) {
    for (int k=0; k<bcs[n].faces.size(); k++) {
      int f = bcs[n].faces[k];
      int c1 = (*upwind_cell_)[f]; 
      int c2 = (*downwind_cell_)[f]; 

      if (c2 >= 0) {
        u2 = (*scalar_field)[c2];

        if (bcs[n].type == TRANSPORT_BC_CONSTANT_TCC) {
          u1 = bcs[n].values[component];
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
      }
      else if (c1 >= 0) {
        //(*limiter)[c1] = 0.0;  // ad-hoc testing (lipnikov@lanl.gov)
      }
    }
  }

#ifdef HAVE_MPI
  const Epetra_BlockMap& source_fmap = (*limiter_).Map();
  const Epetra_BlockMap& target_fmap = (*limiter_).Map();

  Epetra_Import importer(target_fmap, source_fmap);
  (*limiter_).Import(*limiter_, importer, Insert);
#endif
}
 

/* *******************************************************************
 * Routine projects gradient on a plane defined by normal and point p.
 ****************************************************************** */
void Transport_PK::apply_directional_limiter(AmanziGeometry::Point& normal, 
                                             AmanziGeometry::Point& p,
                                             AmanziGeometry::Point& direction, 
                                             AmanziGeometry::Point& gradient)
{
  double norm = L22(direction);
  double a = ((p - gradient) * direction) / norm;
 
  gradient += a * direction;
}
  
}  // namespace AmanziTransport
}  // namespace Amanzi

