/*
This is the transport component of the Amanzi code. 
License: BSD
Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "Transport_PK.hpp"

namespace Amanzi {
namespace AmanziTransport {

void Transport_PK::fun(const double t, const Epetra_Vector& component, Epetra_Vector& f_component)
{
  const Epetra_Vector& darcy_flux = TS_nextBIG->ref_darcy_flux();
  const Epetra_Vector& ws  = TS_nextBIG->ref_water_saturation();
  const Epetra_Vector& phi = TS_nextBIG->ref_porosity();

  // transport routines need RCP pointer
  Epetra_Vector* component_tmp = const_cast<Epetra_Vector*> (&component);  
  const Teuchos::RCP<Epetra_Vector> component_rcp(component_tmp, false);

  lifting.reset_field(mesh_, component_rcp);
  lifting.calculateCellGradient();
 
  Teuchos::RCP<Epetra_MultiVector> gradient = lifting.get_gradient();
  std::vector<double>& field_local_min = lifting.get_field_local_min();
  std::vector<double>& field_local_max = lifting.get_field_local_max();

  if (advection_limiter == TRANSPORT_LIMITER_BARTH_JESPERSEN) {  
    limiterBarthJespersen(
        current_component_, component_rcp, gradient, field_local_min, field_local_max, limiter_);
    lifting.applyLimiter(limiter_);
  }
  else if (advection_limiter == TRANSPORT_LIMITER_TENSORIAL) {
    limiterTensorial(current_component_, component_rcp, gradient);
  }
cout << "passed " << MyPID << endl;

  // ADVECTIVE FLUXES
  // We assume that limiters made their job up to round-off errors. 
  // Min-max condition will enforce robustness for these errors.
  int f, c1, c2;
  double u, u1, u2, umin, umax, upwind_tcc, tcc_flux;

  f_component.PutScalar(0.0);
  for (int f=fmin; f<=fmax; f++) {  // loop over master and slave faces
    c1 = (*upwind_cell_)[f]; 
    c2 = (*downwind_cell_)[f]; 

    if (c1 >= 0 && c2 >= 0) {
      u1 = component[c1];
      u2 = component[c2];
      umin = std::min(u1, u2);
      umax = std::max(u1, u2);
    }
    else if (c1 >= 0) {
      u1 = u2 = umin = umax = component[c1];
    }
    else if (c2 >= 0) {
      u1 = u2 = umin = umax = component[c2];
    }

    u = fabs(darcy_flux[f]);
    const AmanziGeometry::Point& xf = mesh_->face_centroid(f);

    if (c1 >= 0 && c1 <= cmax_owned && c2 >= 0 && c2 <= cmax_owned) {
      upwind_tcc = lifting.getValue(c1, xf);
      upwind_tcc = std::max(upwind_tcc, umin);
      upwind_tcc = std::min(upwind_tcc, umax);

      tcc_flux = u * upwind_tcc;
      f_component[c1] -= tcc_flux;
      f_component[c2] += tcc_flux;
    } 
    else if (c1 >= 0 && c1 <= cmax_owned && (c2 > cmax_owned || c2 < 0)) {
      upwind_tcc = lifting.getValue(c1, xf);
      upwind_tcc = std::max(upwind_tcc, umin);
      upwind_tcc = std::min(upwind_tcc, umax);

      tcc_flux = u * upwind_tcc;
      f_component[c1] -= tcc_flux;
    } 
    else if (c1 > cmax_owned && c2 >= 0 && c2 <= cmax_owned) {
      upwind_tcc = lifting.getValue(c1, xf);
      upwind_tcc = std::max(upwind_tcc, umin);
      upwind_tcc = std::min(upwind_tcc, umax);

      tcc_flux = u * upwind_tcc;
      f_component[c2] += tcc_flux;
    }
  } 

  for (int c=cmin; c<=cmax_owned; c++) {  // calculate conservative quantatity 
    double vol_phi_ws = mesh_->cell_volume(c) * phi[c] * ws[c]; 
    f_component[c] /= vol_phi_ws; 
  }

  // BOUNDARY CONDITIONS for ADVECTION
  for (int n=0; n<bcs.size(); n++) {
    if (current_component_ == bcs_tcc_index[n]) {
      for (BoundaryFunction::Iterator bc=bcs[n]->begin(); bc != bcs[n]->end(); ++bc) {
        f = bc->first;
        c2 = (*downwind_cell_)[f]; 

        if (c2 >= 0) {
          u = fabs(darcy_flux[f]);
          double vol_phi_ws = mesh_->cell_volume(c2) * phi[c2] * ws[c2]; 
          tcc_flux = u * bc->second;
          f_component[c2] += tcc_flux / vol_phi_ws;
        }
      } 
    }
  }   
}

}  // namespace AmanziTransport
}  // namespace Amanzi


