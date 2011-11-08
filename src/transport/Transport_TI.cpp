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

  calculateLimiterBarthJespersen(current_component_, 
                                 component_rcp, 
                                 gradient, 
                                 field_local_min, 
                                 field_local_max, 
                                 limiter_);
  lifting.applyLimiter(limiter_);

  // ADVECTIVE FLUXES
  int f, c1, c2;
  double u, upwind_tcc, tcc_flux;

  f_component.PutScalar(0.0);
  for (int f=fmin; f<=fmax; f++) {  // loop over master and slave faces
    c1 = (*upwind_cell_)[f]; 
    c2 = (*downwind_cell_)[f]; 

    u = fabs(darcy_flux[f]);
    const AmanziGeometry::Point& xf = mesh_->face_centroid(f);

    if (c1 >=0 && c1 <= cmax_owned && c2 >= 0 && c2 <= cmax_owned) {
      upwind_tcc = lifting.getValue(c1, xf);
      tcc_flux = u * upwind_tcc;
      f_component[c1] -= tcc_flux;
      f_component[c2] += tcc_flux;
    } 
    else if (c1 >=0 && c1 <= cmax_owned && (c2 > cmax_owned || c2 < 0)) {
      upwind_tcc = lifting.getValue(c1, xf);
      tcc_flux = u * upwind_tcc;
      f_component[c1] -= tcc_flux;
    } 
    else if (c1 > cmax_owned && c2 >= 0 && c2 <= cmax_owned) {
      upwind_tcc = lifting.getValue(c1, xf);
      tcc_flux = u * upwind_tcc;
      f_component[c2] += tcc_flux;
    }
  } 

  for (int c=cmin; c<=cmax_owned; c++) {  // calculate conservative quantatity 
    double vol_phi_ws = mesh_->cell_volume(c) * phi[c] * ws[c]; 
    f_component[c] /= vol_phi_ws; 
  }

  // BOUNDARY CONDITIONS for ADVECTION
  for (int n=0; n<bcs.size(); n++) {  // analyze boundary sets
    for (int k=0; k<bcs[n].faces.size(); k++) {
      f = bcs[n].faces[k];
      c2 = (*downwind_cell_)[f]; 

      if (c2 >= 0) {
        u = fabs(darcy_flux[f]);

        if (bcs[n].type == TRANSPORT_BC_CONSTANT_TCC) {
          double vol_phi_ws = mesh_->cell_volume(c2) * phi[c2] * ws[c2]; 
          tcc_flux = u * bcs[n].values[current_component_];
          f_component[c2] += tcc_flux / vol_phi_ws;
          bcs[n].influx[current_component_] += tcc_flux;
        }
      } 
    }
  }   
}

}  // namespace AmanziTransport
}  // namespace Amanzi


