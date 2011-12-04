/*
This is the transport component of the Amanzi code. 
License: BSD
Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "errors.hh"

#include "Transport_PK.hpp"

namespace Amanzi {
namespace AmanziTransport {

/* *******************************************************************
 * Routine verifies that the velocity field is divergence free                 
 ****************************************************************** */
void Transport_PK::check_divergence_property()
{
  int i, c, f;
  double div, u, umax, error_max, error_avg;

  Teuchos::RCP<AmanziMesh::Mesh> mesh = TS->get_mesh_maps();
  Epetra_Vector& darcy_flux = TS_nextBIG->ref_darcy_flux();

  error_max = error_avg = 0.0;

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> fdirs;

  for (int c=cmin; c<=cmax_owned; c++) {
    mesh->cell_get_faces(c, &faces);
    mesh->cell_get_face_dirs(c, &fdirs);
    int nfaces = faces.size();

    div = umax = 0;
    for (i=0; i<nfaces; i++) {
      f = faces[i];
      u = darcy_flux[f];
      div += u * fdirs[i];
      umax = std::max(umax, fabs(u) / pow(mesh->face_area(f), 0.5));
    }
    div /= mesh->cell_volume(c);

    if (umax) {
      error_max = std::max(error_max, fabs(div) / umax);
      error_avg += fabs(div) / umax;
    }

    /* verify that divergence complies with the flow model  */
    int flag = 0;
    if (TRANSPORT_AMANZI_VERSION == 1 && fabs(div) > tests_tolerance * umax) {
      cout << "TRANSPORT: The flow violates conservation property." << endl;
      cout << "    Modify either flow convergence criteria or transport tolerance." << endl;
      flag = 1;
    }
    if (TRANSPORT_AMANZI_VERSION == 2 && div < -tests_tolerance * umax) {
      cout << "TRANSPORT: The flow has large artificial sinks."<< endl;
      flag = 1;
    }

    if (flag && verbosity_level > 1) {
      cout << "    MyPID = " << MyPID << endl;
      cout << "    cell  = " << c << endl;
      cout << "    divergence = " << div << endl;
      cout << "    maximal velocity = " << umax << endl;
      Errors::Message msg;
      msg << "Velocity field is not divergence free " << "\n";
      Exceptions::amanzi_throw(msg);
    }
  }
  error_avg /= (cmax_owned + 1);

  if (verbosity_level > 3) {
#ifdef HAVE_MPI
    double global_max;
    const Epetra_Comm& comm = darcy_flux.Comm(); 
 
    comm.MinAll(&error_max, &global_max, 1);
    error_max = global_max;
#endif
    if (! MyPID) {
      cout << "Transport_PK: " << endl; 
      cout << "    maximal (divergence / flux) = " << error_max << endl;
      cout << "    average (divergence / flux) = " << error_avg << endl; 
    }
  }
}


/* *******************************************************************
 * Check for completeness of influx boundary conditions.                        
 ****************************************************************** */
void Transport_PK::check_influx_bc() const
{
  const Epetra_Vector& darcy_flux = TS_nextBIG->ref_darcy_flux();
  std::vector<int> influx_face(fmax + 1);

  for (int i=0; i<number_components; i++) {
    influx_face.assign(fmax + 1, 0);    

    for (int n=0; n<bcs.size(); n++) {
      if (i == bcs_tcc_index[n]) {
        for (BoundaryFunction::Iterator bc=bcs[n]->begin(); bc != bcs[n]->end(); ++bc) {
          int f = bc->first;
          influx_face[f] = 1; 
        }
      }
    }

    for (int f=0; f<=fmax_owned; f++) {
      if (darcy_flux[f] < 0 && influx_face[f] == 0) {
        char component[3];
        std::sprintf(component, "%3d", i);

        Errors::Message msg;
        msg << "No influx boundary condition has been found for component " << component << ".\n";
        Exceptions::amanzi_throw(msg);
      }
    }
  }
}


/* *******************************************************************
 * Check that global extrema diminished                          
 ****************************************************************** */
void Transport_PK::check_GEDproperty(Epetra_MultiVector& tracer) const
{ 
  int i, num_components = tracer.NumVectors();
  double tr_min[num_components];
  double tr_max[num_components];

  tracer.MinValue(tr_min);
  tracer.MaxValue(tr_max);

  if (TRANSPORT_AMANZI_VERSION == 1) {
    for (i=0; i<num_components; i++) {
      if (tr_min[i] < 0) {
        cout << "Transport_PK: concentration violates GED property" << endl; 
        cout << "    Make an Amanzi ticket or turn off internal transport tests" << endl;
        cout << "    MyPID = " << MyPID << endl;
        cout << "    component = " << i << endl;
        cout << "    time = " << T_internal << endl;
        cout << "    min/max values = " << tr_min[i] << " " << tr_max[i] << endl;

        Errors::Message msg;
        msg << "Concentration violates GED property." << "\n";
        Exceptions::amanzi_throw(msg);
      }
    } 
  }
}


/* *******************************************************************
 * Check that te tracer is between 0 and 1                         
 ****************************************************************** */
void Transport_PK::check_tracer_bounds(Epetra_MultiVector& tracer, 
                                       int component,
                                       double lower_bound,
                                       double upper_bound,
                                       double tol) const
{ 
  Teuchos::RCP<Epetra_MultiVector> tcc = TS->get_total_component_concentration();

  for (int c=cmin; c<cmax_owned; c++) {
    double value = tracer[component][c];
    if (value < lower_bound - tol || value > upper_bound + tol) {
      cout << "Transport_PK: tracer violates bounds" << endl; 
      cout << "    Make an Amanzi ticket or turn off internal transport tests" << endl;
      cout << "    MyPID = " << MyPID << endl;
      cout << "    component = " << component << endl;
      cout << "    internal time = " << T_internal << endl;
      cout << "      cell = " << c << endl;
      cout << "      center = " << mesh_->cell_centroid(c) << endl;
      cout << "      limiter = " << (*limiter_)[c] << endl;
      cout << "      value (old) = " << (*tcc)[component][c] << endl;
      cout << "      value (new) = " << value << endl;
 
      Errors::Message msg;
      msg << "Tracer violates bounds." << "\n";
      Exceptions::amanzi_throw(msg);
    }
  }
}

}  // namespace AmanziTransport
}  // namespace Amanzi

