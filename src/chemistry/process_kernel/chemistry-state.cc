/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#include "chemistry_state.hh"

#include "Epetra_SerialDenseVector.h"
#include "Epetra_MultiVector.h"
#include "Teuchos_RCPDecl.hpp"

#include "State.hpp"

#include "cell_geometry.hh"

#include "chemistry_exception.hh"
#include "errors.hh"
#include "exceptions.hh"


Chemistry_State::Chemistry_State(Teuchos::RCP<State> S)
  : total_component_concentration_(S->get_total_component_concentration()),
    porosity_(S->get_porosity()),
    water_density_(S->get_water_density()),
    water_saturation_(S->get_water_saturation()),
    mesh_maps_(S->get_mesh_maps()) {
  // TODO(bandre): can we make this the same type as the other state vectors
  volume_ =
      Teuchos::rcp(new Epetra_SerialDenseVector(
          get_mesh_maps()->count_entities(Mesh_data::CELL, OWNED)));

  ExtractVolumeFromMesh();
}  // end Chemistry_State


Chemistry_State::~Chemistry_State() {
}  // end ~Chemistry_State


void Chemistry_State::ExtractVolumeFromMesh(void) {
  Teuchos::RCP<const Mesh_maps_base> const_mesh = get_mesh_maps();

  // one of the mesh calls below requires removing the const from the
  // mesh pointer....
  Teuchos::RCP<Mesh_maps_base> mesh =
      Teuchos::rcp_const_cast<Mesh_maps_base>(const_mesh);

  int ncell = mesh->count_entities(Mesh_data::CELL, OWNED);

  if (ncell != volume_->Length()) {
    Exceptions::amanzi_throw(
        ChemistryException("Chemistry_State::ExtractVolumeFromMesh() size error."));
  }

  double xdata[24];  // 8 x 3
  Epetra_SerialDenseMatrix xmatrix(View, xdata, 3, 3, 8);
  for (int j = 0; j < ncell; ++j) {
    mesh->cell_to_coordinates((unsigned int) j, xdata, xdata+24);
    (*volume_)[j] = cell_geometry::hex_volume(xmatrix);
  }
}  // end ExtractVolumeFromMesh()
