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


namespace amanzi {
namespace chemistry {

Chemistry_State::Chemistry_State(Teuchos::RCP<State> S)
    : total_component_concentration_(S->get_total_component_concentration()),
      porosity_(S->get_porosity()),
      water_density_(S->get_water_density()),
      water_saturation_(S->get_water_saturation()),
      mesh_maps_(S->get_mesh_maps()),
      volume_(Teuchos::null),
      aqueous_components_(Teuchos::null),
      free_ion_species_(Teuchos::null),
      mineral_volume_fractions_(Teuchos::null),
      mineral_specific_surface_area_(Teuchos::null),
      ion_exchange_sites_(Teuchos::null),
      sorption_sites_(Teuchos::null),
      total_sorbed_(Teuchos::null) {

  ExtractVolumeFromMesh();
}  // end Chemistry_State


Chemistry_State::~Chemistry_State() {
}  // end ~Chemistry_State


void Chemistry_State::AllocateMemory(const int num_aqueous,
                                     const int num_free_ion,
                                     const int num_minerals,
                                     const int num_ion_exchange_sites,
                                     const int num_total_sorbed,
                                     const int num_sorption_sites) {
  Epetra_Map mesh_info = mesh_maps()->cell_map(false);

  aqueous_components_ = Teuchos::rcp(new Epetra_MultiVector(mesh_info, num_aqueous));

  free_ion_species_ = Teuchos::rcp(new Epetra_MultiVector(mesh_info, num_free_ion));
  if (num_minerals > 0) {
    mineral_volume_fractions_ = Teuchos::rcp(new Epetra_MultiVector(mesh_info, num_minerals));
    mineral_specific_surface_area_ = Teuchos::rcp(new Epetra_MultiVector(mesh_info, num_minerals));
  }
  if (num_ion_exchange_sites) {
    ion_exchange_sites_ = Teuchos::rcp(new Epetra_MultiVector(mesh_info, num_ion_exchange_sites));
  }
  if (num_sorption_sites) {
    sorption_sites_ = Teuchos::rcp(new Epetra_MultiVector(mesh_info, num_sorption_sites));
  }
  total_sorbed_ = Teuchos::rcp(new Epetra_MultiVector(mesh_info, num_total_sorbed));

}  // end AllocateMemory()

void Chemistry_State::ExtractVolumeFromMesh(void) {
  // one of the mesh calls below requires removing the const from the
  // mesh pointer....
  Teuchos::RCP<const Mesh> const_mesh = mesh_maps();
  Teuchos::RCP<Mesh> mesh =
      Teuchos::rcp_const_cast<Mesh>(const_mesh);

  Epetra_Map mesh_info = mesh->cell_map(false);

  volume_ = Teuchos::rcp(new Epetra_Vector(mesh_info));


  int ncell = mesh->count_entities(Amanzi::AmanziMesh::CELL,
                                   Amanzi::AmanziMesh::OWNED);

  if (ncell != volume_->MyLength()) {
    Exceptions::amanzi_throw(
        ChemistryException("Chemistry_State::ExtractVolumeFromMesh() size error."));
  }

  double xdata[24];  // 8 x 3
  Epetra_SerialDenseMatrix xmatrix(View, xdata, 3, 3, 8);
  for (int j = 0; j < ncell; ++j) {
    mesh->cell_to_coordinates((unsigned int) j, xdata, xdata + 24);
    (*volume_)[j] = cell_geometry::hex_volume(xmatrix);
  }
}  // end ExtractVolumeFromMesh()

}  // namespace chemistry
}  // namespace amanzi
