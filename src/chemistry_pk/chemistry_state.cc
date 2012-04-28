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
    : simulation_state_(S),
      aqueous_components_(Teuchos::null),
      free_ion_species_(Teuchos::null),
      mineral_volume_fractions_(Teuchos::null),
      mineral_specific_surface_area_(Teuchos::null),
      ion_exchange_sites_(Teuchos::null),
      sorption_sites_(Teuchos::null),
      total_sorbed_(Teuchos::null) {

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


}  // namespace chemistry
}  // namespace amanzi
