/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Ethan Coon

   Base class for advection.
   ------------------------------------------------------------------------- */

namespace Amanzi {
namespace Operators {

void Advection::set_flux(Teuchos::RCP<const CompositeVector>& flux) {
  // check that flux includes FACES and has one dof
  flux_ = flux;
}

void Adevection::set_num_dofs(int num_dofs) {
  if (field_ == Teuchos::null || num_dofs_ != num_dofs) {
    num_dofs_ = num_dofs;

    std::vector<std::string> names(2);
    names[0] = "cell";
    names[1] = "face";

    std::vector<AmanziMesh::Entity_kind> locations(2);
    locations[0] = AmanziMesh::CELL;
    locations[1] = AmanziMesh::FACE;

    field_ = Teuchos::rcp(new CompositeVector(names, locations, num_dofs_, true));
    field_->CreateData();
  }
}

} // namespace Operators
} // namespace Amanzi
