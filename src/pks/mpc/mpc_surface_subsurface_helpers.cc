#include "mpc_surface_subsurface_helpers.hh"

namespace Amanzi {

void
CopySurfaceToSubsurface(const CompositeVector& surf,
                        const Teuchos::Ptr<CompositeVector>& sub) {
  const Epetra_MultiVector& surf_c = *surf.ViewComponent("cell",false);
  Epetra_MultiVector& sub_f = *sub->ViewComponent("face",false);

  for (unsigned int sc=0; sc!=surf_c.MyLength(); ++sc) {
    AmanziMesh::Entity_ID f =
        surf.Mesh()->entity_get_parent(AmanziMesh::CELL, sc);
    sub_f[0][f] = surf_c[0][sc];
  }
}

void
CopySubsurfaceToSurface(const CompositeVector& sub,
                        const Teuchos::Ptr<CompositeVector>& surf) {
  const Epetra_MultiVector& sub_f = *sub.ViewComponent("face",false);
  Epetra_MultiVector& surf_c = *surf->ViewComponent("cell",false);

  for (unsigned int sc=0; sc!=surf_c.MyLength(); ++sc) {
    AmanziMesh::Entity_ID f =
        surf->Mesh()->entity_get_parent(AmanziMesh::CELL, sc);
    surf_c[0][sc] = sub_f[0][f];
  }
}

} // namespace
