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

void
MergeSubsurfaceAndSurfacePressure(const Teuchos::Ptr<CompositeVector>& sub_p,
				  const Teuchos::Ptr<CompositeVector>& sub_T,
				  const Teuchos::Ptr<CompositeVector>& surf_p,
				  const Teuchos::Ptr<CompositeVector>& surf_T) {
  Epetra_MultiVector& sub_p_f = *sub_p->ViewComponent("face",false);
  Epetra_MultiVector& sub_T_f = *sub_T->ViewComponent("face",false);
  Epetra_MultiVector& surf_p_c = *surf_p->ViewComponent("cell",false);
  Epetra_MultiVector& surf_T_c = *surf_T->ViewComponent("cell",false);
  double p_atm = 101325.;

  for (unsigned int sc=0; sc!=surf_p_c.MyLength(); ++sc) {
    AmanziMesh::Entity_ID f =
        surf_p->Mesh()->entity_get_parent(AmanziMesh::CELL, sc);
    if (surf_p_c[0][sc] >= p_atm) {
      sub_p_f[0][f] = surf_p_c[0][sc];
      sub_T_f[0][f] = surf_T_c[0][sc];
    } else {
      surf_p_c[0][sc] = sub_p_f[0][f];
      surf_T_c[0][sc] = sub_T_f[0][f];
    }
  }
}


} // namespace
