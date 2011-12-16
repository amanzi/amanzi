#include "MimeticHex.hpp"
#include "cell_geometry.hh"

namespace Amanzi {

MimeticHex::MimeticHex(const Teuchos::RCP<AmanziMesh::Mesh> &mesh) : mesh_(mesh)
{
  int ncell = mesh->count_entities(AmanziMesh::CELL, AmanziMesh::USED);
  int nface = mesh->count_entities(AmanziMesh::FACE, AmanziMesh::USED);

  { // Cell volumes
    volume_.Size(ncell);
    double xdata[24]; // 8 x 3
    Epetra_SerialDenseMatrix xmatrix(View, xdata, 3, 3, 8);
    for (int j = 0; j < ncell; ++j) {
      mesh->cell_to_coordinates((unsigned int) j, xdata, xdata+24);
      volume_[j] = cell_geometry::hex_volume(xmatrix);
    }
  }

  { // Area-weighted face normals (oriented areas)
    face_normal_.Shape(3, nface);
    double xdata[12]; // 4 x 3
    Epetra_SerialDenseMatrix xmatrix(View, xdata, 3, 3, 4);
    for (int j = 0; j < nface; ++j) {
      mesh->face_to_coordinates((unsigned int) j, xdata, xdata+12);
      cell_geometry::quad_face_normal(xmatrix, face_normal_[j]);
    }
  }

  // Face areas
  face_area_.Size(nface);
  for (int j = 0; j < nface; ++j)
    face_area_[j] = cell_geometry::vector_length(face_normal_[j], 3);
}

}  // namespace Amanzi
