#ifndef __MIMETICHEX_H__
#define __MIMETICHEX_H__

#include "Teuchos_RCP.hpp"

#include "Epetra_Vector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"

#include "Mesh.hh"

namespace Amanzi {

class MimeticHex {
 public:
  MimeticHex(const Teuchos::RCP<Amanzi::AmanziMesh::Mesh> &mesh);
  ~MimeticHex(){}

  const Epetra_Comm& Comm() const { return *(mesh_->get_comm()); }

  const Epetra_Map& CellMap(bool ghost=false) const { return mesh_->cell_map(ghost); }
  const Epetra_Map& FaceMap(bool ghost=false) const { return mesh_->face_map(ghost); }

  void DeriveFluxes(const Epetra_Vector&, Epetra_Vector&) const;

  Epetra_SerialDenseVector face_area_;

 private:
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh_;

  Epetra_SerialDenseVector volume_;
  Epetra_SerialDenseMatrix face_normal_;
};

}  // close namespace Amanzi

#endif
