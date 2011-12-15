/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#ifndef __MIMETICDISCRETIZATION_HH__
#define __MIMETICDISCRETIZATION_HH__

#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Mesh.hh"

#include "MimeticHexLocal.hpp"
#include "MimeticHex.hpp"
#include "DiffusionMatrix.hpp"
#include "DiffusionPrecon.hpp"


//namespace Amanzi {

class MimeticDiscretization {

public:
  MimeticDiscretization(const Teuchos::RCP<AmanziMesh::Mesh>& mesh,
                        Teuchos::ParameterList& discretization_plist);
  ~MimeticDiscretization();

  void CalcCellVolumes(Teuchos::RCP<Epetra_Vector> &cell_volumes) const;

  // other methods for assembly?
private:
  std::vector<MimeticHexLocal>  MD_;
  MimeticHex *md_;

  // private helper methods
  const Epetra_Map& CellMap(bool ghost=false) const { return mesh_->cell_map(ghost); }
  const Epetra_Map& FaceMap(bool ghost=false) const { return mesh_->face_map(ghost); }
  void init_mimetic_disc_(Teuchos::RCP<AmanziMesh::Mesh>&,
                          std::vector<MimeticHexLocal>&) const;

};

// } // namespace

#endif
