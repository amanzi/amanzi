/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
Amanzi Flow

License: see COPYRIGHT
Author: Ethan Coon

Interface layer between Flow and State, this is a harness for
accessing the new state-dev from the old Flow PK.

 ------------------------------------------------------------------------- */

#ifndef AMANZI_PK_STATE_HH_
#define AMANZI_PK_STATE_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Mesh.hh"
#include "State.hh"

namespace Amanzi {

class PK_State {

public:
  enum PKStateConstructMode {
    PK_STATE_CONSTRUCT_MODE_COPY_POINTERS,
    PK_STATE_CONSTRUCT_MODE_VIEW_DATA,
    PK_STATE_CONSTRUCT_MODE_VIEW_DATA_GHOSTED,
    PK_STATE_CONSTRUCT_MODE_COPY_DATA,
    PK_STATE_CONSTRUCT_MODE_COPY_DATA_GHOSTED,
  };

  PK_State(std::string name, Teuchos::RCP<AmanziMesh::Mesh> mesh);
  PK_State(std::string name, Teuchos::RCP<State> S);
  PK_State(std::string name, State& S);
  PK_State(PK_State& other, StateConstructMode mode=STATE_CONSTRUCT_MODE_COPY_POINTERS);

  Teuchos::RCP<const AmanziMesh::Mesh> mesh() { return mesh_; }

  // data management
  void CopyMasterCell2GhostCell(Epetra_Vector& v);
  void CopyMasterCell2GhostCell(const Epetra_Vector& v, Epetra_Vector& vhost);
  void CopyMasterFace2GhostFace(Epetra_Vector& v);
  void CopyMasterFace2GhostFace(const Epetra_Vector& v, Epetra_Vector& vhost);
  void CopyMasterMultiCell2GhostMultiCell(Epetra_MultiVector& v);
  void CopyMasterMultiCell2GhostMultiCell(const Epetra_MultiVector& v,
          Epetra_MultiVector& vv, int parallel_comm = 1);
  void CombineGhostFace2MasterFace(Epetra_Vector& v, Epetra_CombineMode mode = Insert);
  void CombineGhostCell2MasterCell(Epetra_Vector& v, Epetra_CombineMode mode = Insert);

  Epetra_Vector* CreateCellView(const Epetra_Vector& u) const;
  Epetra_Vector* CreateFaceView(const Epetra_Vector& u) const;

  double normLpCell(const Epetra_Vector& v1, double p);
  double normLpCell(const Epetra_Vector& v1, const Epetra_Vector& v2, double p);

  // extension of Trilinos
  void MinValueMasterCells(Epetra_MultiVector& v, double* vmin);
  void MaxValueMasterCells(Epetra_MultiVector& v, double* vmax);

protected:

  Teuchos::RCP<State> S_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  std::string name_;
  bool ghosted_;
};


} // namespace

#endif
