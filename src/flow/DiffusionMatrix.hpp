#ifndef __DIFFUSIONMATRIX_H__
#define __DIFFUSIONMATRIX_H__

#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_SerialSymDenseMatrix.h"

#include "Mesh_maps_base.hh"

class DiffusionMatrix {
public:

  DiffusionMatrix(const Teuchos::RCP<Mesh_maps_base> &mesh, const std::vector<int> &dir_faces);

  ~DiffusionMatrix();

  void Compute(const std::vector<double> &K);

  void Compute(const std::vector<Epetra_SerialSymDenseMatrix> &K);

  void ComputeFaceSchur();

  void ApplyDirichletProjection(Epetra_FECrsMatrix&) const;

  void ApplyDirichletProjection(Epetra_MultiVector&) const;

  void Print(ostream &os=std::cout) const;

  void add_to_celldiag(const Epetra_Vector &celldiag);


  // Accessors

  const Mesh_maps_base& Mesh() const { return *mesh_; }

  const Epetra_Comm& Comm() const { return *(mesh_->get_comm()); }

  const Epetra_Map& CellMap(bool ghost=false) const { return mesh_->cell_map(ghost); }

  const Epetra_Map& FaceMap(bool ghost=false) const { return mesh_->face_map(ghost); }

  const Epetra_CrsGraph& FaceGraph() const { return Dff_->Graph(); }

  const Epetra_Vector& Dcc() const { return *Dcc_; }

  const Epetra_CrsMatrix& Dcf() const { return *Dcf_; }

  const Epetra_FECrsMatrix& Dff() const { return *Dff_; }

  const Epetra_FECrsMatrix& Sff();

  const std::vector<int>& DirFaces() const { return dir_faces_; }

private:

  Teuchos::RCP<Mesh_maps_base> mesh_;
  Epetra_Vector *Dcc_;
  Epetra_CrsMatrix *Dcf_;
  Epetra_FECrsMatrix *Dff_;
  std::vector<int> dir_faces_;
  Epetra_FECrsMatrix *Sff_;

  template <typename T> void Compute(const std::vector<T> &K);
};

#endif
