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
  DiffusionMatrix(Teuchos::RCP<Mesh_maps_base> &mesh, std::vector<int> &dir_faces);
  ~DiffusionMatrix();

  void Compute(const std::vector<double> &K);
  void Compute(const std::vector<Epetra_SerialSymDenseMatrix> &K);
  void Print();

  void ComputeFaceSchur();

  Teuchos::RCP<Mesh_maps_base> &GetMesh() { return mesh_;}
  const Mesh_maps_base& Mesh() { return *mesh_;}

  const Epetra_Comm& Comm() { return *(mesh_->get_comm()); }

  const Epetra_CrsGraph& FaceGraph() { return Dff_->Graph(); }

  void ApplyDirichletProjection(Epetra_CrsMatrix&);
  void ApplyDirichletProjection(Epetra_MultiVector&);

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
};

#endif
