/*
  MPC PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov

  Insertion of fracture into a model changes mesh topology, mesh maps, etc.
  This class provides supporting tools and data structures that can be
  shared between MPC PKs.
*/

#ifndef AMANZI_FRACTURE_INSERTION_HH_
#define AMANZI_FRACTURE_INSERTION_HH_

#include <memory>

#include "Epetra_BlockMap.h"
#include "Epetra_MultiVector.h"
#include "Teuchos_RCP.hpp"

#include "CompositeVectorSpace.hh"
#include "Mesh.hh"

namespace Amanzi {

class FractureInsertion {
 public:
  FractureInsertion(Teuchos::RCP<const AmanziMesh::Mesh>& mesh_matrix,
                    Teuchos::RCP<const AmanziMesh::Mesh>& mesh_fracture);

  // inialization tools
  // -- matrix faces coupled to fracture cells
  void InitMatrixFaceToFractureCell(Teuchos::RCP<const Epetra_BlockMap> mmap,
                                    Teuchos::RCP<const Epetra_BlockMap> gmap);
  // -- matrix cells coupled to fracture cells
  void InitMatrixCellToFractureCell();

  // -- compute coupling coefficients
  void SetValues(const Epetra_MultiVector& values, double scale);
  void SetValues(const CompositeVector& flux);

  // access
  const Teuchos::RCP<CompositeVectorSpace>& get_cvs_matrix() { return cvs_matrix_; }
  const Teuchos::RCP<CompositeVectorSpace>& get_cvs_fracture() { return cvs_fracture_; }

  const std::shared_ptr<std::vector<std::vector<int>>>& get_inds_matrix() { return inds_matrix_; }
  const std::shared_ptr<std::vector<std::vector<int>>>& get_inds_fracture()
  {
    return inds_fracture_;
  }

  const std::shared_ptr<std::vector<double>>& get_values() { return values_; }
  const std::shared_ptr<std::vector<double>>& get_values2() { return values2_; }

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_matrix_, mesh_fracture_;

  Teuchos::RCP<CompositeVectorSpace> cvs_matrix_, cvs_fracture_;
  Teuchos::RCP<const Epetra_BlockMap> mmap_;

  std::shared_ptr<std::vector<std::vector<int>>> inds_matrix_, inds_fracture_;
  std::shared_ptr<std::vector<double>> values_, values2_;
};

} // namespace Amanzi
#endif
