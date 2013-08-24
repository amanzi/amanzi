/*
This is the transport component of the Amanzi code. 

Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Author: Konstantin Lipnikov (lipnikov@lanl.gov)
Usage: 
*/

#ifndef __MATRIX_DISPERSION_HH__
#define __MATRIX_DISPERSION_HH__

#include "Epetra_Vector.h"

namespace Amanzi {
namespace AmanziTransport {

class Matrix_Dispersion {
 public:
  Matrix_Dispersion();
  ~Matrix_Dispersion();

  // primary members
  void Apply();
  void ApplyInverse();

  // I/O methods
  // void ProcessStringDispersionModel(const std::string name, int* method);

 private:
  void CalculateDispersionTensor();
  void ExtractBoundaryConditions(const int component,
                                 std::vector<int>& bc_face_id,
                                 std::vector<double>& bc_face_value);
  void PopulateHarmonicPointsValues(int component,
                                    Teuchos::RCP<Epetra_MultiVector> tcc,
                                    std::vector<int>& bc_face_id,
                                    std::vector<double>& bc_face_values);
  void AddDispersiveFluxes(int component,
                           Teuchos::RCP<Epetra_MultiVector> tcc,
                           std::vector<int>& bc_face_id,
                           std::vector<double>& bc_face_values,
                           Teuchos::RCP<Epetra_MultiVector> tcc_next);

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
};

}  // namespace AmanziTransport
}  // namespace Amanzi

#endif

