/*
This is the flow component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Authors: Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

#ifndef AMANZI_FLOW_MATRIX_BASE_HH_
#define AMANZI_FLOW_MATRIX_BASE_HH_

#include <strings.h>

#include "Epetra_Map.h"
#include "Epetra_Operator.h"
#include "Epetra_MultiVector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_FECrsMatrix.h"

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "Point.hh"
#include "boundary_function.hh"

#include "Flow_typedefs.hh"
#include "RelativePermeability.hh"


namespace Amanzi {
namespace AmanziFlow {

class Matrix : public Epetra_Operator {
 public:
  Matrix(Teuchos::RCP<const AmanziMesh::Mesh> mesh, const Epetra_Map& map) : 
      mesh_(mesh), map_(map) {};
  ~Matrix();

  // main methods
  virtual void CreateMassMatrices(std::vector<WhetStone::Tensor>& K) = 0;
  virtual void CreateStiffnessMatrices(std::vector<WhetStone::Tensor>& K) = 0;

  virtual void SymbolicAssemble();
  virtual void Assemble();

  virtual void ApplyBoundaryConditions(std::vector<int>& bc_model, std::vector<bc_tuple>& bc_values) = 0;

  virtual int Apply(const Epetra_MultiVector& v, Epetra_MultiVector& av) const = 0;
  virtual int ApplyInverse(const Epetra_MultiVector& v, Epetra_MultiVector& hv) const = 0;

  virtual void InitPreconditioner(std::string prec_name, Teuchos::ParameterList& prec_list);
  virtual void UpdatePreconditioner();
  virtual void DestroyPreconditioner();

  virtual void DeriveDarcyMassFlux(const Epetra_Vector& v, Epetra_Vector& flux) = 0;

  virtual Teuchos::RCP<Epetra_FECrsMatrix>& matrix() const = 0;
  virtual Teuchos::RCP<Epetra_Vector>& rhs() const = 0;

  // universal members
  int property() { return property_; }

  bool UseTranspose() const { return false; }
  int SetUseTranspose(bool) { return 1; }

  const Epetra_Comm& Comm() const { return *(mesh_->get_comm()); }
  const Epetra_Map& OperatorDomainMap() const { return map_; }
  const Epetra_Map& OperatorRangeMap() const { return map_; }

  const char* Label() const { return strdup("Flow Matrix"); }
  double NormInf() const { return 0.0; }
  bool HasNormInf() const { return false; }

 protected:
  double ComputeResidual(const Epetra_Vector& v, Epetra_Vector& r) {};
  double ComputeNegativeResidual(const Epetra_Vector& v, Epetra_Vector& r) {};

  void AddProperty(int property) { property_ += property; }

 protected:
  Epetra_Map map_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;

  int property_;
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
