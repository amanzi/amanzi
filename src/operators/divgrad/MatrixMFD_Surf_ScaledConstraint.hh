/*
  License: BSD
  Ethan Coon (ecoon@lanl.gov)

  MatrixMFD_Surf_ScaledConstraint provides for solving the p-lambda system where a subset of
  the lambdas are more densely coupled for overland flow.

*/

#ifndef OPERATORS_MATRIX_MFD_SURF_SCALED_CONSTRAINT_HH_
#define OPERATORS_MATRIX_MFD_SURF_SCALED_CONSTRAINT_HH_

#include "MatrixMFD_TPFA.hh"
#include "MatrixMFD_Surf.hh"
#include "MatrixMFD_ScaledConstraint.hh"

namespace Amanzi {
namespace Operators {

class MatrixMFD_Surf_ScaledConstraint : public MatrixMFD_Surf,
                                        public MatrixMFD_ScaledConstraint {

 public:
  MatrixMFD_Surf_ScaledConstraint(Teuchos::ParameterList& plist,
				  const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);

  virtual void AssembleGlobalMatrices();

  virtual void ApplyBoundaryConditions(const std::vector<MatrixBC>& subsurface_markers,
          const std::vector<double>& subsurface_values);

  virtual void ComputeSchurComplement(const std::vector<MatrixBC>& bc_markers,
          const std::vector<double>& bc_values);

};


} //namespace
} //namespace


#endif
