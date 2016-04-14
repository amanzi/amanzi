/*
  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
  Daniil Svyatskiy (dasvyat@lanl.gov)
         Ethan Coon (ATS version) (ecoon@lanl.gov)

  The class provides a different implementation of solvers than in
  the base class. In particular, Lagrange multipliers are elliminated
  from the DAE system and short vectors are used in the nonlinear solver.

  MatrixMFD_TPFA_ScaledConstraint provides an TPFA implementation where the
  constraint equation can be rescaled (by, for instance, the faces rel perm)
  which allows zero diffusion coefficient.
*/

#ifndef OPERATORS_MATRIX_MFD_TPFA_SCALED_CONSTRAINT_HH_
#define OPERATORS_MATRIX_MFD_TPFA_SCALED_CONSTRAINT_HH_

#include <strings.h>

#include "MatrixMFD.hh"
#include "MatrixMFD_ScaledConstraint.hh"
#include "MatrixMFD_TPFA.hh"

namespace Amanzi {
namespace Operators {

class MatrixMFD_TPFA_ScaledConstraint : public MatrixMFD_TPFA,
                                        public MatrixMFD_ScaledConstraint {
 public:
  MatrixMFD_TPFA_ScaledConstraint(Teuchos::ParameterList& plist,
          const Teuchos::RCP<const AmanziMesh::Mesh> mesh) :
      MatrixMFD_TPFA(plist,mesh),
      MatrixMFD_ScaledConstraint(plist,mesh),
      MatrixMFD(plist,mesh) {}

  virtual void CreateMFDstiffnessMatrices(
      const Teuchos::Ptr<const CompositeVector>& Krel);

};

} // namespace
} // namespace

#endif
