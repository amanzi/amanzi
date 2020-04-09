/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#ifndef AMANZI_OPERATOR_MAGNETIC_DIFFUSION_TM_HH_
#define AMANZI_OPERATOR_MAGNETIC_DIFFUSION_TM_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

// Amanzi::Operators
#include "PDE_MagneticDiffusion.hh"

namespace Amanzi {
namespace Operators {

class PDE_MagneticDiffusion_TM : public PDE_MagneticDiffusion {
 public:
  PDE_MagneticDiffusion_TM(const Teuchos::RCP<Operator>& global_op)
    : PDE_MagneticDiffusion(global_op){};

  PDE_MagneticDiffusion_TM(Teuchos::ParameterList& plist,
                           const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
    : PDE_MagneticDiffusion(plist, mesh){};

  // main virtual members
  // -- before solving the problem
  virtual void
  ModifyMatrices(CompositeVector& E, CompositeVector& B, double dt);
  virtual void ApplyBCs(bool primary, bool eliminate, bool essential_eqn);

  // -- after solving the problem
  virtual void ModifyFields(CompositeVector& E, CompositeVector& B, double dt);

  // -- postprocessing
  virtual double CalculateOhmicHeating(const CompositeVector& E);

 private:
  void ApplyBCs_Node_(const Teuchos::Ptr<const BCs>& bc_f,
                      const Teuchos::Ptr<const BCs>& bc_v, bool primary,
                      bool eliminate, bool essential_eqn);
};

} // namespace Operators
} // namespace Amanzi

#endif
