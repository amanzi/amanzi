/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_MECHANICS_FRACTURED_MATRIX_PK_HH_
#define AMANZI_MECHANICS_FRACTURED_MATRIX_PK_HH_

#include <string>
#include <vector>

// TPLs
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

// Amanzi
#include "Mesh.hh"

// Amanzi::Mechanics
#include "MechanicsSmallStrain_PK.hh"

namespace Amanzi {
namespace Mechanics {

class MechanicsFracturedMatrix_PK : public MechanicsSmallStrain_PK {
 public:
  MechanicsFracturedMatrix_PK(Teuchos::ParameterList& pk_tree,
                              const Teuchos::RCP<Teuchos::ParameterList>& glist,
                              const Teuchos::RCP<State>& S,
                              const Teuchos::RCP<TreeVector>& soln);

  ~MechanicsFracturedMatrix_PK(){};

  // methods required for PK interface
  virtual void Setup() final;
  virtual void Initialize() final;
  virtual void CommitStep(double t_old, double t_new, const Tag& tag) final;

  // -- computes the non-linear functional f = f(t,u,udot) and related norm.
  void FunctionalResidual(const double t_old,
                          double t_new,
                          Teuchos::RCP<const TreeVector> u_old,
                          Teuchos::RCP<TreeVector> u_new,
                          Teuchos::RCP<TreeVector> f) final;

 private:
  void AddFractureMatrices_(CompositeVector& rhs);

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_fracture_;

  Key aperture_key_, ref_aperture_key_, compliance_key_, pressure_key_;

 private:
  static RegisteredPKFactory<MechanicsFracturedMatrix_PK> reg_;
};

} // namespace Mechanics
} // namespace Amanzi

#endif
