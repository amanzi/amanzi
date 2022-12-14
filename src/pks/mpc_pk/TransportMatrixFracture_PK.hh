/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
      Daniil Svyatskiy
*/

/*
  MPC PK

  Process kernel that couples transport in matrix and fracture.
*/

#ifndef AMANZI_TRANSPORT_MATRIX_FRACTURE_PK_HH_
#define AMANZI_TRANSPORT_MATRIX_FRACTURE_PK_HH_

#include "Teuchos_RCP.hpp"

#include "EvaluatorIndependentFunction.hh"
#include "EvaluatorSecondary.hh"
#include "Key.hh"
#include "PK_BDF.hh"
#include "PK_MPCWeak.hh"
#include "PK_Factory.hh"
#include "TreeOperator.hh"

namespace Amanzi {

class TransportMatrixFracture_PK : public PK_MPCWeak {
 public:
  TransportMatrixFracture_PK(Teuchos::ParameterList& pk_tree,
                             const Teuchos::RCP<Teuchos::ParameterList>& glist,
                             const Teuchos::RCP<State>& S,
                             const Teuchos::RCP<TreeVector>& soln);

  // PK methods
  virtual void Setup() override;
  virtual void Initialize() override;
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit = false) override;
  virtual double get_dt() override; // the minimum of sub-pks's dt

  // -- miscaleneous methods
  virtual std::string name() override { return "coupled transport"; }

 private:
  const Teuchos::RCP<Teuchos::ParameterList>& glist_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_domain_, mesh_fracture_;

  Teuchos::RCP<EvaluatorIndependentFunction> matrix_bc;
  Teuchos::RCP<EvaluatorIndependentFunction> fracture_src;

  bool flag_dispersion_;
  Teuchos::RCP<Operators::TreeOperator> op_dispersion_;

  Key tcc_matrix_key_, tcc_fracture_key_;
  Key normal_diffusion_key_;

  // factory registration
  static RegisteredPKFactory<TransportMatrixFracture_PK> reg_;
};

} // namespace Amanzi
#endif
