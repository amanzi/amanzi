/*
  MPC PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
           Daniil Svyatskiy

  Process kernel that couples transport in matrix and fracture.
*/

#ifndef AMANZI_TRANSPORT_MATRIX_FRACTURE_PK_HH_
#define AMANZI_TRANSPORT_MATRIX_FRACTURE_PK_HH_

#include "Teuchos_RCP.hpp"

#include "EvaluatorIndependentFunction.hh"
#include "EvaluatorSecondary.hh"
#include "PK_BDF.hh"
#include "PK_MPCWeak.hh"
#include "PK_Factory.hh"

namespace Amanzi {

class TransportMatrixFracture_PK : public PK_MPCWeak {
 public:
  TransportMatrixFracture_PK(Teuchos::ParameterList& pk_tree,
                             const Teuchos::RCP<Teuchos::ParameterList>& glist,
                             const Teuchos::RCP<State>& S,
                             const Teuchos::RCP<TreeVector>& soln);

  // PK methods
  // -- setup
  virtual void Setup(const Teuchos::Ptr<State>& S);

  // -- dt is the minimum of the sub pks
  virtual double get_dt();
  // virtual void set_dt(double dt);

  // -- advance each sub pk from t_old to t_new.
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit = false);
  // virtual void CommitStep(double t_old, double t_new);

  // -- miscaleneous methods
  virtual std::string name() { return "coupled transport"; } 
  // virtual void CalculateDiagnostics() {};

 private:
  const Teuchos::RCP<Teuchos::ParameterList>& glist_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_domain_, mesh_fracture_;

  Teuchos::RCP<EvaluatorIndependentFunction> matrix_bc;
  Teuchos::RCP<EvaluatorIndependentFunction> fracture_src;

  // factory registration
  static RegisteredPKFactory<TransportMatrixFracture_PK> reg_;
};

}  // namespace Amanzi
#endif
