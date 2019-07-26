/*
  MPC PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
           Daniil Svyatskiy

  Process kernel that couples transport in matrix and fracture
  using implicit scheme.
*/

#ifndef AMANZI_TRANSPORT_MATRIX_FRACTURE_IMPLICIT_PK_HH_
#define AMANZI_TRANSPORT_MATRIX_FRACTURE_IMPLICIT_PK_HH_

#include "Teuchos_RCP.hpp"

#include "PK_MPCStrong.hh"
#include "PDE_AdvectionUpwind.hh"
#include "PDE_CouplingFlux.hh"
#include "TreeOperator.hh"

#include "TransportMatrixFracture_PK.hh"

namespace Amanzi {

class TransportMatrixFractureImplicit_PK : public PK_MPCStrong<PK_BDF> {
 public:
  TransportMatrixFractureImplicit_PK(Teuchos::ParameterList& pk_tree,
                                     const Teuchos::RCP<Teuchos::ParameterList>& glist,
                                     const Teuchos::RCP<State>& S,
                                     const Teuchos::RCP<TreeVector>& soln);

  // PK methods
  // -- setup
  virtual void Setup(const Teuchos::Ptr<State>& S);
  virtual void Initialize(const Teuchos::Ptr<State>& S);

  // -- advance each sub pk from t_old to t_new.
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit = false);

  // miscaleneous methods
  virtual std::string name() { return "coupled transport implicit"; } 

 private:
  const Teuchos::RCP<Teuchos::ParameterList> glist_;
  Teuchos::RCP<VerboseObject> vo_;

  Teuchos::RCP<const AmanziMesh::Mesh> mesh_domain_, mesh_fracture_;

  Teuchos::RCP<Operators::PDE_CouplingFlux> op_coupling00_, op_coupling01_;
  Teuchos::RCP<Operators::PDE_CouplingFlux> op_coupling10_, op_coupling11_;
  Teuchos::RCP<Operators::TreeOperator> op_tree_;

  // factory registration
  static RegisteredPKFactory<TransportMatrixFractureImplicit_PK> reg_;
};

}  // namespace Amanzi
#endif
