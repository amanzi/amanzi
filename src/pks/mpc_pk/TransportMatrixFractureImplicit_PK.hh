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

#include "Epetra_BlockMap.h"
#include "Teuchos_RCP.hpp"

#include "PK_MPCStrong.hh"
#include "PDE_AdvectionUpwind.hh"
#include "PDE_CouplingFlux.hh"
#include "TreeOperator.hh"

#include "FractureInsertion.hh"
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
  virtual void Setup() override;
  virtual void Initialize() override;

  // -- advance each sub pk from t_old to t_new.
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit = false) override;
  virtual void CommitStep(double t_old, double t_new, const Tag& tag) override;

  // miscaleneous methods
  virtual std::string name() override { return "coupled transport implicit"; } 

 private:
  bool AdvanceStepLO_(double t_old, double t_new, int* tot_itrs);
  bool AdvanceStepHO_(double t_old, double t_new, int* tot_itrs);

 private:
  const Teuchos::RCP<Teuchos::ParameterList> glist_;
  Teuchos::RCP<Teuchos::ParameterList> tp_list_;

  Teuchos::RCP<const AmanziMesh::Mesh> mesh_domain_, mesh_fracture_;

  Key matrix_vol_flowrate_key_, fracture_vol_flowrate_key_;

  int num_aqueous_;
  Teuchos::RCP<BDF1_TI<TreeVector, TreeVectorSpace> > bdf1_dae_;

  Teuchos::RCP<Operators::PDE_CouplingFlux> op_coupling00_, op_coupling01_;
  Teuchos::RCP<Operators::PDE_CouplingFlux> op_coupling10_, op_coupling11_;

  bool flag_dispersion_;
  Teuchos::RCP<FractureInsertion> fid_;
  Teuchos::RCP<Operators::PDE_CouplingFlux> op_coupling00d_, op_coupling01d_;
  Teuchos::RCP<Operators::PDE_CouplingFlux> op_coupling10d_, op_coupling11d_;

  Teuchos::RCP<TimestepController> ts_control_;

  // factory registration
  static RegisteredPKFactory<TransportMatrixFractureImplicit_PK> reg_;
};

}  // namespace Amanzi
#endif
