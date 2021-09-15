/*
  This is the mpc_pk component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy

  Process kernel for coupling Transport and Chemistry PKs in the 
  matrix and fracture network.
*/

#ifndef AMANZI_REACTIVE_TRANSPORT_MATRIX_FRACTURE_PK_HH_
#define AMANZI_REACTIVE_TRANSPORT_MATRIX_FRACTURE_PK_HH_

// TPLs
#include "Teuchos_RCP.hpp"

// Amanzi
#include "PDE_DiffusionFracturedMatrix.hh"
#include "PK.hh"
#include "PK_Factory.hh"
#include "PK_MPCAdditive.hh"
#include "PK_MPCWeak.hh"

// Amanzi::MPC
#include "ChemistryMatrixFracture_PK.hh"
#include "TransportMatrixFracture_PK.hh"

namespace Amanzi {

class ReactiveTransportMatrixFracture_PK : public PK_MPCAdditive<PK> {
 public:
  ReactiveTransportMatrixFracture_PK(
      Teuchos::ParameterList& pk_tree,
      const Teuchos::RCP<Teuchos::ParameterList>& global_list,
      const Teuchos::RCP<State>& S,
      const Teuchos::RCP<TreeVector>& soln);

  ~ReactiveTransportMatrixFracture_PK() {};

  // PK methods
  // -- dt is the minimum of the sub pks
  virtual double get_dt() override;
  virtual void set_dt(double dt) override;

  virtual void Setup(const Teuchos::Ptr<State>& S) override;
  virtual void Initialize(const Teuchos::Ptr<State>& S) override;  

  // -- advance each sub pk from t_old to t_new.
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit = false) override;

  std::string name() override { return "reactive transport matrix fracture";}

 protected:  
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_domain_, mesh_fracture_;

  Teuchos::RCP<ChemistryMatrixFracture_PK> coupled_chemistry_pk_;

  Key tcc_matrix_key_, tcc_fracture_key_;

  // factory registration
  static RegisteredPKFactory<ReactiveTransportMatrixFracture_PK> reg_;
};

}  // namespace Amanzi
#endif
