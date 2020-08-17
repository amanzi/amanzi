/*
  Flow PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Neil Carlson (version 1) 
           Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)

  This is a virtual base class for two flow models: Darcy and Richards.
*/

#ifndef AMANZI_FLOW_PK_HH_
#define AMANZI_FLOW_PK_HH_

#include <vector>

// TPLs
#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_FECrsMatrix.h"
#include "Teuchos_RCP.hpp"

// Amanzi
#include "BCs.hh"
#include "BDFFnBase.hh"
#include "Checkpoint.hh"
#include "CompositeVectorSpace.hh"
#include "independent_variable_field_evaluator_fromfunction.hh"
#include "Key.hh"
#include "Operator.hh"
#include "PK_DomainFunction.hh"
#include "PK_PhysicalBDF.hh"
#include "primary_variable_field_evaluator.hh"
#include "Tensor.hh"
#include "Units.hh"
#include "VerboseObject.hh"

// Flow
#include "FlowBoundaryFunction.hh"
#include "FlowDefs.hh"
#include "FlowTypeDefs.hh"

namespace Amanzi {
namespace Flow {

class Flow_PK : public PK_PhysicalBDF {
 public:
  Flow_PK();
  Flow_PK(Teuchos::ParameterList& pk_tree,
                 const Teuchos::RCP<Teuchos::ParameterList>& glist,
                 const Teuchos::RCP<State>& S,
                 const Teuchos::RCP<TreeVector>& soln);
  virtual ~Flow_PK() {};

  // members required by PK interface
  virtual void Setup(const Teuchos::Ptr<State>& S) override;
  virtual void Initialize(const Teuchos::Ptr<State>& S) override;

  // other members of this PK.
  // -- initialize simple fields common for both flow models.
  void UpdateLocalFields_(const Teuchos::Ptr<State>& S);

  // --- management of boundary and source terms
  void UpdateSourceBoundaryData(double t_old, double t_new, const CompositeVector& u);
  void ComputeOperatorBCs(const CompositeVector& u);
  void SeepageFacePFloTran(const CompositeVector& u, int* nseepage, double* area_seepage);
  void SeepageFaceFACT(const CompositeVector& u, int* nseepage, double* area_seepage);

  void AddSourceTerms(CompositeVector& rhs);
  void ComputeWellIndex(Teuchos::ParameterList& spec);
  bool IsWellIndexRequire(Teuchos::ParameterList& spec);

  // -- absolute permeability and derived quantities.
  void SetAbsolutePermeabilityTensor();

  // -- miscallenous members
  void DeriveFaceValuesFromCellValues(const Epetra_MultiVector& ucells, Epetra_MultiVector& ufaces);

  // -- io members
  void OutputTimeHistory(const Teuchos::ParameterList& plist, std::vector<dt_tuple>& dt_history);

  // -- utilities
  double WaterVolumeChangePerSecond(const std::vector<int>& bc_model,
                                    const Epetra_MultiVector& darcy_flux) const;

  // -- V&V
  void VV_ValidateBCs() const;
  void VV_ReportWaterBalance(const Teuchos::Ptr<State>& S) const;
  void VV_ReportSeepageOutflow(const Teuchos::Ptr<State>& S, double dT) const;
  void VV_PrintHeadExtrema(const CompositeVector& pressure) const;
  void VV_PrintSourceExtrema() const;

  // -- extensions 
  int BoundaryFaceGetCell(int f) const;  // of AmanziMesh
  void VerticalNormals(int c, AmanziGeometry::Point& n1, AmanziGeometry::Point& n2);
  virtual double BoundaryFaceValue(int f, const CompositeVector& u);

  // access
  Teuchos::RCP<Operators::BCs> op_bc() { return op_bc_; }
  double seepage_mass() { return seepage_mass_; } // support of unit tests

 private:
  void InitializeFields_();

 protected:
  void InitializeBCsSources_(Teuchos::ParameterList& list);

 public:
  int ncells_owned, ncells_wghost;
  int nfaces_owned, nfaces_wghost;

  double dt_, dt_next_;

  int MyPID;  // parallel information: will be moved to private
  int missed_bc_faces_, dirichlet_bc_faces_;
  int ti_phase_counter;

 public:
  Teuchos::RCP<const Teuchos::ParameterList> linear_operator_list_;
  Teuchos::RCP<const Teuchos::ParameterList> preconditioner_list_;
  Teuchos::RCP<Teuchos::ParameterList> ti_list_;

 protected:
  int dim;

  std::string passwd_;
  bool peaceman_model_;

  // Stationary physical quantatities
  std::vector<WhetStone::Tensor> K; 
  AmanziGeometry::Point gravity_;
  double g_, rho_, molar_rho_, atm_pressure_;
  double flux_units_;  // scaling for flux units from kg to moles.

  Teuchos::RCP<Epetra_Vector> Kxy;
  std::string coordinate_system_;

  // boundary conditions
  std::vector<Teuchos::RCP<FlowBoundaryFunction> > bcs_; 
  int nseepage_prev;

  Teuchos::RCP<Operators::BCs> op_bc_;

  // source terms and liquid balance
  std::vector<Teuchos::RCP<PK_DomainFunction> > srcs;
  mutable double mass_bc, seepage_mass_, mass_initial;

  // field evaluators (MUST GO AWAY lipnikov@lanl.gov)
  Teuchos::RCP<PrimaryVariableFieldEvaluator> darcy_flux_eval_;
  Teuchos::RCP<PrimaryVariableFieldEvaluator> pressure_eval_, pressure_matrix_eval_;

  // DFN model
  bool flow_on_manifold_;  // true for the DFN model
  bool coupled_to_matrix_, coupled_to_fracture_;

  // names of state fields 
  Key pressure_key_;
  Key darcy_flux_key_, specific_storage_key_, specific_yield_key_;
  Key saturation_liquid_key_, prev_saturation_liquid_key_;
  Key porosity_key_, hydraulic_head_key_, pressure_head_key_;
  Key permeability_key_;
  Key darcy_velocity_key_;

  // io
  Utils::Units units_;
  Teuchos::RCP<Teuchos::ParameterList> fp_list_;
};

}  // namespace Flow
}  // namespace Amanzi

#endif
