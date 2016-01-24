/*
  Chemistry PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.
*/
 
#ifndef CHEMISTRY_AMANZI_PK_HH_
#define CHEMISTRY_AMANZI_PK_HH_

#include <string>
#include <vector>

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

// Amanzi
#include "beaker.hh"
#include "chemistry_exception.hh"
#include "chemistry_verbosity.hh"
#include "Chemistry_PK.hh"
#include "Mesh.hh"

// forward declarations
class Epetra_MultiVector;
class Epetra_Vector;
class Epetra_SerialDenseVector;

namespace Amanzi {
namespace AmanziChemistry {

// Trilinos based chemistry process kernel for the unstructured mesh
class Amanzi_PK : public Chemistry_PK {
 public:
  Amanzi_PK(const Teuchos::RCP<Teuchos::ParameterList>& glist,
            Teuchos::RCP<State> S,
            Teuchos::RCP<const AmanziMesh::Mesh> mesh);

  ~Amanzi_PK();

  // members required by PK interface
  virtual void Setup();
  virtual void Initialize();

  void Advance(const double& delta_time,
               Teuchos::RCP<Epetra_MultiVector> total_component_concentration);
  void CommitState(const double& time);

  // modifiers
  void set_max_time_step(const double mts) { this->max_time_step_ = mts; }
  double time_step(void) const { return this->max_time_step_; }

  // The following two routines provide the interface for
  // output of auxillary cellwise data from chemistry
  Teuchos::RCP<Epetra_MultiVector> extra_chemistry_output_data();
  void set_chemistry_output_names(std::vector<std::string>* names);

 private:
  void AllocateAdditionalChemistryStorage_(const Beaker::BeakerComponents& components);

  void XMLParameters();
  void SetupAuxiliaryOutput();
  void SizeBeakerStructures_();

  void CopyCellStateToBeakerStructures(
      int cell_id, Teuchos::RCP<Epetra_MultiVector> total_component_concentration);
  void CopyBeakerStructuresToCellState(
      int cell_id, Teuchos::RCP<Epetra_MultiVector> total_component_concentration);

 private:
  Teuchos::RCP<Teuchos::ParameterList> cp_list_;

  Beaker* chem_;
  Beaker::BeakerParameters beaker_parameters_;
  Beaker::BeakerComponents beaker_components_;
  Beaker::BeakerComponents beaker_components_copy_;

  double max_time_step_;
  double current_time_;
  double saved_time_;

  std::vector<std::string> aux_names_;
  std::vector<int> aux_index_;
  Teuchos::RCP<Epetra_MultiVector> aux_data_;
};

}  // namespace AmanziChemistry
}  // namespace Amanzi
#endif
