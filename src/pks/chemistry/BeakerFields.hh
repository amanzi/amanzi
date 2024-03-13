/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  Chemistry PK

  A list of beaker fields to avoid search
*/

#ifndef CHEMISTRY_BEAKER_FIELDS_HH_
#define CHEMISTRY_BEAKER_FIELDS_HH_

namespace Amanzi {
namespace AmanziChemistry {

struct BeakerFields {
  Teuchos::RCP<const Epetra_MultiVector> porosity;
  Teuchos::RCP<const Epetra_MultiVector> density;
  Teuchos::RCP<const Epetra_MultiVector> saturation, prev_saturation;
  Teuchos::RCP<const Epetra_MultiVector> temperature;

  Teuchos::RCP<Epetra_MultiVector> free_ion;
  Teuchos::RCP<Epetra_MultiVector> activity, secondary_activity;
  Teuchos::RCP<Epetra_MultiVector> mineral_vf, mineral_ssa;

  Teuchos::RCP<Epetra_MultiVector> sorbed;
  Teuchos::RCP<Epetra_MultiVector> isotherm_kd, isotherm_freundlich_n, isotherm_langmuir_b;

  Teuchos::RCP<Epetra_MultiVector> sorption_sites;
  Teuchos::RCP<Epetra_MultiVector> surface_complex;

  Teuchos::RCP<Epetra_MultiVector> ion_exchange_sites, ion_exchange_ref_cation_conc;
};

} // namespace AmanziChemistry
} // namespace Amanzi

#endif
