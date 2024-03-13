/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#ifndef CHEMISTRY_DATA__
#define CHEMISTRY_DATA__

#include "Epetra_BlockMap.h"
#include "Epetra_MultiVector.h"
#include "Teuchos_RCP.hpp"

struct chemistry_data {
  chemistry_data(const Epetra_BlockMap& map,
                 int number_of_components,
                 int number_of_minerals,
                 int number_of_secondaries,
                 int number_of_ion_exchange_sites,
                 int number_of_sorption_sites,
                 bool using_sorption,
                 bool using_sorption_isotherms);

  void store(Teuchos::RCP<Epetra_MultiVector> free_ion_concentrations,
             Teuchos::RCP<Epetra_MultiVector> primary_activity_coeff,
             Teuchos::RCP<Epetra_MultiVector> secondary_activity_coeff,
             Teuchos::RCP<Epetra_MultiVector> mineral_volume_fractions,
             Teuchos::RCP<Epetra_MultiVector> mineral_specific_surface_area,
             Teuchos::RCP<Epetra_MultiVector> total_sorbed,
             Teuchos::RCP<Epetra_MultiVector> sorption_sites,
             Teuchos::RCP<Epetra_MultiVector> surface_complex_free_site_conc,
             Teuchos::RCP<Epetra_MultiVector> ion_exchange_sites,
             Teuchos::RCP<Epetra_MultiVector> ion_exchange_ref_cation_conc,
             Teuchos::RCP<Epetra_MultiVector> isotherm_kd,
             Teuchos::RCP<Epetra_MultiVector> isotherm_freundlich_n,
             Teuchos::RCP<Epetra_MultiVector> isotherm_langmuir_b);

  void retrieve(Teuchos::RCP<Epetra_MultiVector> free_ion_concentrations,
                Teuchos::RCP<Epetra_MultiVector> primary_activity_coeff,
                Teuchos::RCP<Epetra_MultiVector> secondary_activity_coeff,
                Teuchos::RCP<Epetra_MultiVector> mineral_volume_fractions,
                Teuchos::RCP<Epetra_MultiVector> mineral_specific_surface_area,
                Teuchos::RCP<Epetra_MultiVector> total_sorbed,
                Teuchos::RCP<Epetra_MultiVector> sorption_sites,
                Teuchos::RCP<Epetra_MultiVector> surface_complex_free_site_conc,
                Teuchos::RCP<Epetra_MultiVector> ion_exchange_sites,
                Teuchos::RCP<Epetra_MultiVector> ion_exchange_ref_cation_conc,
                Teuchos::RCP<Epetra_MultiVector> isotherm_kd,
                Teuchos::RCP<Epetra_MultiVector> isotherm_freundlich_n,
                Teuchos::RCP<Epetra_MultiVector> isotherm_langmuir_b);

  Teuchos::RCP<Epetra_MultiVector> free_ion_concentrations_;
  Teuchos::RCP<Epetra_MultiVector> primary_activity_coeff_;
  Teuchos::RCP<Epetra_MultiVector> secondary_activity_coeff_;
  Teuchos::RCP<Epetra_MultiVector> mineral_volume_fractions_;      // [cell][mineral]
  Teuchos::RCP<Epetra_MultiVector> mineral_specific_surface_area_; // [cell][mineral]
  Teuchos::RCP<Epetra_MultiVector> total_sorbed_;                  // [cell][species]
  Teuchos::RCP<Epetra_MultiVector>
    sorption_sites_; // [cell][site], eventually [cell][mineral][site]
  Teuchos::RCP<Epetra_MultiVector> surface_complex_free_site_conc_;
  Teuchos::RCP<Epetra_MultiVector> ion_exchange_sites_; // CEC, [cell][mineral]
  Teuchos::RCP<Epetra_MultiVector> ion_exchange_ref_cation_conc_;
  Teuchos::RCP<Epetra_MultiVector> isotherm_kd_;           // [cell][species]
  Teuchos::RCP<Epetra_MultiVector> isotherm_freundlich_n_; // [cell][species]
  Teuchos::RCP<Epetra_MultiVector> isotherm_langmuir_b_;   // [cell][species]
};

#endif
