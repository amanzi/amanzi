/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include "chemistry_data.hh"

chemistry_data::chemistry_data(const Epetra_BlockMap& map,
                               int number_of_components,
                               int number_of_minerals,
                               int number_of_secondaries,
                               int number_of_ion_exchange_sites,
                               int number_of_sorption_sites,
                               bool using_sorption,
                               bool using_sorption_isotherms)
{
  free_ion_concentrations_ = Teuchos::rcp(new Epetra_MultiVector(map, number_of_components));
  primary_activity_coeff_ = Teuchos::rcp(new Epetra_MultiVector(map, number_of_components));

  if (number_of_secondaries > 0) {
    secondary_activity_coeff_ = Teuchos::rcp(new Epetra_MultiVector(map, number_of_secondaries));
  } else {
    secondary_activity_coeff_ = Teuchos::null;
  }

  if (number_of_minerals > 0) {
    mineral_volume_fractions_ = Teuchos::rcp(new Epetra_MultiVector(map, number_of_minerals));
    mineral_specific_surface_area_ = Teuchos::rcp(new Epetra_MultiVector(map, number_of_minerals));
  } else {
    mineral_volume_fractions_ = Teuchos::null;
    mineral_specific_surface_area_ = Teuchos::null;
  }

  if (using_sorption) {
    total_sorbed_ = Teuchos::rcp(new Epetra_MultiVector(map, number_of_components));
  } else {
    total_sorbed_ = Teuchos::null;
  }

  if (number_of_ion_exchange_sites > 0) {
    ion_exchange_sites_ = Teuchos::rcp(new Epetra_MultiVector(map, number_of_ion_exchange_sites));
    ion_exchange_ref_cation_conc_ =
      Teuchos::rcp(new Epetra_MultiVector(map, number_of_ion_exchange_sites));
  } else {
    ion_exchange_sites_ = Teuchos::null;
    ion_exchange_ref_cation_conc_ = Teuchos::null;
  }

  if (number_of_sorption_sites > 0) {
    sorption_sites_ = Teuchos::rcp(new Epetra_MultiVector(map, number_of_sorption_sites));
    surface_complex_free_site_conc_ =
      Teuchos::rcp(new Epetra_MultiVector(map, number_of_sorption_sites));
  } else {
    sorption_sites_ = Teuchos::null;
    surface_complex_free_site_conc_ = Teuchos::null;
  }

  if (using_sorption_isotherms) {
    isotherm_kd_ = Teuchos::rcp(new Epetra_MultiVector(map, number_of_components));
    isotherm_freundlich_n_ = Teuchos::rcp(new Epetra_MultiVector(map, number_of_components));
    isotherm_langmuir_b_ = Teuchos::rcp(new Epetra_MultiVector(map, number_of_components));
  } else {
    isotherm_kd_ = Teuchos::null;
    isotherm_freundlich_n_ = Teuchos::null;
    isotherm_langmuir_b_ = Teuchos::null;
  }
}

void
chemistry_data::store(Teuchos::RCP<Epetra_MultiVector> free_ion_concentrations,
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
                      Teuchos::RCP<Epetra_MultiVector> isotherm_langmuir_b)
{
  if (free_ion_concentrations != Teuchos::null && free_ion_concentrations_ != Teuchos::null) {
    *free_ion_concentrations_ = *free_ion_concentrations;
  }
  if (primary_activity_coeff != Teuchos::null && primary_activity_coeff_ != Teuchos::null) {
    *primary_activity_coeff_ = *primary_activity_coeff;
  }
  if (secondary_activity_coeff != Teuchos::null && secondary_activity_coeff_ != Teuchos::null) {
    *secondary_activity_coeff_ = *secondary_activity_coeff;
  }
  if (mineral_volume_fractions != Teuchos::null && mineral_volume_fractions_ != Teuchos::null) {
    *mineral_volume_fractions_ = *mineral_volume_fractions;
  }
  if (mineral_specific_surface_area != Teuchos::null &&
      mineral_specific_surface_area_ != Teuchos::null) {
    *mineral_specific_surface_area_ = *mineral_specific_surface_area;
  }
  if (total_sorbed != Teuchos::null && total_sorbed_ != Teuchos::null) {
    *total_sorbed_ = *total_sorbed;
  }
  if (sorption_sites != Teuchos::null && sorption_sites_ != Teuchos::null) {
    *sorption_sites_ = *sorption_sites;
  }
  if (surface_complex_free_site_conc != Teuchos::null &&
      surface_complex_free_site_conc_ != Teuchos::null) {
    *surface_complex_free_site_conc_ = *surface_complex_free_site_conc;
  }
  if (ion_exchange_sites != Teuchos::null && ion_exchange_sites_ != Teuchos::null) {
    *ion_exchange_sites_ = *ion_exchange_sites;
  }
  if (ion_exchange_ref_cation_conc != Teuchos::null &&
      ion_exchange_ref_cation_conc_ != Teuchos::null) {
    *ion_exchange_ref_cation_conc_ = *ion_exchange_ref_cation_conc;
  }
  if (isotherm_kd != Teuchos::null && isotherm_kd_ != Teuchos::null) {
    *isotherm_kd_ = *isotherm_kd;
  }
  if (isotherm_freundlich_n != Teuchos::null && isotherm_freundlich_n_ != Teuchos::null) {
    *isotherm_freundlich_n_ = *isotherm_freundlich_n;
  }
  if (isotherm_langmuir_b != Teuchos::null && isotherm_langmuir_b_ != Teuchos::null) {
    *isotherm_langmuir_b_ = *isotherm_langmuir_b;
  }
}


void
chemistry_data::retrieve(Teuchos::RCP<Epetra_MultiVector> free_ion_concentrations,
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
                         Teuchos::RCP<Epetra_MultiVector> isotherm_langmuir_b)
{
  if (free_ion_concentrations != Teuchos::null && free_ion_concentrations_ != Teuchos::null) {
    *free_ion_concentrations = *free_ion_concentrations_;
  }
  if (primary_activity_coeff != Teuchos::null && primary_activity_coeff_ != Teuchos::null) {
    *primary_activity_coeff = *primary_activity_coeff_;
  }
  if (secondary_activity_coeff != Teuchos::null && secondary_activity_coeff_ != Teuchos::null) {
    *secondary_activity_coeff = *secondary_activity_coeff_;
  }
  if (mineral_volume_fractions != Teuchos::null && mineral_volume_fractions_ != Teuchos::null) {
    *mineral_volume_fractions = *mineral_volume_fractions_;
  }
  if (mineral_specific_surface_area != Teuchos::null &&
      mineral_specific_surface_area_ != Teuchos::null) {
    *mineral_specific_surface_area = *mineral_specific_surface_area_;
  }
  if (total_sorbed != Teuchos::null && total_sorbed_ != Teuchos::null) {
    *total_sorbed = *total_sorbed_;
  }
  if (sorption_sites != Teuchos::null && sorption_sites_ != Teuchos::null) {
    *sorption_sites = *sorption_sites_;
  }
  if (surface_complex_free_site_conc != Teuchos::null &&
      surface_complex_free_site_conc_ != Teuchos::null) {
    *surface_complex_free_site_conc = *surface_complex_free_site_conc_;
  }
  if (ion_exchange_sites != Teuchos::null && ion_exchange_sites_ != Teuchos::null) {
    *ion_exchange_sites = *ion_exchange_sites_;
  }
  if (ion_exchange_ref_cation_conc != Teuchos::null &&
      ion_exchange_ref_cation_conc_ != Teuchos::null) {
    *ion_exchange_ref_cation_conc = *ion_exchange_ref_cation_conc_;
  }
  if (isotherm_kd != Teuchos::null && isotherm_kd_ != Teuchos::null) {
    *isotherm_kd = *isotherm_kd_;
  }
  if (isotherm_freundlich_n != Teuchos::null && isotherm_freundlich_n_ != Teuchos::null) {
    *isotherm_freundlich_n = *isotherm_freundlich_n_;
  }
  if (isotherm_langmuir_b != Teuchos::null && isotherm_langmuir_b_ != Teuchos::null) {
    *isotherm_langmuir_b = *isotherm_langmuir_b_;
  }
}
