/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#ifndef _AlquimiaHelper_Structured_H_
#define _AlquimiaHelper_Structured_H_

#include <ChemistryHelper_Structured.H>

#ifdef ALQUIMIA_ENABLED
#include <ChemistryEngine.hh>

class AlquimiaHelper_Structured
  : public ChemistryHelper_Structured
{
public:
  AlquimiaHelper_Structured(Amanzi::AmanziChemistry::ChemistryEngine* engine);

  virtual ~AlquimiaHelper_Structured();

  virtual void Advance(const FArrayBox& aqueous_saturation,       int sSat,
                       const FArrayBox& aqueous_pressure,         int sPress,
                       const FArrayBox& porosity,                 int sPhi,
                       const FArrayBox& volume,                   int sVol,
                       FArrayBox&       primary_species_mobile,   int sPrimMob,
                       FArrayBox&       fcnCnt,                   int sFunc,
                       FArrayBox&       aux_data, Real water_density, Real temperature,
                       const Box& box, Real dt, int chem_verbose);

  virtual void EnforceCondition(FArrayBox& primary_species_mobile,   int sPrimMob,
                                FArrayBox& auxiliary_data, Real water_density, Real temperature,
                                const Box& box, const std::string& condition_name, Real time,
				int chem_verbose);

protected:
  void BL_to_Alquimia(const FArrayBox& aqueous_saturation,       int sSat,
                      const FArrayBox& aqueous_pressure,         int sPress,
                      const FArrayBox& porosity,                 int sPhi,
                      const FArrayBox& volume,                   int sVol,
                      FArrayBox&       primary_species_mobile,   int sPrimMob,
                      FArrayBox&                   aux_data,
                      const IntVect&               iv,
                      Real                         water_density,
                      Real                         temperature,
                      AlquimiaProperties&  mat_props,
                      AlquimiaState&               chem_state,
                      AlquimiaAuxiliaryData&       aux_input,
                      AlquimiaAuxiliaryOutputData& aux_output);

  void Alquimia_to_BL(FArrayBox& primary_species_mobile,   int sPrimMob,
                      FArrayBox&                   aux_data,
                      const IntVect&               iv,
                      AlquimiaProperties&  mat_props,
                      AlquimiaState&               chem_state,
                      AlquimiaAuxiliaryData&       aux_input,
                      AlquimiaAuxiliaryOutputData& aux_output);

  void DumpAlquimiaStructures(std::ostream&                os,
			      AlquimiaProperties&  mat_props,
			      AlquimiaState&               chem_state,
			      AlquimiaAuxiliaryData&       aux_input,
			      AlquimiaAuxiliaryOutputData& aux_output);

  Amanzi::AmanziChemistry::ChemistryEngine* engine;
  int nthreads;

  std::vector<AlquimiaState> alquimia_state;
  std::vector<AlquimiaProperties> alquimia_properties;
  std::vector<AlquimiaAuxiliaryData> alquimia_aux_in;
  std::vector<AlquimiaAuxiliaryOutputData> alquimia_aux_out;
  const AlquimiaSizes& alq_sizes;
  int NauxInts, NauxDoubles;
};

#endif

#endif
