/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#ifdef ALQUIMIA_ENABLED

#ifdef _OPENMP
#include "omp.h"
#endif

#include <AlquimiaHelper_Structured.H>
#include <PorousMedia.H> // To get atm -> Pa conversion factor

#include <cmath>
#include <Utility.H>
#include <ParallelDescriptor.H>

static Real small_primary_free_ion_guess = 1.e-30;

AlquimiaHelper_Structured::AlquimiaHelper_Structured(Amanzi::AmanziChemistry::ChemistryEngine* _engine)
  : engine(_engine),   alq_sizes(engine->Sizes())
{
  using_sorption = false;
  using_isotherms = false;

  engine->GetPrimarySpeciesNames(primarySpeciesNames);
  engine->GetMineralNames(mineralNames);
  engine->GetSurfaceSiteNames(surfSiteNames);
  engine->GetIonExchangeNames(ionExchangeNames);
  Nimmobile = engine->NumSorbedSpecies();
  Nmobile = engine->NumPrimarySpecies();
  Nminerals = mineralNames.size();
  Nisotherms = engine->NumIsothermSpecies();
  NionExchange = ionExchangeNames.size();
  NsorptionSites = surfSiteNames.size();
  NfreeIonSpecies = engine->NumFreeIonSpecies();

  NauxInts = alq_sizes.num_aux_integers;
  NauxDoubles = alq_sizes.num_aux_doubles;

  if (Nimmobile > 0) {
    BL_ASSERT(Nimmobile == Nmobile);
    using_sorption = true;
  }

  if (NionExchange > 0) {
    BL_ASSERT(Nimmobile > 0);
    using_sorption = true;
  }

  if (Nisotherms > 0) {
    BL_ASSERT(Nimmobile > 0);
    BL_ASSERT(Nisotherms == Nmobile);
    using_isotherms = true;
    using_sorption = true;
  }

  SetupAuxVariables();

  if (NauxInts > 0) {
    int ndigits_ints = std::log(NauxInts + 1) + 1;
    for (int i=0; i<NauxInts; ++i) {
      const std::string label=BoxLib::Concatenate("Auxiliary_Integers_",i,ndigits_ints);
      int n = aux_chem_variables.size();
      aux_chem_variables[label] = n;
    }
  }

  if (NauxDoubles > 0) {
    int ndigits_doubles = std::log(NauxDoubles + 1) + 1;
    for (int i=0; i<NauxDoubles; ++i) {
      const std::string label=BoxLib::Concatenate("Auxiliary_Doubles_",i,ndigits_doubles);
      int n = aux_chem_variables.size();
      aux_chem_variables[label] = n;
    }
  }

  int nthreads = 1;
#ifdef _OPENMP
  nthreads = omp_get_max_threads();
#endif
  alquimia_properties.resize(nthreads);
  alquimia_state.resize(nthreads);
  alquimia_aux_in.resize(nthreads);
  alquimia_aux_out.resize(nthreads);


  for (int ithread = 0; ithread < nthreads; ithread++) {
    engine->InitState(alquimia_properties[ithread],
                      alquimia_state[ithread],
                      alquimia_aux_in[ithread],
                      alquimia_aux_out[ithread]);
  }

  bool dump = ParallelDescriptor::IOProcessor();
  if (dump) {
    std::cout << "AlquimiaHelper_Structured:" << std::endl;
    std::cout << "  Primary Species Names: " << std::endl;
    for (int i=0; i<primarySpeciesNames.size(); ++i) {
      std::cout << "      " << primarySpeciesNames[i] << std::endl;
    }
    std::cout << "  Mineral Names:" << std::endl;
    for (int i=0; i<Nminerals; ++i) {
      std::cout << "      " << mineralNames[i] << std::endl;
    }
    std::cout << "  Surface Site Names:" << std::endl;
    for (int i=0; i<surfSiteNames.size(); ++i) {
      std::cout << "      " << surfSiteNames[i] << std::endl;
    }
    std::cout << "  Ion Exchange Names:" << std::endl;
    for (int i=0; i<ionExchangeNames.size(); ++i) {
      std::cout << "      " << ionExchangeNames[i] << std::endl;
    }
    std::cout << "  Isotherm Species Names:" << std::endl;
    for (int i=0; i<Nisotherms; ++i) {
      std::cout << "      " << primarySpeciesNames[i] << std::endl;
    }
    std::cout << "  Number of sorbed species: " << Nimmobile << std::endl;
    std::cout << "  Number of aux ints: " << NauxInts << std::endl;
    std::cout << "  Number of aux doubles: " << NauxDoubles << std::endl;
    std::cout << "  Number of free ion species: " << NfreeIonSpecies << std::endl;

    std::cout << "  AuxChemMap: " << std::endl;
    Array<std::string> tmp(aux_chem_variables.size());
    for (std::map<std::string,int>::const_iterator it=aux_chem_variables.begin(); it!=aux_chem_variables.end(); ++it) {
      tmp[it->second] = it->first;
    }
    for (int i=0; i<tmp.size(); ++i) {
      std::cout << "    " << i << ":  " << tmp[i] << std::endl;
    }
  }

}

AlquimiaHelper_Structured::~AlquimiaHelper_Structured()
{
  int num_threads = 1;
#ifdef _OPENMP
  num_threads = omp_get_max_threads();
#endif
  for (int ithread = 0; ithread < num_threads; ithread++) {
    engine->FreeState(alquimia_properties[ithread],
                      alquimia_state[ithread],
                      alquimia_aux_in[ithread],
                      alquimia_aux_out[ithread]);
  }
}

static int
GetIndex(const std::map<std::string,int>& label_map,
	 const std::string&               key)
{
  std::map<std::string,int>::const_iterator it = label_map.find(key);
  if (it == label_map.end()) {
    BoxLib::Abort("Unable to find key in label map");
  }
  return it->second;
}


void
AlquimiaHelper_Structured::BL_to_Alquimia(const FArrayBox& aqueous_saturation,       int sSat,
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
                                          AlquimiaAuxiliaryOutputData& aux_output)
{
  mat_props.volume = volume(iv,sVol);
  mat_props.saturation = aqueous_saturation(iv,sSat);

  chem_state.water_density = water_density;
  chem_state.porosity = porosity(iv,sPhi);
  chem_state.temperature = temperature;
  chem_state.aqueous_pressure = aqueous_pressure(iv,sPress) * BL_ONEATM;

  for (int i=0; i<Nmobile; ++i) {
    chem_state.total_mobile.data[i] = std::max(0.,primary_species_mobile(iv,sPrimMob+i));
  }

  if (using_sorption) {
    for (int i=0; i<primarySpeciesNames.size(); ++i) {
      const std::string label=primarySpeciesNames[i] + "_Sorbed_Concentration";
      //chem_state.total_immobile.data[i] = aux_data(iv,aux_chem_variables[label]);
      int index = GetIndex(aux_chem_variables,label);
      chem_state.total_immobile.data[i] = std::max(1.e-20,aux_data(iv,index));
    }
  }

  if (Nminerals > 0) {
    for (int i=0; i<Nminerals; ++i) {
      const std::string label=mineralNames[i] + "_Volume_Fraction";
      int index = GetIndex(aux_chem_variables,label);
      chem_state.mineral_volume_fraction.data[i] = aux_data(iv,index);
    }
    for (int i=0; i<mineralNames.size(); ++i) {
      const std::string label=mineralNames[i] + "_Specific_Surface_Area";
      int index = GetIndex(aux_chem_variables,label);
      chem_state.mineral_specific_surface_area.data[i] = aux_data(iv,index);
    }
  }

  if (NionExchange > 0) {
    int ndigIES = std::log(NionExchange+1);
    for (int i=0; i<NionExchange; ++i) {
      const std::string label = BoxLib::Concatenate("Ion_Exchange_Site_Density_",i,ndigIES);
      int index = GetIndex(aux_chem_variables,label);
      chem_state.cation_exchange_capacity.data[i] = aux_data(iv,index);
    }

    for (int i=0; i<NionExchange; ++i) {
      const std::string label = BoxLib::Concatenate("Ion_Exchange_Reference_Cation_Concentration_",i,ndigIES);
      //int index = GetIndex(aux_chem_variables,label);
      //chem_state.ion_exchange_ref_cation_conc.data[i] = aux_data(iv,index);
    }
  }

  if (NsorptionSites > 0) {
    for (int i=0; i<surfSiteNames.size(); ++i) {
      const std::string label=surfSiteNames[i] + "_Surface_Site_Density";
      int index = GetIndex(aux_chem_variables,label);
      chem_state.surface_site_density.data[i] = aux_data(iv,index);
    }
  }

  for (int i=0; i<primarySpeciesNames.size(); ++i) {
    const std::string label=primarySpeciesNames[i] + "_Activity_Coefficient";
    int index = GetIndex(aux_chem_variables,label);
    aux_output.primary_activity_coeff.data[i] = aux_data(iv,index);

    //aux_output.primary_activity_coeff.data[i] = 0;

  }

  if (NfreeIonSpecies > 0) {
    for (int i=0; i<primarySpeciesNames.size(); ++i) {
      const std::string label=primarySpeciesNames[i] + "_Free_Ion_Guess";
      int index = GetIndex(aux_chem_variables,label);
      aux_output.primary_free_ion_concentration.data[i] = aux_data(iv,index);

      //aux_output.primary_free_ion_concentration.data[i] = 0;

    }
  }

  if (using_isotherms) {
    for (int i=0; i<Nisotherms; ++i) {
      const std::string label=primarySpeciesNames[i] + "_Isotherm_Kd";
      int index = GetIndex(aux_chem_variables,label);
      mat_props.isotherm_kd.data[i] = aux_data(iv,index);
    }
    for (int i=0; i<Nisotherms; ++i) {
      const std::string label=primarySpeciesNames[i] + "_Isotherm_Freundlich_n";
      int index = GetIndex(aux_chem_variables,label);
      mat_props.freundlich_n.data[i] = aux_data(iv,index);
    }
    for (int i=0; i<Nisotherms; ++i) {
      const std::string label=primarySpeciesNames[i] + "_Isotherm_Langmuir_b";
      int index = GetIndex(aux_chem_variables,label);
      mat_props.langmuir_b.data[i] = aux_data(iv,index);
    }
  }

  if (NauxInts > 0) {
    int ndigits_ints = std::log(NauxInts + 1) + 1;
    for (int i=0; i<NauxInts; ++i) {
      const std::string label=BoxLib::Concatenate("Auxiliary_Integers_",i,ndigits_ints);
      int index = GetIndex(aux_chem_variables,label);
      // safe way to cast integers, in case something is wrong with initialization (e.g. val=NaN)
      double val =  aux_data(iv,index);
      if (val == val) {
          aux_input.aux_ints.data[i] = (int) val;
      }
      else {
          aux_input.aux_ints.data[i] = 0;
      }
    }
  }
  if (NauxDoubles > 0) {
    int ndigits_doubles = std::log(NauxDoubles + 1) + 1;
    for (int i=0; i<NauxDoubles; ++i) {
      const std::string label=BoxLib::Concatenate("Auxiliary_Doubles_",i,ndigits_doubles);
      int index = GetIndex(aux_chem_variables,label);
      aux_input.aux_doubles.data[i] = aux_data(iv,index);
    }

    if (NfreeIonSpecies > 0) {
      for (int i=0; i<primarySpeciesNames.size(); ++i) {
	const std::string label=BoxLib::Concatenate("Auxiliary_Doubles_",i,ndigits_doubles);
	aux_input.aux_doubles.data[i] = std::max(aux_input.aux_doubles.data[i],
						 small_primary_free_ion_guess);
      }
    }

  }
}

void
AlquimiaHelper_Structured::EnforceCondition(FArrayBox& primary_species_mobile,   int sPrimMob,
                                            FArrayBox& auxiliary_data, Real water_density, Real temperature,
                                            const Box& box, const std::string& condition_name, Real time,
					    int chem_verbose)
{
#if (BL_SPACEDIM == 3) && defined(_OPENMP)
#pragma omp parallel for schedule(dynamic,1)
#endif

  int thread_outer_lo = box.smallEnd()[BL_SPACEDIM-1];
  int thread_outer_hi = box.bigEnd()[BL_SPACEDIM-1];
  bool chem_ok = true;

  FArrayBox dumS(box,1); dumS.setVal(1); int sDumS=0;
  FArrayBox dumP(box,1); dumP.setVal(101325); int sDumP=0;
  FArrayBox dumPhi(box,1); dumPhi.setVal(1); int sDumPhi=0;
  FArrayBox dumV(box,1); dumV.setVal(1); int sDumV=0;

  for (int tli=thread_outer_lo; tli<=thread_outer_hi && chem_ok; tli++) {
#if (BL_SPACEDIM == 3) && defined(_OPENMP)
    int threadid = omp_get_thread_num();
#else
    int threadid = 0;
#endif

    Box thread_box(box);
    thread_box.setSmall(BL_SPACEDIM-1,tli);
    thread_box.setBig(BL_SPACEDIM-1,tli);

    for (IntVect iv=thread_box.smallEnd(), End=thread_box.bigEnd(); iv<=End; thread_box.next(iv)) {

      BL_to_Alquimia(dumS,sDumS,dumP,sDumP,dumPhi,sDumPhi,dumV,sDumV,
                     primary_species_mobile,sPrimMob,
		     auxiliary_data,iv,water_density,temperature,
                     alquimia_properties[threadid],
                     alquimia_state[threadid],
                     alquimia_aux_in[threadid],
                     alquimia_aux_out[threadid]);

      if (NauxDoubles > 0) {
	for (int i=0; i<NauxDoubles; ++i) {
          alquimia_aux_in[0].aux_doubles.data[i] = 0;
	}
      }
      engine->EnforceCondition(condition_name,time,
                               alquimia_properties[threadid],
                               alquimia_state[threadid],
                               alquimia_aux_in[threadid],
                               alquimia_aux_out[threadid]);

      Alquimia_to_BL(primary_species_mobile,   sPrimMob,
                     auxiliary_data, iv,
                     alquimia_properties[threadid],
                     alquimia_state[threadid],
                     alquimia_aux_in[threadid],
                     alquimia_aux_out[threadid]);
    }
  }
}

void
AlquimiaHelper_Structured::Advance(const FArrayBox& aqueous_saturation,       int sSat,
                                   const FArrayBox& aqueous_pressure,         int sPress,
                                   const FArrayBox& porosity,                 int sPhi,
                                   const FArrayBox& volume,                   int sVol,
                                   FArrayBox&       primary_species_mobile,   int sPrimMob,
                                   FArrayBox&       fcnCnt,                   int sFunc,
                                   FArrayBox&       aux_data, Real water_density, Real temperature,
                                   const Box& box, Real dt, int chem_verbose)
{
#if (BL_SPACEDIM == 3) && defined(_OPENMP)
#pragma omp parallel for schedule(dynamic,1)
#endif

  int thread_outer_lo = box.smallEnd()[BL_SPACEDIM-1];
  int thread_outer_hi = box.bigEnd()[BL_SPACEDIM-1];
  bool chem_ok = true;
  for (int tli=thread_outer_lo; tli<=thread_outer_hi && chem_ok; tli++) {
#if (BL_SPACEDIM == 3) && defined(_OPENMP)
    int threadid = omp_get_thread_num();
#else
    int threadid = 0;
#endif

    Box thread_box(box);
    thread_box.setSmall(BL_SPACEDIM-1,tli);
    thread_box.setBig(BL_SPACEDIM-1,tli);

    for (IntVect iv=thread_box.smallEnd(), End=thread_box.bigEnd(); iv<=End; thread_box.next(iv)) {

      BL_to_Alquimia(aqueous_saturation,       sSat,
                     aqueous_pressure,         sPress,
                     porosity,                 sPhi,
                     volume,                   sVol,
                     primary_species_mobile,   sPrimMob,
                     aux_data, iv, water_density, temperature,
                     alquimia_properties[threadid],
                     alquimia_state[threadid],
                     alquimia_aux_in[threadid],
                     alquimia_aux_out[threadid]);

      if (ParallelDescriptor::IOProcessor() && chem_verbose>0) {
	std::cout << "************* GOING INTO engine->Advance: " << std::endl;
	std::cout << "iv: " << iv << std::endl;
	std::cout << "dt: " << dt << std::endl;
	DumpAlquimiaStructures(std::cout,alquimia_properties[threadid],
			       alquimia_state[threadid],
			       alquimia_aux_in[threadid],
			       alquimia_aux_out[threadid]);
      }

      int newton_iters;
      int natural_id;
      engine->Advance(dt,alquimia_properties[threadid],alquimia_state[threadid],
                      alquimia_aux_in[threadid],alquimia_aux_out[threadid],newton_iters,natural_id);

      if (ParallelDescriptor::IOProcessor() && chem_verbose>0) {
	std::cout << "************* COMING OUT OF engine->Advance: " << std::endl;
	std::cout << "iv: " << iv << std::endl;
	std::cout << "dt: " << dt << std::endl;
	DumpAlquimiaStructures(std::cout,alquimia_properties[threadid],
			       alquimia_state[threadid],
			       alquimia_aux_in[threadid],
			       alquimia_aux_out[threadid]);
      }

      Alquimia_to_BL(primary_species_mobile,   sPrimMob,
                     aux_data, iv,
                     alquimia_properties[threadid],
                     alquimia_state[threadid],
                     alquimia_aux_in[threadid],
                     alquimia_aux_out[threadid]);

      fcnCnt(iv,sFunc) = newton_iters;
    }
  }
}

void
AlquimiaHelper_Structured::Alquimia_to_BL(FArrayBox& primary_species_mobile,   int sPrimMob,
                                          FArrayBox&                   aux_data,
                                          const IntVect&               iv,
                                          AlquimiaProperties&  mat_props,
                                          AlquimiaState&               chem_state,
                                          AlquimiaAuxiliaryData&       aux_input,
                                          AlquimiaAuxiliaryOutputData& aux_output)
{
  for (int i=0; i<Nmobile; ++i) {
    primary_species_mobile(iv,sPrimMob+i) = chem_state.total_mobile.data[i];
  }

  if (using_sorption) {
    for (int i=0; i<primarySpeciesNames.size(); ++i) {
      const std::string label=primarySpeciesNames[i] + "_Sorbed_Concentration";
      aux_data(iv,aux_chem_variables[label]) = chem_state.total_immobile.data[i];
    }
  }

#if 0
  for (int i=0; i<primarySpeciesNames.size(); ++i) {
    const std::string label=primarySpeciesNames[i] + "_Activity_Coefficient";
    aux_data(iv,aux_chem_variables[label]) = aux_output.primary_activity_coeff.data[i];
  }
#endif

  if (NfreeIonSpecies > 0) {
    for (int i=0; i<primarySpeciesNames.size(); ++i) {
      const std::string label=primarySpeciesNames[i] + "_Free_Ion_Guess";
      aux_data(iv,aux_chem_variables[label]) = std::max(aux_output.primary_free_ion_concentration.data[i],
							small_primary_free_ion_guess);
    }
  }

  if (Nminerals > 0) {
    for (int i=0; i<Nminerals; ++i) {
      const std::string label=mineralNames[i] + "_Volume_Fraction";
      aux_data(iv,aux_chem_variables[label]) = chem_state.mineral_volume_fraction.data[i];
    }
    for (int i=0; i<mineralNames.size(); ++i) {
      const std::string label=mineralNames[i] + "_Specific_Surface_Area";
      aux_data(iv,aux_chem_variables[label]) = chem_state.mineral_specific_surface_area.data[i];
    }
  }

  if (NsorptionSites > 0) {
    for (int i=0; i<surfSiteNames.size(); ++i) {
      const std::string label=surfSiteNames[i] + "_Surface_Site_Density";
      aux_data(iv,aux_chem_variables[label]) = chem_state.surface_site_density.data[i];
    }
  }

  if (NionExchange > 0) {
    int ndigIES = std::log(NionExchange+1);
    for (int i=0; i<NionExchange; ++i) {
      const std::string label = BoxLib::Concatenate("Ion_Exchange_Site_Density_",i,ndigIES);
      aux_data(iv,aux_chem_variables[label]) = chem_state.cation_exchange_capacity.data[i];
    }
    for (int i=0; i<NionExchange; ++i) {
      const std::string label = BoxLib::Concatenate("Ion_Exchange_Reference_Cation_Concentration_",i,ndigIES);
      //aux_data(iv,aux_chem_variables[label]) = chem_state.ion_exchange_ref_cation_conc.data[i];
    }
  }

#if 0
  if (using_isotherms) {
    for (int i=0; i<Nisotherms; ++i) {
      const std::string label=primarySpeciesNames[i] + "_Isotherm_Kd";
      aux_data(iv,aux_chem_variables[label]) = mat_props.isotherm_kd.data[i];
    }
    for (int i=0; i<Nisotherms; ++i) {
      const std::string label=primarySpeciesNames[i] + "_Isotherm_Freundlich_n";
      aux_data(iv,aux_chem_variables[label]) = mat_props.freundlich_n.data[i];
    }
    for (int i=0; i<Nisotherms; ++i) {
      const std::string label=primarySpeciesNames[i] + "_Isotherm_Langmuir_b";
      aux_data(iv,aux_chem_variables[label]) = mat_props.langmuir_b.data[i];
    }
  }
#endif

  for (int i=0; i<NauxInts; ++i) {
    int ndigits_ints = std::log(NauxInts + 1) + 1;
    const std::string label=BoxLib::Concatenate("Auxiliary_Integers_",i,ndigits_ints);
    aux_data(iv,aux_chem_variables[label]) = (int) aux_input.aux_ints.data[i];
  }
  for (int i=0; i<NauxDoubles; ++i) {
    int ndigits_doubles = std::log(NauxDoubles + 1) + 1;
    const std::string label=BoxLib::Concatenate("Auxiliary_Doubles_",i,ndigits_doubles);
    aux_data(iv,aux_chem_variables[label]) = aux_input.aux_doubles.data[i];
  }

  if (NfreeIonSpecies > 0) {
    for (int i=0; i<primarySpeciesNames.size(); ++i) {
      int ndigits_doubles = std::log(NauxDoubles + 1) + 1;
      const std::string label=BoxLib::Concatenate("Auxiliary_Doubles_",i,ndigits_doubles);
      //aux_data(iv,aux_chem_variables[label]) = aux_input.aux_doubles.data[i];
      aux_data(iv,aux_chem_variables[label]) = std::max(aux_data(iv,aux_chem_variables[label]),
							small_primary_free_ion_guess);
    }
  }
}


void
AlquimiaHelper_Structured::DumpAlquimiaStructures(std::ostream&                os,
						  AlquimiaProperties&  mat_props,
						  AlquimiaState&               chem_state,
						  AlquimiaAuxiliaryData&       aux_input,
						  AlquimiaAuxiliaryOutputData& aux_output)
{
  os << "Volume " << mat_props.volume << std::endl;
  os << "Saturation " << mat_props.saturation << std::endl;
  os << "Water density " << chem_state.water_density << std::endl;
  os << "Porosity " << chem_state.porosity << std::endl;
  os << "Temperature " << chem_state.temperature << std::endl;
  os << "Aqueous Pressure " << chem_state.aqueous_pressure << std::endl;
  for (int i=0; i<Nmobile; ++i) {
    const std::string label=primarySpeciesNames[i] + "_Primary_Species_Mobile";
    os << label << "  " << chem_state.total_mobile.data[i] << std::endl;
  }

  if (using_sorption) {
    for (int i=0; i<Nmobile; ++i) {
      const std::string label=primarySpeciesNames[i] + "_Primary_Species_Sorbed";
      os << label << "  " << chem_state.total_immobile.data[i] << std::endl;
    }
  }

  if (Nminerals > 0) {
    for (int i=0; i<Nminerals; ++i) {
      const std::string label=mineralNames[i] + "_Volume_Fraction";
      os << label << "  " << chem_state.mineral_volume_fraction.data[i] << std::endl;
    }
    for (int i=0; i<Nminerals; ++i) {
      const std::string label=mineralNames[i] + "_Specific_Surface_Area";
      os << label << "  " << chem_state.mineral_specific_surface_area.data[i] << std::endl;
    }
  }

  if (NionExchange > 0) {
    int ndigIES = std::log(NionExchange+1);
    for (int i=0; i<NionExchange; ++i) {
      const std::string label = BoxLib::Concatenate("Ion_Exchange_Site_Density_",i,ndigIES);
      os << label << "  " << chem_state.cation_exchange_capacity.data[i] << std::endl;
    }
  }

  if (NsorptionSites > 0) {
    for (int i=0; i<surfSiteNames.size(); ++i) {
      const std::string label=surfSiteNames[i] + "_Surface_Site_Density";
      os << label << "  " << chem_state.surface_site_density.data[i] << std::endl;
    }
  }

  for (int i=0; i<primarySpeciesNames.size(); ++i) {
    const std::string label=primarySpeciesNames[i] + "_Activity_Coefficient";
    os << label << "  " << aux_output.primary_activity_coeff.data[i] << std::endl;
  }

  if (NfreeIonSpecies > 0) {
    for (int i=0; i<primarySpeciesNames.size(); ++i) {
      const std::string label=primarySpeciesNames[i] + "_Free_Ion_Guess";
      os << label << "  " << aux_output.primary_free_ion_concentration.data[i] << std::endl;
    }
  }

  if (using_isotherms) {
    for (int i=0; i<Nisotherms; ++i) {
      const std::string label=primarySpeciesNames[i] + "_Isotherm_Kd";
      os << label << "  " << mat_props.isotherm_kd.data[i] << std::endl;
    }
    for (int i=0; i<Nisotherms; ++i) {
      const std::string label=primarySpeciesNames[i] + "_Isotherm_Freundlich_n";
      os << label << "  " << mat_props.freundlich_n.data[i] << std::endl;
    }
    for (int i=0; i<Nisotherms; ++i) {
      const std::string label=primarySpeciesNames[i] + "_Isotherm_Langmuir_b";
      os << label << "  " << mat_props.langmuir_b.data[i] << std::endl;
    }
  }

  if (NauxInts > 0) {
    int ndigits_ints = std::log(NauxInts + 1) + 1;
    for (int i=0; i<NauxInts; ++i) {
      const std::string label=BoxLib::Concatenate("Auxiliary_Integers_",i,ndigits_ints);
      os << label << "  " << aux_input.aux_ints.data[i] << std::endl;
    }
  }
  if (NauxDoubles > 0) {
    int ndigits_doubles = std::log(NauxDoubles + 1) + 1;
    for (int i=0; i<NauxDoubles; ++i) {
      const std::string label=BoxLib::Concatenate("Auxiliary_Doubles_",i,ndigits_doubles);
      os << label << "  " << aux_input.aux_doubles.data[i] << std::endl;
    }
  }
}

#endif
