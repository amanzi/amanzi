
#include <AlquimiaHelper_Structured.H>

#include <cmath>
#include <Utility.H>

AlquimiaHelper_Structured::AlquimiaHelper_Structured(Amanzi::AmanziChemistry::ChemistryEngine* _engine)
  : engine(_engine),   alq_sizes(engine->Sizes())

{
  engine->GetMineralNames(mineralNames);
  Nminerals = mineralNames.size();
  engine->GetSurfaceSiteNames(surfSiteNames);
  engine->GetIonExchangeNames(ionExchangeNames);
  engine->GetIsothermSpeciesNames(isothermSpeciesNames);
  Nimmobile = engine->NumSorbedSpecies();
  Nmobile = engine->NumPrimarySpecies() - Nimmobile;
  num_aux_ints = alq_sizes.num_aux_integers;
  ndigits_int = std::log(num_aux_ints + 1) + 1;
  num_aux_doubles = alq_sizes.num_aux_doubles;
  ndigits_doubles = std::log(num_aux_doubles + 1) + 1;
  num_free_ion_species = engine->NumFreeIonSpecies();
  ndigits_freeion = std::log(num_free_ion_species + 1) + 1;

  std::vector<std::string> primarySpeciesNames;
  engine->GetPrimarySpeciesNames(primarySpeciesNames);
  std::vector<std::string> mineralNames;
  engine->GetMineralNames(mineralNames);

  // Setup input maps
  aux_chem_variables.clear();
  for (int i=0; i<mineralNames.size(); ++i) {
    const std::string label=mineralNames[i] + "_Volume_Fraction"; 
    aux_chem_variables[label] = aux_chem_variables.size()-1;
  }
  for (int i=0; i<mineralNames.size(); ++i) {
    const std::string label=mineralNames[i] + "_Specific_Surface_Area"; 
    aux_chem_variables[label] = aux_chem_variables.size()-1;
  }
  std::vector<std::string> surfSiteNames; engine->GetSurfaceSiteNames(surfSiteNames);
  for (int i=0; i<surfSiteNames.size(); ++i) {
    const std::string label=surfSiteNames[i] + "_Surface_Site_Density"; 
    aux_chem_variables[label] = aux_chem_variables.size()-1;
  }
  std::vector<std::string> ionExchangeNames; engine->GetIonExchangeNames(ionExchangeNames);
  for (int i=0; i<ionExchangeNames.size(); ++i) {
    const std::string label=ionExchangeNames[i] + "_Ion_Exchange_Capacity"; 
    aux_chem_variables[label] = aux_chem_variables.size()-1;
  }
  std::vector<std::string> isothermSpeciesNames; engine->GetIsothermSpeciesNames(isothermSpeciesNames);
  for (int i=0; i<isothermSpeciesNames.size(); ++i) {
    const std::string label=isothermSpeciesNames[i] + "_Isotherm_Kd"; 
    aux_chem_variables[label] = aux_chem_variables.size()-1;
  }
  for (int i=0; i<isothermSpeciesNames.size(); ++i) {
    const std::string label=isothermSpeciesNames[i] + "_Freundlich_n"; 
    aux_chem_variables[label] = aux_chem_variables.size()-1;
  }
  for (int i=0; i<isothermSpeciesNames.size(); ++i) {
    const std::string label=isothermSpeciesNames[i] + "_Langmuir_b"; 
    aux_chem_variables[label] = aux_chem_variables.size()-1;
  }
  int num_free_ion_species = engine->NumFreeIonSpecies();
  int ndigits = std::log(num_free_ion_species + 1) + 1;
  for (int i=0; i<num_free_ion_species; ++i) {
    const std::string label=BoxLib::Concatenate("Free_Ion_Species_",i,ndigits); 
    aux_chem_variables[label] = aux_chem_variables.size()-1;
  }
  
  const AlquimiaSizes& alq_sizes = engine->Sizes();
  int num_aux_ints = alq_sizes.num_aux_integers;
  int ndigits_ints = std::log(num_aux_ints + 1) + 1;
  for (int i=0; i<num_aux_ints; ++i) {
    const std::string label=BoxLib::Concatenate("Auxiliary_Integers_",i,ndigits_ints); 
    aux_chem_variables[label] = aux_chem_variables.size()-1;
  }
  int num_aux_doubles = alq_sizes.num_aux_doubles;
  int ndigits_doubles = std::log(num_aux_doubles + 1) + 1;
  for (int i=0; i<num_aux_doubles; ++i) {
    const std::string label=BoxLib::Concatenate("Auxiliary_Doubles_",i,ndigits_doubles); 
    aux_chem_variables[label] = aux_chem_variables.size()-1;
  }
  
  int nthreads = 1;
#ifdef _OPENMP
  nthreads = omp_get_max_threads();
#endif
  alquimia_material_properties.resize(nthreads);
  alquimia_state.resize(nthreads);
  alquimia_aux_in.resize(nthreads);
  alquimia_aux_out.resize(nthreads);


  for (int ithread = 0; ithread < nthreads; ithread++) {
    engine->InitState(alquimia_material_properties[ithread],
                      alquimia_state[ithread],
                      alquimia_aux_in[ithread],
                      alquimia_aux_out[ithread]);
  }

  bool dump = false;
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
    for (int i=0; i<isothermSpeciesNames.size(); ++i) {
      std::cout << "      " << isothermSpeciesNames[i] << std::endl;
    }
    std::cout << "  Number of sorbed species: " << Nimmobile << std::endl;
    std::cout << "  Number of aux ints: " << num_aux_ints << std::endl;
    std::cout << "  Number of aux doubles: " << num_aux_doubles << std::endl;
    std::cout << "  Number of free ion species: " << num_free_ion_species << std::endl;

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
  int nthreads = 1;
#ifdef _OPENMP
  nthreads = omp_get_max_threads();
#endif
  for (int ithread = 0; ithread < nthreads; ithread++) {
    engine->FreeState(alquimia_material_properties[ithread],
                      alquimia_state[ithread],
                      alquimia_aux_in[ithread],
                      alquimia_aux_out[ithread]);
  }
}

void
AlquimiaHelper_Structured::BL_to_Alquimia(const FArrayBox& aqueous_saturation,       int sSat,
                                          const FArrayBox& aqueous_pressure,         int sPress,
                                          const FArrayBox& porosity,                 int sPhi,
                                          const FArrayBox& volume,                   int sVol,
                                          FArrayBox&       primary_species_mobile,   int sPrimMob,
                                          FArrayBox&       primary_species_immobile, int sPrimImmob,
                                          FArrayBox&                   aux_data,
                                          bool                         auxilary_data_initialized,
                                          const IntVect&               iv,
                                          Real                         water_density,
                                          Real                         temperature,
                                          AlquimiaMaterialProperties&  mat_props,
                                          AlquimiaState&               chem_state,
                                          AlquimiaAuxiliaryData&       aux_input,
                                          AlquimiaAuxiliaryOutputData& aux_output)
{
  mat_props.volume = volume(iv,sVol);
  mat_props.saturation = aqueous_saturation(iv,sSat);

  chem_state.water_density = water_density;
  chem_state.porosity = porosity(iv,sPhi);
  chem_state.temperature = temperature;
  chem_state.aqueous_pressure = aqueous_pressure(iv,sPress);
  for (int i=0; i<Nmobile; ++i) {
    chem_state.total_mobile.data[i] = primary_species_mobile(iv,sPrimMob+i);
  }
  for (int i=0; i<Nimmobile; ++i) {
    chem_state.total_immobile.data[i] = primary_species_immobile(iv,sPrimImmob+i);
  }

  if (auxilary_data_initialized) {
    for (int i=0; i<Nminerals; ++i) {
      const std::string label=mineralNames[i] + "_Volume_Fraction"; 
      chem_state.mineral_volume_fraction.data[i] = aux_data(iv,aux_chem_variables[label]);
    }
    for (int i=0; i<mineralNames.size(); ++i) {
      const std::string label=mineralNames[i] + "_Specific_Surface_Area"; 
      chem_state.mineral_specific_surface_area.data[i] = aux_data(iv,aux_chem_variables[label]);
    }
    for (int i=0; i<surfSiteNames.size(); ++i) {
      const std::string label=surfSiteNames[i] + "_Surface_Site_Density"; 
      chem_state.surface_site_density.data[i] = aux_data(iv,aux_chem_variables[label]);
    }
    for (int i=0; i<ionExchangeNames.size(); ++i) {
      const std::string label=ionExchangeNames[i] + "_Ion_Exchange_Capacity"; 
      chem_state.cation_exchange_capacity.data[i] = aux_data(iv,aux_chem_variables[label]);
    }

    for (int i=0; i<isothermSpeciesNames.size(); ++i) {
      const std::string label=isothermSpeciesNames[i] + "_Isotherm_Kd"; 
      mat_props.isotherm_kd.data[i] = aux_data(iv,aux_chem_variables[label]);
    }
    for (int i=0; i<isothermSpeciesNames.size(); ++i) {
      const std::string label=isothermSpeciesNames[i] + "_Freundlich_n"; 
      mat_props.freundlich_n.data[i] = aux_data(iv,aux_chem_variables[label]);
    }
    for (int i=0; i<isothermSpeciesNames.size(); ++i) {
      const std::string label=isothermSpeciesNames[i] + "_Langmuir_b"; 
      mat_props.langmuir_b.data[i] = aux_data(iv,aux_chem_variables[label]);
    }

    for (int i=0; i<num_free_ion_species; ++i) {
      const std::string label=BoxLib::Concatenate("Free_Ion_Species_",i,ndigits_freeion); 
      aux_output.primary_free_ion_concentration.data[i] = aux_data(iv,aux_chem_variables[label]);
    }

    for (int i=0; i<num_aux_ints; ++i) {
      const std::string label=BoxLib::Concatenate("Auxiliary_Integers_",i,ndigits_int); 
      aux_input.aux_ints.data[i] = (int) aux_data(iv,aux_chem_variables[label]);
    }
    for (int i=0; i<num_aux_doubles; ++i) {
      const std::string label=BoxLib::Concatenate("Auxiliary_Doubles_",i,ndigits_doubles); 
      aux_input.aux_doubles.data[i] = aux_data(iv,aux_chem_variables[label]);
    }
  }
}

void
AlquimiaHelper_Structured::EnforceCondition(FArrayBox& primary_species_mobile,   int sPrimMob,
                                            FArrayBox& primary_species_immobile, int sPrimImmob,
                                            FArrayBox& auxiliary_data, bool initialize_auxiliary_data,
                                            const Box& box, const std::string& condition_name, Real time)
{
  bool auxiliary_data_initialized = !initialize_auxiliary_data;

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
  Real dumRho = 998.2;
  Real dumTemp = 295;
  
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
                     primary_species_mobile,sPrimMob,primary_species_immobile,sPrimImmob,
		     auxiliary_data,auxiliary_data_initialized,iv,dumRho,dumTemp,
                     alquimia_material_properties[threadid],
                     alquimia_state[threadid],
                     alquimia_aux_in[threadid],
                     alquimia_aux_out[threadid]);

      engine->EnforceCondition(condition_name,time,
                               alquimia_material_properties[threadid],
                               alquimia_state[threadid],
                               alquimia_aux_in[threadid],
                               alquimia_aux_out[threadid]);

      Alquimia_to_BL(primary_species_mobile,   sPrimMob,
                     primary_species_immobile, sPrimImmob,
                     auxiliary_data, initialize_auxiliary_data, iv,
                     alquimia_material_properties[threadid],
                     alquimia_state[threadid],
                     alquimia_aux_out[threadid]);
      // {
      //   Box cbox = BoxLib::adjCellLo(box,0,1); cbox.shift(0,1);
      //   FArrayBox tc(cbox,Nmobile);
      //   tc.copy(primary_species_mobile,sPrimMob,0,Nmobile);
      //   std::cout << "Constraint state" << std::endl;
      //   std::cout << tc << std::endl;
      //   FArrayBox ac(cbox,auxiliary_data.nComp());
      //   ac.copy(auxiliary_data);
      //   std::cout << ac << std::endl;
      //   FArrayBox vc(cbox,1); vc.copy(dumV,sDumV,0,1);
      //   FArrayBox sc(cbox,1); sc.copy(dumS,sDumS,0,1);
      //   FArrayBox pc(cbox,1); pc.copy(dumPhi,sDumPhi,0,1);
      //   std::cout << vc << std::endl;
      //   std::cout << sc << std::endl;
      //   std::cout << pc << std::endl;
      // }
    }
  }
}

void
AlquimiaHelper_Structured::Advance(const FArrayBox& aqueous_saturation,       int sSat,
                                   const FArrayBox& aqueous_pressure,         int sPress,
                                   const FArrayBox& porosity,                 int sPhi,
                                   const FArrayBox& volume,                   int sVol,
                                   FArrayBox&       primary_species_mobile,   int sPrimMob,
                                   FArrayBox&       primary_species_immobile, int sPrimImmob,
                                   FArrayBox&       fcnCnt,                   int sFunc,
                                   FArrayBox&       aux_data, Real water_density, Real temperature,
                                   const Box& box, Real dt)
{
  bool auxiliary_data_initialized = true;
  bool initialize_auxiliary_data = false;

#if (BL_SPACEDIM == 3) && defined(_OPENMP)
#pragma omp parallel for schedule(dynamic,1) 
#endif

  int thread_outer_lo = box.smallEnd()[BL_SPACEDIM-1];
  int thread_outer_hi = box.bigEnd()[BL_SPACEDIM-1];
  bool chem_ok = true;
  // {
  //   Box cbox = BoxLib::adjCellLo(box,0,1); cbox.shift(0,1);
  //   FArrayBox tc(cbox,Nmobile);
  //   tc.copy(primary_species_mobile,sPrimMob,0,Nmobile);
  //   std::cout << "Chem State before" << std::endl;
  //   std::cout << tc << std::endl;
  //   FArrayBox ac(cbox,aux_data.nComp());
  //   ac.copy(aux_data);
  //   std::cout << ac << std::endl;
  //   FArrayBox vc(cbox,1); vc.copy(volume,sVol,0,1);
  //   FArrayBox sc(cbox,1); sc.copy(aqueous_saturation,sSat,0,1);
  //   FArrayBox pc(cbox,1); pc.copy(porosity,sPhi,0,1);
  //   std::cout << vc << std::endl;
  //   std::cout << sc << std::endl;
  //   std::cout << pc << std::endl;
  //   std::cout << water_density << " " << temperature << " dt=" << dt << std::endl;
  // }
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
	
      // Fill local struct
      bool is_neg = false;

      BL_to_Alquimia(aqueous_saturation,       sSat,
                     aqueous_pressure,         sPress,
                     porosity,                 sPhi,
                     volume,                   sVol,
                     primary_species_mobile,   sPrimMob,
                     primary_species_immobile, sPrimImmob,
                     aux_data, auxiliary_data_initialized, iv, water_density, temperature,
                     alquimia_material_properties[threadid],
                     alquimia_state[threadid],
                     alquimia_aux_in[threadid],
                     alquimia_aux_out[threadid]);
        
      int newton_iters;
      engine->Advance(dt,alquimia_material_properties[threadid],alquimia_state[threadid],
                      alquimia_aux_in[threadid],alquimia_aux_out[threadid],newton_iters);
        
      Alquimia_to_BL(primary_species_mobile,   sPrimMob,
                     primary_species_immobile, sPrimImmob,
                     aux_data, initialize_auxiliary_data, iv,
                     alquimia_material_properties[threadid],
                     alquimia_state[threadid],
                     alquimia_aux_out[threadid]);

      fcnCnt(iv,sFunc) = newton_iters;
    }
  }
  // {
  //   Box cbox = BoxLib::adjCellLo(box,0,1); cbox.shift(0,1);
  //   FArrayBox tc(cbox,Nmobile);
  //   tc.copy(primary_species_mobile,sPrimMob,0,Nmobile);
  //   std::cout << "Chem State after" << std::endl;
  //   std::cout << tc << std::endl;
  //   FArrayBox ac(cbox,aux_data.nComp());
  //   ac.copy(aux_data);
  //   std::cout << ac << std::endl;
  //   FArrayBox vc(cbox,1); vc.copy(volume,sVol,0,1);
  //   FArrayBox sc(cbox,1); sc.copy(aqueous_saturation,sSat,0,1);
  //   FArrayBox pc(cbox,1); pc.copy(porosity,sPhi,0,1);
  //   std::cout << vc << std::endl;
  //   std::cout << sc << std::endl;
  //   std::cout << pc << std::endl;
  //   std::cout << water_density << " " << temperature << std::endl;
  // }
}

void
AlquimiaHelper_Structured::Alquimia_to_BL(FArrayBox& primary_species_mobile,   int sPrimMob,
                                          FArrayBox& primary_species_immobile, int sPrimImmob,
                                          FArrayBox&                   aux_data,
                                          bool                         initialize_auxiliary_data,
                                          const IntVect&               iv,
                                          AlquimiaMaterialProperties&  mat_props,
                                          AlquimiaState&               chem_state,
                                          AlquimiaAuxiliaryOutputData& aux_output)
{
  for (int i=0; i<Nmobile; ++i) {
    primary_species_mobile(iv,sPrimMob+i) = chem_state.total_mobile.data[i];
  }
  for (int i=0; i<Nimmobile; ++i) {
    primary_species_immobile(iv,sPrimImmob+i) = chem_state.total_immobile.data[i];
  }

  if (initialize_auxiliary_data) {

    for (int n=0; n<aux_data.nComp(); ++n) {
      aux_data(iv,n) = 1.e-8;
    }

    for (int i=0; i<Nminerals; ++i) {
      const std::string label=mineralNames[i] + "_Volume_Fraction"; 
      aux_data(iv,aux_chem_variables[label]) = chem_state.mineral_volume_fraction.data[i];
    }
    for (int i=0; i<mineralNames.size(); ++i) {
      const std::string label=mineralNames[i] + "_Specific_Surface_Area"; 
      aux_data(iv,aux_chem_variables[label]) = chem_state.mineral_specific_surface_area.data[i];
    }
    for (int i=0; i<surfSiteNames.size(); ++i) {
      const std::string label=surfSiteNames[i] + "_Surface_Site_Density"; 
      aux_data(iv,aux_chem_variables[label]) = chem_state.surface_site_density.data[i];
    }
    for (int i=0; i<ionExchangeNames.size(); ++i) {
      const std::string label=ionExchangeNames[i] + "_Ion_Exchange_Capacity"; 
      aux_data(iv,aux_chem_variables[label]) = chem_state.cation_exchange_capacity.data[i];
    }
    for (int i=0; i<num_free_ion_species; ++i) {
      const std::string label=BoxLib::Concatenate("Free_Ion_Species_",i,ndigits_freeion); 
      aux_data(iv,aux_chem_variables[label]) = aux_output.primary_free_ion_concentration.data[i]; 
    }
  }
}
