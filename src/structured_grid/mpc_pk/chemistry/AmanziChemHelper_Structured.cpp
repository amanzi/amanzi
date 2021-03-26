
#include <AmanziChemHelper_Structured.H>

#ifdef _OPENMP
#include "omp.h"
#endif

#include <cmath>
#include <Utility.H>
#include <ParallelDescriptor.H>

#include <chemistry_exception.hh>

static bool abort_on_chem_fail = true;
//#define DEBUG_NO_CHEM 
#undef DEBUG_NO_CHEM 

AmanziChemHelper_Structured::AmanziChemHelper_Structured(const std::vector<std::string>& _primarySpeciesNames,
                                                         const std::vector<std::string>& _sorbedPrimarySpeciesNames,
                                                         const std::vector<std::string>& _mineralNames,
                                                         const std::vector<std::string>& _surfaceComplexationSiteNames,
                                                         bool                            _hasCationExchangeCapacity,
                                                         const std::vector<std::string>& _isothermSpeciesNames,
                                                         const std::vector<std::string>& _freeIonSpeciesNames,
                                                         const std::string&              _thermo_database_filename,
                                                         const std::string&              _thermo_database_format,
                                                         const std::string&              _activity_model)
  : thermo_database_file(_thermo_database_filename),
    thermo_database_format(_thermo_database_format),
    activity_model(_activity_model)
{
  using_sorption = false;
  using_isotherms = false;

  primarySpeciesNames = _primarySpeciesNames;
  mineralNames = _mineralNames;
  surfSiteNames = _surfaceComplexationSiteNames;
  Nmobile = primarySpeciesNames.size();
  Nimmobile = _sorbedPrimarySpeciesNames.size();
  Nminerals = mineralNames.size();
  Nisotherms = _isothermSpeciesNames.size();
  NionExchange = (_hasCationExchangeCapacity ? 1 : 0);
  NsorptionSites = _surfaceComplexationSiteNames.size();
  NfreeIonSpecies = _freeIonSpeciesNames.size();

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

  if (ParallelDescriptor::IOProcessor()) {
    std::cout << "AmanziChemHelper_Structured: Iniitialized" << std::endl;
    std::cout << "   Auxiliary data structure will contain the following components:" << std::endl;

    Array<std::string> tmp(aux_chem_variables.size());
    for (std::map<std::string,int>::const_iterator it = aux_chem_variables.begin(); it!=aux_chem_variables.end(); ++it) {
      BL_ASSERT(it->second<tmp.size());
      tmp[it->second] = it->first;
    }
    for (int i=0; i<aux_chem_variables.size(); ++i) {
      std::cout << "      " << i << ": " << tmp[i] << std::endl; 
    }
  }

  // Dummy parameters
  Real total_conc = 1.e-20;
  Real free_ion_conc = 1.e-9;
  Real mineral_vol_frac = 0;
  Real mineral_specific_surf_area = 0;
  Real total_sorbed = 0;
  Real ion_exchange_sites = 1;
  Real ion_exchange_ref_cation_conc = 1;
  Real primary_activity_coef = 1;
  Real isotherm_Kd = -1;
  Real langmuir_b = 1;
  Real freundlich_n = 1;
  Real surface_site_density = 1;

  int ntracers = primarySpeciesNames.size();
  int nminerals = mineralNames.size();
  //
  // In order to thread the AMANZI chemistry, we had to give each thread 
  // its own chemSolve and components object.
  //
  int tnum = 1;
#ifdef _OPENMP
  tnum = omp_get_max_threads();
#endif
  chemSolve.resize(tnum, PArrayManage);
  components.resize(tnum);
  parameters.resize(tnum);
      
  for (int ithread = 0; ithread < tnum; ithread++) {
    Teuchos::ParameterList vo_list;
    vo_ = Teuchos::rcp(new Amanzi::VerboseObject("Chemistry PK", vo_list));
    chemSolve.set(ithread, new Amanzi::AmanziChemistry::SimpleThermoDatabase(vo_));
	  
    parameters[ithread] = chemSolve[ithread].GetDefaultParameters();
    parameters[ithread].thermo_database_file = thermo_database_file;
    parameters[ithread].activity_model_name  = activity_model;
 
    components[ithread].total.resize(Nmobile,total_conc);
    components[ithread].free_ion.resize(NfreeIonSpecies, free_ion_conc); 
    components[ithread].primary_activity_coeff.resize(Nmobile, primary_activity_coef); 
    components[ithread].mineral_volume_fraction.resize(Nminerals, mineral_vol_frac);
    components[ithread].mineral_specific_surface_area.resize(Nminerals, mineral_specific_surf_area);
    if (using_sorption) { 
      components[ithread].total_sorbed.resize(Nimmobile, total_sorbed);
    }
    if (NionExchange>0) {
      components[ithread].ion_exchange_sites.resize(NionExchange, ion_exchange_sites);
      components[ithread].ion_exchange_ref_cation_conc.resize(NionExchange, ion_exchange_ref_cation_conc);
    }
    if (using_isotherms) {
      components[ithread].isotherm_kd.resize(Nisotherms,isotherm_Kd);
      components[ithread].isotherm_freundlich_n.resize(Nisotherms,freundlich_n);
      components[ithread].isotherm_langmuir_b.resize(Nisotherms,langmuir_b);
    }
    if (NsorptionSites > 0) {
      components[ithread].surface_site_density.resize(NsorptionSites,surface_site_density);
    }
    
    chemSolve[ithread].Setup(components[ithread], parameters[ithread]);
  }
}

AmanziChemHelper_Structured::~AmanziChemHelper_Structured()
{
  int tnum = 1;
#ifdef _OPENMP
  tnum = omp_get_max_threads();
#endif
      
  for (int ithread = 0; ithread < tnum; ithread++) {
	  
    components[ithread].total.clear();
    components[ithread].free_ion.clear();
    components[ithread].primary_activity_coeff.clear();
    components[ithread].mineral_volume_fraction.clear();
    if (using_sorption) { 
      components[ithread].total_sorbed.clear();
    }
    if (NionExchange>0) {
      components[ithread].ion_exchange_sites.clear();
      components[ithread].ion_exchange_ref_cation_conc.clear();
    }
    if (using_isotherms) {
      components[ithread].isotherm_kd.clear();
      components[ithread].isotherm_freundlich_n.clear();
      components[ithread].isotherm_langmuir_b.clear();
    }
    chemSolve.clear(ithread);
  }  // for(threads)

  chemSolve.clear();
  components.clear();
  parameters.clear();
}

void
AmanziChemHelper_Structured::EnforceCondition(FArrayBox& primary_species_mobile,   int sPrimMob,
                                              FArrayBox& auxiliary_data, Real water_density, Real temperature,
                                              const Box& box, const std::string& condition_name, Real time,
					      int chem_verbose)
{
  BoxLib::Abort("Geochemical conditions/constraints not currently support in Amanzi native chemistry engine");
}

void
AmanziChemHelper_Structured::Advance(const FArrayBox& aqueous_saturation,       int sSat,
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

    for (IntVect iv=thread_box.smallEnd(), End=thread_box.bigEnd(); iv<=End && chem_ok; thread_box.next(iv)) {
      Amanzi::AmanziChemistry::SimpleThermoDatabase&     TheChemSolve = chemSolve[threadid];
      Amanzi::AmanziChemistry::Beaker::BeakerState&      TheComponent = components[threadid];
      Amanzi::AmanziChemistry::Beaker::BeakerParameters& TheParameter = parameters[threadid];
      
      TheParameter.volume     = volume(iv,sVol);
      TheParameter.saturation = std::min(1., std::max(0., aqueous_saturation(iv,sSat)));

      TheParameter.water_density = water_density;
      TheParameter.porosity   = porosity(iv,sPhi);

      bool is_neg = false;
      for (int i = 0; i < Nmobile; ++i) {
        TheComponent.total[i] = primary_species_mobile(iv,sPrimMob+i);
        if (TheComponent.total[i] < 0) is_neg = true;
        TheComponent.total[i] = std::max(0.,TheComponent.total[i]);
      }
      if (is_neg) {
        break;
      }

      if (using_sorption) {
        for (int i=0; i<primarySpeciesNames.size(); ++i) {
          const std::string label=primarySpeciesNames[i] + "_Sorbed_Concentration"; 
          TheComponent.total_sorbed[i] = aux_data(iv,aux_chem_variables[label]);
        }
      }

      if (Nminerals > 0) {
        for (int i=0; i<Nminerals; ++i) {
          const std::string label=mineralNames[i] + "_Volume_Fraction"; 
          TheComponent.mineral_volume_fraction[i] = aux_data(iv,aux_chem_variables[label]);
        }
        for (int i=0; i<mineralNames.size(); ++i) {
          const std::string label=mineralNames[i] + "_Specific_Surface_Area"; 
          TheComponent.mineral_specific_surface_area[i] = aux_data(iv,aux_chem_variables[label]);
        }
      }

      if (NionExchange > 0) {
        int ndigIES = std::log(NionExchange+1);
        for (int i=0; i<NionExchange; ++i) {
          const std::string label = BoxLib::Concatenate("Ion_Exchange_Site_Density_",i,ndigIES);
          TheComponent.ion_exchange_sites[i] = aux_data(iv,aux_chem_variables[label]);
        }
        
        for (int i=0; i<NionExchange; ++i) {
          const std::string label = BoxLib::Concatenate("Ion_Exchange_Reference_Cation_Concentration_",i,ndigIES);
          TheComponent.ion_exchange_ref_cation_conc[i] = aux_data(iv,aux_chem_variables[label]);
        }
      }
      
      if (NsorptionSites > 0) {
        for (int i=0; i<surfSiteNames.size(); ++i) {
          const std::string label=surfSiteNames[i] + "_Surface_Site_Density"; 
          TheComponent.surface_site_density[i] = aux_data(iv,aux_chem_variables[label]);
        }
      }

      for (int i = 0; i < Nmobile; ++i) {
        const std::string label=primarySpeciesNames[i] + "_Activity_Coefficient"; 
        TheComponent.primary_activity_coeff[i] = aux_data(iv,aux_chem_variables[label]);
      }

      if (NfreeIonSpecies > 0) {
        for (int i=0; i<primarySpeciesNames.size(); ++i) {
          const std::string label=primarySpeciesNames[i] + "_Free_Ion_Guess"; 
          //TheComponent.free_ion[i] = aux_data(iv,aux_chem_variables[label]);
          TheComponent.free_ion[i] = 1.e-20;
        }
      }
      
      if (using_isotherms) {
        for (int i=0; i<Nisotherms; ++i) {
          const std::string label=primarySpeciesNames[i] + "_Isotherm_Kd"; 
          TheComponent.isotherm_kd[i] = aux_data(iv,aux_chem_variables[label]);
        }
        for (int i=0; i<Nisotherms; ++i) {
          const std::string label=primarySpeciesNames[i] + "_Isotherm_Freundlich_n"; 
          TheComponent.isotherm_freundlich_n[i] = aux_data(iv,aux_chem_variables[label]);
        }
        for (int i=0; i<Nisotherms; ++i) {
          const std::string label=primarySpeciesNames[i] + "_Isotherm_Langmuir_b"; 
          TheComponent.isotherm_langmuir_b[i] = aux_data(iv,aux_chem_variables[label]);
        }
      }

      chem_ok = true;

      if (ParallelDescriptor::IOProcessor() && chem_verbose>0) {
	std::cout << "iv: " << iv << std::endl;
	std::cout << "dt: " << dt << std::endl;
	DumpChemStructures(std::cout,TheChemSolve,TheComponent,TheParameter);
      }

      Amanzi::AmanziChemistry::Beaker::SolverStatus stat;
      try
      { 
#ifndef DEBUG_NO_CHEM 
        TheChemSolve.ReactionStep(&TheComponent,TheParameter,dt);
#endif
        stat = TheChemSolve.status();
        fcnCnt(iv,sFunc) = stat.num_rhs_evaluations;        
      }
      catch (const Amanzi::AmanziChemistry::ChemistryException& geochem_error)
      {
        if (chem_verbose>=0) {
          std::cout << "CHEMISTRY FAILED on level at " << iv << " : ";
          TheComponent.Display("components: ", vo_);
          if (abort_on_chem_fail) {
            BoxLib::Abort(geochem_error.what());
          }
        } else {
          chem_ok = false;
        }
      }
      
      // If successful update the state variables.
      if (chem_ok) {

        for (int i = 0; i < Nmobile; ++i) {
          primary_species_mobile(iv,sPrimMob+i) = TheComponent.total[i];
        }

        if (using_sorption) {
          for (int i=0; i<Nimmobile; ++i) {
            const std::string label=primarySpeciesNames[i] + "_Sorbed_Concentration"; 
            aux_data(iv,aux_chem_variables[label]) = TheComponent.total_sorbed[i];
          }
        }
#if 0
        for (int i=0; i<Nmobile; ++i) {
          const std::string label=primarySpeciesNames[i] + "_Activity_Coefficient"; 
          aux_data(iv,aux_chem_variables[label]) = TheComponent.primary_activity_coeff[i];
        }

        if (NfreeIonSpecies > 0) {
          for (int i=0; i<NfreeIonSpecies; ++i) {
            const std::string label=primarySpeciesNames[i] + "_Free_Ion_Guess"; 
            aux_data(iv,aux_chem_variables[label]) = TheComponent.free_ion[i];
          }
        }
#endif
        if (Nminerals > 0) {
          for (int i=0; i<Nminerals; ++i) {
            const std::string label=mineralNames[i] + "_Volume_Fraction"; 
            aux_data(iv,aux_chem_variables[label]) = TheComponent.mineral_volume_fraction[i];
          }
          for (int i=0; i<mineralNames.size(); ++i) {
            const std::string label=mineralNames[i] + "_Specific_Surface_Area"; 
            aux_data(iv,aux_chem_variables[label]) = TheComponent.mineral_specific_surface_area[i];
          }
        }

        if (NsorptionSites > 0) {
          for (int i=0; i<surfSiteNames.size(); ++i) {
            const std::string label=surfSiteNames[i] + "_Surface_Site_Density"; 
            aux_data(iv,aux_chem_variables[label]) = TheComponent.surface_site_density[i];
          }
        }
        
        if (NionExchange > 0) {
          int ndigIES = std::log(NionExchange+1);
          for (int i=0; i<NionExchange; ++i) {
            const std::string label = BoxLib::Concatenate("Ion_Exchange_Site_Density_",i,ndigIES);
            aux_data(iv,aux_chem_variables[label]) = TheComponent.ion_exchange_sites[i];
          }
          for (int i=0; i<NionExchange; ++i) {
            const std::string label = BoxLib::Concatenate("Ion_Exchange_Reference_Cation_Concentration_",i,ndigIES);
            aux_data(iv,aux_chem_variables[label]) = TheComponent.ion_exchange_ref_cation_conc[i];
          }
        }

        if (using_isotherms) {
          for (int i=0; i<Nisotherms; ++i) {
            const std::string label=primarySpeciesNames[i] + "_Isotherm_Kd"; 
            aux_data(iv,aux_chem_variables[label]) = TheComponent.isotherm_kd[i];
          }
          for (int i=0; i<Nisotherms; ++i) {
            const std::string label=primarySpeciesNames[i] + "_Isotherm_Freundlich_n"; 
            aux_data(iv,aux_chem_variables[label]) = TheComponent.isotherm_freundlich_n[i];
          }
          for (int i=0; i<Nisotherms; ++i) {
            const std::string label=primarySpeciesNames[i] + "_Isotherm_Langmuir_b"; 
            aux_data(iv,aux_chem_variables[label]) = TheComponent.isotherm_langmuir_b[i];
          }
        }
      } // chem_ok

    }
  }
}

void
AmanziChemHelper_Structured::Initialize(const FArrayBox& aqueous_saturation,       int sSat,
                                        const FArrayBox& aqueous_pressure,         int sPress,
                                        const FArrayBox& porosity,                 int sPhi,
                                        const FArrayBox& volume,                   int sVol,
                                        FArrayBox&       primary_species_mobile,   int sPrimMob,
                                        FArrayBox&       fcnCnt,                   int sFunc,
                                        FArrayBox&       aux_data, Real water_density, Real temperature,
                                        const Box& box)
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

    for (IntVect iv=thread_box.smallEnd(), End=thread_box.bigEnd(); iv<=End && chem_ok; thread_box.next(iv)) {
      Amanzi::AmanziChemistry::SimpleThermoDatabase&     TheChemSolve = chemSolve[threadid];
      Amanzi::AmanziChemistry::Beaker::BeakerState&      TheComponent = components[threadid];
      Amanzi::AmanziChemistry::Beaker::BeakerParameters& TheParameter = parameters[threadid];
      
      TheParameter.volume     = volume(iv,sVol);
      TheParameter.saturation = std::min(1., std::max(0., aqueous_saturation(iv,sSat)));

      TheParameter.water_density = water_density;
      TheParameter.porosity   = porosity(iv,sPhi);

      bool is_neg = false;
      for (int i = 0; i < Nmobile; ++i) {
        TheComponent.total[i] = primary_species_mobile(iv,sPrimMob+i);
        if (TheComponent.total[i] < 0) is_neg = true;
        TheComponent.total[i] = std::max(0.,TheComponent.total[i]);
      }
      if (is_neg) {
        break;
      }

      if (NfreeIonSpecies > 0) {
        for (int i=0; i<primarySpeciesNames.size(); ++i) {
          const std::string label=primarySpeciesNames[i] + "_Free_Ion_Guess"; 
          //TheComponent.free_ion[i] = aux_data(iv,aux_chem_variables[label]);
          TheComponent.free_ion[i] = 1.e-20;
        }
      }
      
      for (int i = 0; i < Nmobile; ++i) {
        const std::string label=primarySpeciesNames[i] + "_Activity_Coefficient"; 
        TheComponent.primary_activity_coeff[i] = aux_data(iv,aux_chem_variables[label]);
      }

      if (Nminerals > 0) {
        for (int i=0; i<Nminerals; ++i) {
          const std::string label=mineralNames[i] + "_Volume_Fraction"; 
          TheComponent.mineral_volume_fraction[i] = aux_data(iv,aux_chem_variables[label]);
        }
        for (int i=0; i<mineralNames.size(); ++i) {
          const std::string label=mineralNames[i] + "_Specific_Surface_Area"; 
          TheComponent.mineral_specific_surface_area[i] = aux_data(iv,aux_chem_variables[label]);
        }
      }

      if (using_sorption) {
        for (int i=0; i<primarySpeciesNames.size(); ++i) {
          const std::string label=primarySpeciesNames[i] + "_Sorbed_Concentration"; 
          TheComponent.total_sorbed[i] = aux_data(iv,aux_chem_variables[label]);
        }
      }

      if (NionExchange > 0) {
        int ndigIES = std::log(NionExchange+1);
        for (int i=0; i<NionExchange; ++i) {
          const std::string label = BoxLib::Concatenate("Ion_Exchange_Site_Density_",i,ndigIES);
          TheComponent.ion_exchange_sites[i] = aux_data(iv,aux_chem_variables[label]);
        }
        
        for (int i=0; i<NionExchange; ++i) {
          const std::string label = BoxLib::Concatenate("Ion_Exchange_Reference_Cation_Concentration_",i,ndigIES);
          TheComponent.ion_exchange_ref_cation_conc[i] = aux_data(iv,aux_chem_variables[label]);
        }
      }
      
      if (NsorptionSites > 0) {
        for (int i=0; i<surfSiteNames.size(); ++i) {
          const std::string label=surfSiteNames[i] + "_Surface_Site_Density"; 
          TheComponent.surface_site_density[i] = aux_data(iv,aux_chem_variables[label]);
        }
      }

      if (using_isotherms) {
        for (int i=0; i<Nisotherms; ++i) {
          const std::string label=primarySpeciesNames[i] + "_Isotherm_Kd"; 
          TheComponent.isotherm_kd[i] = aux_data(iv,aux_chem_variables[label]);
        }
        for (int i=0; i<Nisotherms; ++i) {
          const std::string label=primarySpeciesNames[i] + "_Isotherm_Freundlich_n"; 
          TheComponent.isotherm_freundlich_n[i] = aux_data(iv,aux_chem_variables[label]);
        }
        for (int i=0; i<Nisotherms; ++i) {
          const std::string label=primarySpeciesNames[i] + "_Isotherm_Langmuir_b"; 
          TheComponent.isotherm_langmuir_b[i] = aux_data(iv,aux_chem_variables[label]);
        }
      }

      chem_ok = true;

      Amanzi::AmanziChemistry::Beaker::SolverStatus stat;
      try
      {  
#ifndef DEBUG_NO_CHEM 
        TheChemSolve.Speciate(&TheComponent,TheParameter);
#endif
        stat = TheChemSolve.status();
        fcnCnt(iv,sFunc) = stat.num_rhs_evaluations;        
      }
      catch (const Amanzi::AmanziChemistry::ChemistryException& geochem_error)
      {
	std::cout << "CHEMISTRY SPECIATION FAILED on level at " << iv << " : ";
	TheComponent.Display("components: ", vo_);
	if (abort_on_chem_fail) {
	  BoxLib::Abort(geochem_error.what());
	}
	chem_ok = false;
      }

      // If successful update the state variables.
      if (chem_ok) {

        for (int i = 0; i < Nmobile; ++i) {
          primary_species_mobile(iv,sPrimMob+i) = TheComponent.total[i];
        }

        if (NfreeIonSpecies > 0) {
          for (int i=0; i<NfreeIonSpecies; ++i) {
            const std::string label=primarySpeciesNames[i] + "_Free_Ion_Guess"; 
            aux_data(iv,aux_chem_variables[label]) = TheComponent.free_ion[i];
          }
        }

        for (int i=0; i<Nmobile; ++i) {
          const std::string label=primarySpeciesNames[i] + "_Activity_Coefficient"; 
          aux_data(iv,aux_chem_variables[label]) = TheComponent.primary_activity_coeff[i];
        }

        if (Nminerals > 0) {
          for (int i=0; i<Nminerals; ++i) {
            const std::string label=mineralNames[i] + "_Volume_Fraction"; 
            aux_data(iv,aux_chem_variables[label]) = TheComponent.mineral_volume_fraction[i];
          }
          for (int i=0; i<mineralNames.size(); ++i) {
            const std::string label=mineralNames[i] + "_Specific_Surface_Area"; 
            aux_data(iv,aux_chem_variables[label]) = TheComponent.mineral_specific_surface_area[i];
          }
        }

        if (using_sorption) {
          for (int i=0; i<Nimmobile; ++i) {
            const std::string label=primarySpeciesNames[i] + "_Sorbed_Concentration"; 
            aux_data(iv,aux_chem_variables[label]) = TheComponent.total_sorbed[i];
          }
        }

        if (NsorptionSites > 0) {
          for (int i=0; i<surfSiteNames.size(); ++i) {
            const std::string label=surfSiteNames[i] + "_Surface_Site_Density"; 
            aux_data(iv,aux_chem_variables[label]) = TheComponent.surface_site_density[i];
          }
        }
        
        if (NionExchange > 0) {
          int ndigIES = std::log(NionExchange+1);
          for (int i=0; i<NionExchange; ++i) {
            const std::string label = BoxLib::Concatenate("Ion_Exchange_Site_Density_",i,ndigIES);
            aux_data(iv,aux_chem_variables[label]) = TheComponent.ion_exchange_sites[i];
          }
          for (int i=0; i<NionExchange; ++i) {
            const std::string label = BoxLib::Concatenate("Ion_Exchange_Reference_Cation_Concentration_",i,ndigIES);
            aux_data(iv,aux_chem_variables[label]) = TheComponent.ion_exchange_ref_cation_conc[i];
          }
        }

        if (using_isotherms) {
          for (int i=0; i<Nisotherms; ++i) {
            const std::string label=primarySpeciesNames[i] + "_Isotherm_Kd"; 
            aux_data(iv,aux_chem_variables[label]) = TheComponent.isotherm_kd[i];
          }
          for (int i=0; i<Nisotherms; ++i) {
            const std::string label=primarySpeciesNames[i] + "_Isotherm_Freundlich_n"; 
            aux_data(iv,aux_chem_variables[label]) = TheComponent.isotherm_freundlich_n[i];
          }
          for (int i=0; i<Nisotherms; ++i) {
            const std::string label=primarySpeciesNames[i] + "_Isotherm_Langmuir_b"; 
            aux_data(iv,aux_chem_variables[label]) = TheComponent.isotherm_langmuir_b[i];
          }
        }
      } // chem_ok
    }
  }
}

void
AmanziChemHelper_Structured::DumpChemStructures(std::ostream&                                      os,
						Amanzi::AmanziChemistry::SimpleThermoDatabase&     TheChemSolve,
						Amanzi::AmanziChemistry::Beaker::BeakerState&      TheComponent,
						Amanzi::AmanziChemistry::Beaker::BeakerParameters& TheParameter)
{
  os << "Volume " << TheParameter.volume << std::endl;
  os << "Saturation " << TheParameter.saturation << std::endl;
  os << "Water density " << TheParameter.water_density << std::endl;
  os << "Porosity " << TheParameter.porosity << std::endl;

  for (int i=0; i<Nmobile; ++i) {
    const std::string label=primarySpeciesNames[i] + "_Primary_Species_Mobile"; 
    os << label << "  " << TheComponent.total[i] << std::endl;
  }

  if (using_sorption) {
    for (int i=0; i<Nmobile; ++i) {
      const std::string label=primarySpeciesNames[i] + "_Sorbed_Concentration"; 
      os << label << "  " << TheComponent.total_sorbed[i] << std::endl;
    }
  }

  if (Nminerals > 0) {
    for (int i=0; i<Nminerals; ++i) {
      const std::string label=mineralNames[i] + "_Volume_Fraction"; 
      os << label << "  " << TheComponent.mineral_volume_fraction[i] << std::endl;
    }
    for (int i=0; i<mineralNames.size(); ++i) {
      const std::string label=mineralNames[i] + "_Specific_Surface_Area"; 
      os << label << "  " << TheComponent.mineral_specific_surface_area[i] << std::endl;
    }
  }

  if (NionExchange > 0) {
    int ndigIES = std::log(NionExchange+1);
    for (int i=0; i<NionExchange; ++i) {
      const std::string label = BoxLib::Concatenate("Ion_Exchange_Site_Density_",i,ndigIES);
      os << label << "  " << TheComponent.ion_exchange_sites[i] << std::endl;
    }

    for (int i=0; i<NionExchange; ++i) {
      const std::string label = BoxLib::Concatenate("Ion_Exchange_Reference_Cation_Concentration_",i,ndigIES);
      os << label << "  " << TheComponent.ion_exchange_ref_cation_conc[i] << std::endl;
    }
  }
      
  if (NsorptionSites > 0) {
    for (int i=0; i<surfSiteNames.size(); ++i) {
      const std::string label=surfSiteNames[i] + "_Surface_Site_Density"; 
      os << label << "  " << TheComponent.surface_site_density[i] << std::endl;
    }
  }

  for (int i = 0; i < Nmobile; ++i) {
    const std::string label=primarySpeciesNames[i] + "_Activity_Coefficient"; 
    os << label << "  " << TheComponent.primary_activity_coeff[i] << std::endl;
  }

  if (NfreeIonSpecies > 0) {
    for (int i=0; i<primarySpeciesNames.size(); ++i) {
      const std::string label=primarySpeciesNames[i] + "_Free_Ion_Guess"; 
      os << label << "  " << TheComponent.free_ion[i] << std::endl;
    }
  }

  if (using_isotherms) {
    for (int i=0; i<Nisotherms; ++i) {
      const std::string label=primarySpeciesNames[i] + "_Isotherm_Kd"; 
      os << label << "  " << TheComponent.isotherm_kd[i] << std::endl;
    }
    for (int i=0; i<Nisotherms; ++i) {
      const std::string label=primarySpeciesNames[i] + "_Isotherm_Freundlich_n"; 
      os << label << "  " << TheComponent.isotherm_freundlich_n[i] << std::endl;
    }
    for (int i=0; i<Nisotherms; ++i) {
      const std::string label=primarySpeciesNames[i] + "_Isotherm_Langmuir_b"; 
      os << label << "  " << TheComponent.isotherm_langmuir_b[i] << std::endl;
    }
  }
}
