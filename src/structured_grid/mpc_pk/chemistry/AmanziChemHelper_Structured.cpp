
#include <AmanziChemHelper_Structured.H>

#include <cmath>
#include <Utility.H>
#include <ParallelDescriptor.H>

#include <chemistry_exception.hh>

static
int CopyStrArray(Array<std::string>& dest, const std::vector<std::string>& src) 
{
  int Nsize = src.size();
  dest.resize(Nsize);
  for (int i=0; i<Nsize; ++i) {
    dest[i] = src[i];
  }
  return Nsize;
}

static bool abort_on_chem_fail = true;

AmanziChemHelper_Structured::AmanziChemHelper_Structured(const std::vector<std::string>& _primarySpeciesNames,
                                                         const std::vector<std::string>& _sorbedPrimarySpeciesNames,
                                                         const std::vector<std::string>& _mineralNames,
                                                         const std::vector<std::string>& _surfaceComplexationSiteNames,
                                                         bool                            _hasCationExchangeCapacity,
                                                         const std::vector<std::string>& _isothermSpeciesNames,
                                                         const std::vector<std::string>& _freeIonSpeciesNames,
                                                         const std::string&              _thermo_database_filename,
                                                         const std::string&              _thermo_database_format,
                                                         const std::string&              _activity_model,
                                                         int                             _verbose)
  : using_sorption(false),
    using_isotherms(false),
    NionExchangeSites(0),
    hasCationExchangeCapacity(_hasCationExchangeCapacity),
    thermo_database_file(_thermo_database_filename),
    thermo_database_format(_thermo_database_format),
    activity_model(_activity_model),
    verbose(_verbose)
{
  if (hasCationExchangeCapacity) {
    using_sorption = true;
    NionExchangeSites = 1;
  }
  
  NprimarySpecies           = CopyStrArray(primarySpeciesNames, _primarySpeciesNames);
  NsorbedPrimarySpecies     = CopyStrArray(sorbedPrimarySpeciesNames, _sorbedPrimarySpeciesNames);
  Nminerals                 = CopyStrArray(mineralNames,_mineralNames);
  NsurfaceComplexationSites = CopyStrArray(surfaceComplexationSiteNames,_surfaceComplexationSiteNames);
  NfreeIonSpecies           = CopyStrArray(freeIonSpeciesNames,_freeIonSpeciesNames);
  NisothermSpecies          = CopyStrArray(isothermSpeciesNames,_isothermSpeciesNames);
  BL_ASSERT(NisothermSpecies == 0 || NisothermSpecies == NprimarySpecies);

  if (NsorbedPrimarySpecies > 0) {
    BL_ASSERT(NsorbedPrimarySpecies == NprimarySpecies);
    using_sorption = true;
  }

  if (NisothermSpecies > 0) {
    BL_ASSERT(NisothermSpecies == NprimarySpecies);
    using_isotherms = true;
    using_sorption = true;
  }

  // Create the list of variables that will be managed by Amanzi
  aux_chem_variables.clear();
  SmineralVolumeFraction = aux_chem_variables.size();
  for (int i=0; i<Nminerals; ++i) {
    aux_chem_variables[BuildPropertyParameterName(mineralNames[i],"Volume_Fraction")] = aux_chem_variables.size()-1;
  }
  SmineralSpecificSurfaceArea = aux_chem_variables.size();
  for (int i=0; i<Nminerals; ++i) {
    aux_chem_variables[BuildPropertyParameterName(mineralNames[i],"Specific_Surface_Area")] = aux_chem_variables.size()-1;
  }
  SsurfaceComplexationSiteDensity = aux_chem_variables.size();
  for (int i=0; i<NsurfaceComplexationSites; ++i) {
    aux_chem_variables[BuildPropertyParameterName(surfaceComplexationSiteNames[i],"Surface_Complexation","Site_Density")] = aux_chem_variables.size()-1;
  }
  SsurfaceComplexationFreeSiteConcentration = aux_chem_variables.size();
  for (int i=0; i<NsurfaceComplexationSites; ++i) {
    aux_chem_variables[BuildPropertyParameterName(surfaceComplexationSiteNames[i],"Surface_Complexation","Free_Site_Concentration")] = aux_chem_variables.size()-1;
  }
  ScationExchangeCapacity = aux_chem_variables.size();
  NionExchangeSites = (hasCationExchangeCapacity ? 1 : 0);

  int ndigIES = std::log(NionExchangeSites+1);
  SionExchangeSiteDensity = aux_chem_variables.size();
  for (int i=0; i<NionExchangeSites; ++i) {
    aux_chem_variables[BoxLib::Concatenate("Ion_Exchange_Site_Density_",i,ndigIES)] = aux_chem_variables.size()-1;
  }
  SionExchangeReferenceCationConcentration = aux_chem_variables.size();
  for (int i=0; i<NionExchangeSites; ++i) {
    aux_chem_variables[BoxLib::Concatenate("Ion_Exchange_Reference_Cation_Concentration_",i,ndigIES)] = aux_chem_variables.size()-1;
  }
  SfreeIonConcentration = aux_chem_variables.size();
  for (int i=0; i<NfreeIonSpecies; ++i) {
    aux_chem_variables[BuildPropertyParameterName(freeIonSpeciesNames[i],"Free_Ion","Guess")] = aux_chem_variables.size()-1;
  }
  SactivityCoefficient = aux_chem_variables.size();
  for (int i=0; i<NprimarySpecies; ++i) {
    aux_chem_variables[BuildPropertyParameterName(primarySpeciesNames[i] ,"Activity","Coefficient")] = aux_chem_variables.size()-1;
  }
  if (using_sorption) {
    Ssorbed = aux_chem_variables.size();
    for (int i=0; i<NsorbedPrimarySpecies; ++i) {
      aux_chem_variables[BuildPropertyParameterName(sorbedPrimarySpeciesNames[i],"Sorbed","Concentration")] = aux_chem_variables.size()-1;
    }
  }
  if (using_isotherms) {
    SisothermKd = aux_chem_variables.size();
    for (int i=0; i<isothermSpeciesNames.size(); ++i) {
      aux_chem_variables[BuildPropertyParameterName(isothermSpeciesNames[i],"Isotherm","Kd")] = aux_chem_variables.size()-1;
    }
    SisothermFreundlichN = aux_chem_variables.size();
    for (int i=0; i<isothermSpeciesNames.size(); ++i) {
      aux_chem_variables[BuildPropertyParameterName(isothermSpeciesNames[i],"Isotherm","Freundlich_n")] = aux_chem_variables.size()-1;
    }
    SisothermLangmuirB = aux_chem_variables.size();
    for (int i=0; i<isothermSpeciesNames.size(); ++i) {
      aux_chem_variables[BuildPropertyParameterName(isothermSpeciesNames[i],"Isotherm","Langmuir_b")] = aux_chem_variables.size()-1;
    }
  }

  if (ParallelDescriptor::IOProcessor() && verbose>1) {
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
  Real water_density = 1000;
  Real porosity = 1;
  Real saturation = 1;
  Real volume = 1;
  Real total_conc = 0;
  Real free_ion_conc = 1.e-9;
  Real mineral_vol_frac = 0;
  Real mineral_specific_surf_area = 0;
  Real total_sorbed = 0;
  Real ion_exchange_sites = 1.e-9;
  Real ion_exchange_ref_cation_conc = 1.e-9;
  Real primary_activity_coef = 1;
  Real isotherm_Kd = 1.e8;
  Real langmuir_b = 1.e-11;
  Real freundlich_n = 5.e-12;

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
    chemSolve.set(ithread, new Amanzi::AmanziChemistry::SimpleThermoDatabase());
	  
    parameters[ithread] = chemSolve[ithread].GetDefaultParameters();
    parameters[ithread].thermo_database_file = thermo_database_file;
    parameters[ithread].activity_model_name  = activity_model;
    parameters[ithread].porosity             = porosity; 
    parameters[ithread].saturation           = saturation;
    parameters[ithread].volume               = volume;
    parameters[ithread].water_density        = water_density;
	  
    components[ithread].total.resize(NprimarySpecies,total_conc);
    components[ithread].free_ion.resize(NfreeIonSpecies, free_ion_conc); 
    components[ithread].primary_activity_coeff.resize(NprimarySpecies, primary_activity_coef); 
    components[ithread].mineral_volume_fraction.resize(Nminerals, mineral_vol_frac);
    components[ithread].mineral_specific_surface_area.resize(Nminerals, mineral_specific_surf_area);
    if (using_sorption) { 
      components[ithread].total_sorbed.resize(NsorbedPrimarySpecies, total_sorbed);
    }
    if (NionExchangeSites>0) {
      components[ithread].ion_exchange_sites.resize(NionExchangeSites, ion_exchange_sites);
      components[ithread].ion_exchange_ref_cation_conc.resize(NionExchangeSites, ion_exchange_ref_cation_conc);
    }
    if (using_isotherms) {
      components[ithread].isotherm_kd.resize(NisothermSpecies,isotherm_Kd);
      components[ithread].isotherm_freundlich_n.resize(NisothermSpecies,freundlich_n);
      components[ithread].isotherm_langmuir_b.resize(NisothermSpecies,langmuir_b);
    }

    chemSolve[ithread].verbosity(Amanzi::AmanziChemistry::kTerse);	  
    chemSolve[ithread].Setup(components[ithread], parameters[ithread]);
    chemSolve[ithread].CopyBeakerToComponents(&(components[ithread]));

    if (ParallelDescriptor::IOProcessor() && ithread == 0) {
      chemSolve[ithread].Display();
      chemSolve[ithread].DisplayComponents(components[ithread]);
    }
  }  // for(threads)
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
    if (NionExchangeSites>0) {
      components[ithread].ion_exchange_sites.clear();
      components[ithread].ion_exchange_sites.clear();
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
                                              FArrayBox& auxiliary_data, bool initialize_auxiliary_data,
                                              const Box& box, const std::string& condition_name, Real time)
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
                                     FArrayBox&       auxData, Real water_density, Real temperature,
                                     const Box& box, Real dt, bool initialize)
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
      Amanzi::AmanziChemistry::Beaker::BeakerComponents& TheComponent = components[threadid];
      Amanzi::AmanziChemistry::Beaker::BeakerParameters& TheParameter = parameters[threadid];
      
      bool is_neg = false;
      for (int i = 0; i < NprimarySpecies; ++i) {
        TheComponent.total[i] = primary_species_mobile(iv,sPrimMob+i);
        if (std::abs(TheComponent.total[i]) < 1.e-16) TheComponent.total[i] = 0;
        if (TheComponent.total[i] < 0) is_neg = true;
      }
      if (is_neg) {
        break;
      }

      for (int i = 0; i < NfreeIonSpecies; ++i) {
        TheComponent.free_ion[i] = auxData(iv,SfreeIonConcentration+i);
      }

      for (int i = 0; i < NprimarySpecies; ++i) {
        TheComponent.primary_activity_coeff[i] = auxData(iv,SactivityCoefficient+i);
      }

      for (int i = 0; i < Nminerals; ++i) {
        TheComponent.mineral_volume_fraction[i] = auxData(iv,SmineralVolumeFraction+i);
        TheComponent.mineral_specific_surface_area[i] = auxData(iv,SmineralSpecificSurfaceArea+i);
      }

      for (int i = 0; i < NionExchangeSites; ++i) {
        TheComponent.ion_exchange_sites[i] = auxData(iv,SionExchangeSiteDensity+i);
        TheComponent.ion_exchange_ref_cation_conc[i] = auxData(iv,SionExchangeReferenceCationConcentration+i);
      }

      for (int i = 0; i < NsurfaceComplexationSites; ++i) {
        TheComponent.surface_site_density[i] = auxData(iv,SsurfaceComplexationSiteDensity+i);
      }

      for (int i = 0; i < NsurfaceComplexationSites; ++i) {
        TheComponent.surface_complex_free_site_conc[i] = auxData(iv,SsurfaceComplexationFreeSiteConcentration+i);
      }

      if (using_sorption) {
        for (int i = 0; i < NsorbedPrimarySpecies; ++i) {
          TheComponent.total_sorbed[i] = auxData(iv,Ssorbed+i);
        }
      }

      if (using_isotherms) {
        for (int i = 0; i < NisothermSpecies; ++i) {
          TheComponent.isotherm_kd[i] = auxData(iv,SisothermKd+i);
          TheComponent.isotherm_freundlich_n[i] = auxData(iv,SisothermFreundlichN+i);
          TheComponent.isotherm_langmuir_b[i] = auxData(iv,SisothermLangmuirB+i);
        }
      }

      TheParameter.porosity   = porosity(iv,sPhi);
      TheParameter.saturation = std::min(1., std::max(0., aqueous_saturation(iv,sSat)));
      TheParameter.volume     = volume(iv,sVol);
      TheParameter.water_density = water_density;
	
      chem_ok = true;

      if (initialize)
      {
        TheChemSolve.Speciate(&TheComponent,TheParameter);
      }
      else {
        Amanzi::AmanziChemistry::Beaker::SolverStatus stat;
        try
        {  
          TheChemSolve.ReactionStep(&TheComponent,TheParameter,dt);

          stat = TheChemSolve.status();
          fcnCnt(iv,sFunc) = stat.num_rhs_evaluations;        
        }
        catch (const Amanzi::AmanziChemistry::ChemistryException& geochem_error)
        {
          if (verbose>-1) {
            std::cout << "CHEMISTRY FAILED on level at " << iv << " : ";
            TheComponent.Display("components: ");
            if (abort_on_chem_fail) {
              BoxLib::Abort(geochem_error.what());
            }
          } else {
            chem_ok = false;
          }
        }
      }
      
      // If successful update the state variables.
      if (chem_ok) {

        for (int i = 0; i < NprimarySpecies; ++i) {
          primary_species_mobile(iv,sPrimMob+i) = TheComponent.total[i];
        }
        for (int i = 0; i < NfreeIonSpecies; ++i) {
          auxData(iv,SfreeIonConcentration+i) = TheComponent.free_ion[i];
        }

        for (int i = 0; i < NprimarySpecies; ++i) {
          auxData(iv,SactivityCoefficient+i) = TheComponent.primary_activity_coeff[i];
        }

        for (int i = 0; i < Nminerals; ++i) {
          auxData(iv,SmineralVolumeFraction+i) = TheComponent.mineral_volume_fraction[i];
          auxData(iv,SmineralSpecificSurfaceArea+i) = TheComponent.mineral_specific_surface_area[i];
        }

        for (int i = 0; i < NionExchangeSites; ++i) {
          auxData(iv,SionExchangeSiteDensity+i) = TheComponent.ion_exchange_sites[i];
          auxData(iv,SionExchangeReferenceCationConcentration+i) = TheComponent.ion_exchange_ref_cation_conc[i];
        }
        
        for (int i = 0; i < NsurfaceComplexationSites; ++i) {
          auxData(iv,SsurfaceComplexationSiteDensity+i) = TheComponent.surface_site_density[i];
        }

        if (using_sorption) {
          for (int i = 0; i < NsorbedPrimarySpecies; ++i) {
            auxData(iv,Ssorbed+i) = TheComponent.total_sorbed[i];
          }
        }

        if (using_isotherms) {
          for (int i = 0; i < NisothermSpecies; ++i) {
            auxData(iv,SisothermKd+i) = TheComponent.isotherm_kd[i];
            auxData(iv,SisothermFreundlichN+i) = TheComponent.isotherm_freundlich_n[i];
            auxData(iv,SisothermLangmuirB+i) = TheComponent.isotherm_langmuir_b[i];
          }
        }
      } // chem_ok

    }
  }
}

