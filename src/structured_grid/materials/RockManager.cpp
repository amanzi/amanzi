
#include <RockManager.H>

#include <ParmParse.H>

static std::string CapillaryPressureName    = "capillary_pressure";
static std::string PorosityName             = "porosity";
static std::string PermeabilityName         = "permeability";
static std::string RelativePermeabilityName = "relative_permeability";

static std::string CP_model_None = "None";
static std::string CP_model_vG = "VanGenuchten";
static std::string CP_model_BC = "BrooksCorey";

static std::string Kr_model_None = "None";
static std::string Kr_model_Mualem = "Mualem";
static std::string Kr_model_Burdine = "Burdine";

static std::string Kr_model_vG_Mualem  = CP_model_vG + "_" + Kr_model_Mualem;
static std::string Kr_model_vG_Burdine = CP_model_vG + "_" + Kr_model_Burdine;
static std::string Kr_model_BC_Mualem  = CP_model_BC + "_" + Kr_model_Mualem;
static std::string Kr_model_BC_Burdine = CP_model_BC + "_" + Kr_model_Burdine;

static int MAX_CPL_PARAMS = 6; // Must be set to accommodate the model with the most parameters
static int CPL_MODEL_ID = 0;

static int VG_M           = 1;
static int VG_ALPHA       = 2;
static int VG_SR          = 3;
static int VG_ELL         = 4;
static int VG_KR_MODEL_ID = 5;

static int BC_LAMBDA      = 1;
static int BC_ALPHA       = 2;
static int BC_SR          = 3;
static int BC_ELL         = 4;
static int BC_KR_MODEL_ID = 5;

// Interpolators
static int NUM_INIT_INTERP_EVAL_PTS = 1001;
static Real krel_smoothing_interval = 1.e-3;
static Real pc_at_Sr = 1.e11;

RockManager::RockManager(const RegionManager*   _region_manager,
                         const Array<Geometry>& geomArray,
                         const Array<IntVect>&  refRatio)
  : region_manager(_region_manager)
{
  is_initialized = false;
  Initialize(geomArray,refRatio);
  BuildInterpolators();
  is_initialized = true;
}

void
RockManager::BuildInterpolators()
{
  CP_s_interps.resize(rock.size(),PArrayManage);

  int nComp = materialFiller->NComp(CapillaryPressureName);
  static IntVect iv(D_DECL(0,0,0));
  static Box bx(iv,iv);
  FArrayBox pc_params(bx,nComp);
  int level=0; //not really used
  int dComp=0;
  Real time = 0;

  for (int n=0; n<rock.size(); ++n) {

    Array<Real> s(NUM_INIT_INTERP_EVAL_PTS);
    int Npts = s.size();
    Array<Real> pc(Npts);
    Array<int> mat(Npts,n);

    pc[0] = pc_at_Sr;
    InverseCapillaryPressure(pc.dataPtr(),&n,time,&(s[0]),1);
    Real ds = 1 - s[0];
    for (int i=1; i<s.size(); ++i) {
      s[i] = s[0] + ds*Real(i)/(NUM_INIT_INTERP_EVAL_PTS - 1);
    }
    CapillaryPressure(s.dataPtr(),mat.dataPtr(),time,pc.dataPtr(),Npts);

    //for (int i=1; i<s.size(); ++i) {
    //  std::cout << s[i] << " " << pc[i] << std::endl;
    //}

    CP_s_interps.set(n, new MonotCubicInterpolator(std::vector<Real>(s),std::vector<Real>(pc)));
  }
}

void
RockManager::Initialize(const Array<Geometry>& geomArray,
                        const Array<IntVect>&  refRatio)
{
  is_saturated = false;
  is_diffusive = true;
  tensor_diffusion = false;
  use_shifted_Kr_eval = false;
  saturation_threshold_for_vg_Kr = 1;

  int cnt = 0; // Note, incompatible with previous version of structured Amanzi!!!
  CP_models[CP_model_None] = cnt++;
  CP_models[CP_model_vG] = cnt++;
  CP_models[CP_model_BC] = cnt++;

  cnt = 0;
  Kr_models[Kr_model_None] = cnt++;
  Kr_models[Kr_model_vG_Mualem] = cnt++;
  Kr_models[Kr_model_vG_Burdine] = cnt++;
  Kr_models[Kr_model_BC_Mualem] = cnt++;
  Kr_models[Kr_model_BC_Burdine] = cnt++;

  ParmParse pp;
  int nrock = pp.countval("rock");
  if (nrock <= 0) {
    BoxLib::Abort("At least one rock type must be defined.");
  }
  Array<std::string> r_names;  pp.getarr("rock",r_names,0,nrock);

  rock.clear();
  rock.resize(nrock,PArrayManage);
  Array<std::string> material_regions;

  // Scan rock for properties that must be defined for all
  //   if defined for one. 
  bool user_specified_molecular_diffusion_coefficient = false;
  bool user_specified_dispersivity = false;
  bool user_specified_tortuosity = false;
  bool user_specified_specific_storage = false;
  for (int i = 0; i<nrock; i++) {
    const std::string& rname = r_names[i];
    const std::string prefix("rock." + rname);
    ParmParse ppr(prefix.c_str());
    user_specified_molecular_diffusion_coefficient = ppr.countval("molecular_diffusion.val");
    user_specified_tortuosity = ppr.countval("tortuosity.val");
    user_specified_dispersivity = ppr.countval("dispersivity.alphaL");
    user_specified_specific_storage = ppr.countval("specific_storage.val");
  }

  if (is_saturated) {
    user_specified_specific_storage = true; // Will use default if not specified
  }

  bool enable_diffusion = is_diffusive 
    && ( user_specified_molecular_diffusion_coefficient || user_specified_dispersivity);
  
  bool enable_tensor_diffusion = enable_diffusion && user_specified_dispersivity;

  for (int i = 0; i<nrock; i++) {

    const std::string& rname = r_names[i];
    const std::string prefix("rock." + rname);
    ParmParse ppr(prefix.c_str());
        
    static Property::CoarsenRule arith_crsn = Property::Arithmetic;
    static Property::CoarsenRule harm_crsn = Property::ComponentHarmonic;
    static Property::RefineRule pc_refine = Property::PiecewiseConstant;

    Real rdensity = -1; // ppr.get("density",rdensity); // not actually used anywhere

    Real rDmolec = 0;
    Property* Dmolec_func = 0;
    if (user_specified_molecular_diffusion_coefficient) {
      ppr.query("molecular_diffusion.val",rDmolec);
      std::string Dmolec_str = "molecular_diffusion_coefficient";
      Dmolec_func = new ConstantProperty(Dmolec_str,rDmolec,harm_crsn,pc_refine);
    }

    Array<Real> rDispersivity(2,0);
    Property* Dispersivity_func = 0;
    if (user_specified_dispersivity) {
      ppr.query("dispersivity.alphaL",rDispersivity[0]);
      ppr.query("dispersivity.alphaT",rDispersivity[1]);
      std::string Dispersivity_str = "dispersivity";
      Dispersivity_func = new ConstantProperty(Dispersivity_str,rDispersivity,harm_crsn,pc_refine);
    }

    Real rTortuosity = 1;
    Property* Tortuosity_func = 0;
    if (user_specified_tortuosity) {
      ppr.query("tortuosity.val",rTortuosity);
      std::string Tortuosity_str = "tortuosity";
      Tortuosity_func = new ConstantProperty(Tortuosity_str,rTortuosity,harm_crsn,pc_refine);
    }

    Real rSpecificStorage = 0;
    Property* SpecificStorage_func = 0;
    if (user_specified_specific_storage) {
      ppr.query("specific_storage.val",rSpecificStorage);
      std::string SpecificStorage_str = "specific_storage";
      SpecificStorage_func = new ConstantProperty(SpecificStorage_str,rSpecificStorage,arith_crsn,pc_refine);
    }

    Property* phi_func = 0;
    Array<Real> rpvals(1), rptimes;
    Array<std::string> rpforms;
    std::string PorosityValsName = PorosityName+".vals";
    std::string PorosityTimesName = PorosityName+".times";
    std::string PorosityFormsName = PorosityName+".forms";
    if (ppr.countval(PorosityValsName.c_str())) {
      ppr.getarr(PorosityValsName.c_str(),rpvals,0,ppr.countval(PorosityValsName.c_str()));
      int nrpvals = rpvals.size();
      if (nrpvals>1) {
        ppr.getarr(PorosityTimesName.c_str(),rptimes,0,nrpvals);
        ppr.getarr(PorosityFormsName.c_str(),rpforms,0,nrpvals-1);
        TabularFunction pft(rptimes,rpvals,rpforms);
        phi_func = new TabularInTimeProperty(PorosityName.c_str(),pft,arith_crsn,pc_refine);
      }
      else {
        phi_func = new ConstantProperty(PorosityName.c_str(),rpvals[0],arith_crsn,pc_refine);
      }
    } else if (ppr.countval(PorosityName.c_str())) {
      ppr.get(PorosityName.c_str(),rpvals[0]); // FIXME: For backward compatibility
      phi_func = new ConstantProperty(PorosityName.c_str(),rpvals[0],arith_crsn,pc_refine);
    } else {
      BoxLib::Abort(std::string("No porosity function specified for rock: \""+rname).c_str());
    }

    Property* kappa_func = 0;
    Array<Real> rvpvals(1), rhpvals(1), rvptimes(1), rhptimes(1);
    Array<std::string> rvpforms, rhpforms;

    Array<Real> rperm_in(2);
    if (ppr.countval(PermeabilityName.c_str())) {
      ppr.getarr(PermeabilityName.c_str(),rperm_in,0,2);
      rhpvals[0] = rperm_in[0];
      rvpvals[0] = rperm_in[1];
    }
    else {

      std::string PermeabilityVertValName = PermeabilityName+".vertical.vals";
      std::string PermeabilityVertTimesName = PermeabilityName+".vertical.times";
      std::string PermeabilityVertFormsName = PermeabilityName+".vertical.forms";
      std::string PermeabilityHoriValName = PermeabilityName+".horizontal.vals";
      std::string PermeabilityHoriTimesName = PermeabilityName+".horizontal.times";
      std::string PermeabilityHoriFormsName = PermeabilityName+".horizontal.forms";
      int nrvpvals = ppr.countval(PermeabilityVertValName.c_str());
      int nrhpvals = ppr.countval(PermeabilityHoriValName.c_str());
      if (nrvpvals>0 && nrhpvals>0) {
        ppr.getarr(PermeabilityVertValName.c_str(),rvpvals,0,nrvpvals);
        if (nrvpvals>1) {
          ppr.getarr(PermeabilityVertTimesName.c_str(),rvptimes,0,nrvpvals);
          ppr.getarr(PermeabilityVertFormsName.c_str(),rvpforms,0,nrvpvals-1);
        }

        ppr.getarr(PermeabilityHoriValName.c_str(),rhpvals,0,nrhpvals);
        if (nrhpvals>1) {
          ppr.getarr(PermeabilityHoriTimesName.c_str(),rhptimes,0,nrhpvals);
          ppr.getarr(PermeabilityHoriFormsName.c_str(),rhpforms,0,nrhpvals-1);
        }

      } else {
        BoxLib::Abort(std::string("No permeability function specified for rock: \""+rname).c_str());
      }
    }

    // The permeability is specified in mDa.  
    // This needs to be multiplied with 1e-10 to be consistent 
    // with the other units in the code.  What this means is that
    // we will be evaluating the darcy velocity as:
    //
    //  u_Darcy [m/s] = ( kappa [X . mD] / mu [Pa.s] ).Grad(p) [atm/m]
    //
    // where X is the factor necessary to have this formula be dimensionally
    // consistent.  X here is 1.e-10, and can be combined with kappa for the 
    // the moment because no other derived quantities depend directly on the 
    // value of kappa  (NOTE: We will have to know that this is done however
    // if kappa is used as a diagnostic or in some way for a derived quantity).
    //
    for (int j=0; j<rhpvals.size(); ++j) {
      rhpvals[j] *= 1.e-10;
    }
    for (int j=0; j<rvpvals.size(); ++j) {
      rvpvals[j] *= 1.e-10;
    }

    if (rvpvals.size()>1 || rhpvals.size()>1) {
      Array<TabularFunction> pft(2);
      pft[0] = TabularFunction(rhptimes,rhpvals,rhpforms);
      pft[1] = TabularFunction(rvptimes,rvpvals,rvpforms);
      kappa_func = new TabularInTimeProperty(PermeabilityName,pft,harm_crsn,pc_refine);
    }
    else {
      Array<Real> vals(2); vals[0] = rhpvals[0]; vals[1] = rvpvals[0];
      kappa_func = new ConstantProperty(PermeabilityName,vals,harm_crsn,pc_refine);
    }

    // Set old-style values
    Array<Real> rpermeability(BL_SPACEDIM,rvpvals[0]);
    for (int j=0;j<BL_SPACEDIM-1;j++) rpermeability[j] = rhpvals[0];


    // capillary pressure: include cpl_coef, residual_saturation, sigma
    const std::string cpl_prefix(prefix+".cpl");
    ParmParse pp_cpl(cpl_prefix.c_str());
    std::string cpl_model; pp_cpl.get("type",cpl_model); 
    std::map<std::string,int>::const_iterator it = CP_models.find(cpl_model);
    int rcplType = -1;
    int rKrType = -1;

    Array<Real> rcplParam(MAX_CPL_PARAMS);

    if (it != CP_models.end()) {
      rcplType = it->second;

      bool is_vG = (rcplType == CP_models[CP_model_vG]);
      bool is_BC = (rcplType == CP_models[CP_model_BC]);

      if (rcplType == CP_models[CP_model_None]) {
        rKrType = Kr_models[Kr_model_None];
      }
      else if (is_vG || is_BC)
      {
        Real m, lambda;
        if (is_vG) {
          pp_cpl.get("m",m);
          if (m <= 0) {
            if (ParallelDescriptor::IOProcessor()) {
              std::cerr << "Invalid m (= " << m << " ) for Capillary Pressure model in material: \"" << rname << "\"" << std::endl;
            } BoxLib::Abort();
          }
        } else {
          pp_cpl.get("lambda",lambda);
          if (lambda <= 0) {
            if (ParallelDescriptor::IOProcessor()) {
              std::cerr << "Invalid lambda (= " << lambda << " ) for Capillary Pressure model in material: \"" << rname << "\"" << std::endl;
            } BoxLib::Abort();
          }
        }          
        Real Sr; pp_cpl.get("Sr",Sr);
        if (Sr < 0 || Sr > 1) {
          if (ParallelDescriptor::IOProcessor()) {
            std::cerr << "Invalid Sr (= " << Sr << " ) for Capillary Pressure model in material: \"" << rname << "\"" << std::endl;
          } BoxLib::Abort();
        }

        Real alpha; pp_cpl.get("alpha",alpha);
        if (alpha < 0 ) {
          if (ParallelDescriptor::IOProcessor()) {
            std::cerr << "Invalid alpha (= " << m << " ) for Capillary Pressure model in material: \"" << rname << "\"" << std::endl;
          } BoxLib::Abort();
        }

        Real Kr_ell; ppr.get("Kr_ell",Kr_ell);
        std::string Kr_model; ppr.get("Kr_model", Kr_model);
        std::string Kr_full_model_name = cpl_model + "_" + Kr_model;
        std::map<std::string,int>::const_iterator itKr = Kr_models.find(Kr_full_model_name);
        if (it != Kr_models.end()) {
          rKrType = itKr->second;
        }
        else {
          if (ParallelDescriptor::IOProcessor()) {
            std::cerr << "Invalid Kr model (= \"" << Kr_model
                      << "\") for Relative Permeability with Capillary Pressure model (\"" << cpl_model
                      << "\") in material: \"" << rname << "\""<< std::endl;
          } BoxLib::Abort();
        }

        // Finally, load array of Real numbers for this model
        rcplParam[CPL_MODEL_ID] = (Real)rcplType;
        if (is_vG) {
          rcplParam[VG_M]     = m;
          rcplParam[VG_ALPHA] = alpha;
          rcplParam[VG_SR]    = Sr;
          rcplParam[VG_ELL]   = Kr_ell;
          rcplParam[VG_KR_MODEL_ID]  = (Real)rKrType;
        }
        else {
          rcplParam[BC_LAMBDA] = lambda;
          rcplParam[BC_ALPHA]  = alpha;
          rcplParam[BC_SR]     = Sr;
          rcplParam[BC_ELL]    = Kr_ell;
          rcplParam[BC_KR_MODEL_ID]  = (Real)rKrType;
        }
      }
      else {
        if (ParallelDescriptor::IOProcessor()) {
          std::cerr << "Unknown capillary pressure (" << cpl_model << ") model for " << rname << std::endl;
        }
      }
    }

    Property* cpl_func = new ConstantProperty(CapillaryPressureName,rcplParam,arith_crsn,pc_refine);

    // relative permeability: include kr_coef, residual_saturation
    int rkrType = 0;  ppr.query("kr_type",rkrType);
    Array<Real> rkrParam;
    if (rkrType > 0) {
      ppr.getarr("kr_param",rkrParam,0,ppr.countval("kr_param"));
    }

    Array<Real> krPt(rkrParam.size()+1);
    krPt[0] = Real(rkrType);
    for (int j=0; j<rkrParam.size(); ++j) {
      krPt[j+1] = rkrParam[j];
    }
    std::string krel_str = "relative_permeability";
    Property* krel_func = new ConstantProperty(krel_str,krPt,arith_crsn,pc_refine);

    Array<std::string> region_names;
    ppr.getarr("regions",region_names,0,ppr.countval("regions"));
    Array<const Region*> rregions = region_manager->RegionPtrArray(region_names);
    for (int j=0; j<region_names.size(); ++j) {
      material_regions.push_back(region_names[j]);
    }

    if (ppr.countval("porosity_dist_param")>0) {
      BoxLib::Abort("porosity_dist_param not currently supported");
    }

    std::string porosity_dist="uniform"; ppr.query("porosity_dist",porosity_dist);
    Array<Real> rporosity_dist_param;
    if (porosity_dist!="uniform") {
      BoxLib::Abort("porosity_dist != uniform not currently supported");
      ppr.getarr("porosity_dist_param",rporosity_dist_param,
                 0,ppr.countval("porosity_dist_param"));
    }
        
    std::string permeability_dist="uniform"; ppr.get("permeability_dist",permeability_dist);
    Array<Real> rpermeability_dist_param;
    if (permeability_dist != "uniform") {
      ppr.getarr("permeability_dist_param",rpermeability_dist_param,
                 0,ppr.countval("permeability_dist_param"));
    }

    std::vector<Property*> properties;
    properties.push_back(phi_func);
    properties.push_back(kappa_func);
    if (Dmolec_func) {
      properties.push_back(Dmolec_func);
    }
    if (Dispersivity_func) {
      properties.push_back(Dispersivity_func);
    }
    if (Tortuosity_func) {
      properties.push_back(Tortuosity_func);
    }
    if (SpecificStorage_func) {
      properties.push_back(SpecificStorage_func);
    }
    properties.push_back(krel_func);
    properties.push_back(cpl_func);
    rock.set(i,new Material(rname,rregions,properties));
    delete phi_func;
    delete kappa_func;
    delete Dmolec_func;
    delete Tortuosity_func;
    delete Dispersivity_func;
    delete SpecificStorage_func;
    delete krel_func;
    delete cpl_func;
  }

#if 0
  // Read rock parameters associated with chemistry
  using_sorption = false;
  aux_chem_variables.clear();

  if (do_tracer_chemistry>0)
  {
    ParmParse ppm("mineral");
    nminerals = ppm.countval("minerals");
    minerals.resize(nminerals);
    if (nminerals>0) {
      ppm.getarr("minerals",minerals,0,nminerals);
    }

    ParmParse pps("sorption_site");
    nsorption_sites = pps.countval("sorption_sites");
    sorption_sites.resize(nsorption_sites);
    if (nsorption_sites>0) {
      pps.getarr("sorption_sites",sorption_sites,0,nsorption_sites);
    }

    ICParmPair sorption_isotherm_options;
    sorption_isotherm_options[          "Kd"] = 0;
    sorption_isotherm_options[  "Langmuir_b"] = 0;
    sorption_isotherm_options["Freundlich_n"] = 1;
	
    for (int k=0; k<tNames.size(); ++k) {
      for (ICParmPair::const_iterator it=sorption_isotherm_options.begin();
           it!=sorption_isotherm_options.end(); ++it) {
        const std::string& str = it->first;
        bool found = false;
        for (int i=0; i<nrock; ++i) {
          const std::string prefix("rock."+r_names[i]+".Sorption_Isotherms."+tNames[k]);
          ParmParse pprs(prefix.c_str());
          if (pprs.countval(str.c_str())) {
            pprs.get(str.c_str(),sorption_isotherm_ics[r_names[i]][tNames[k]][str]);
            found = true;
          }
        }
	    
        if (found) {
          using_sorption = true;
          nsorption_isotherms = ntracers;
          for (int i=0; i<nrock; ++i) {
            if (sorption_isotherm_ics[r_names[i]][tNames[k]].count(str) == 0) {
              sorption_isotherm_ics[r_names[i]][tNames[k]][str] = it->second; // set to default value
            }
          }
          const std::string label = str+"_"+tNames[k];
          if (aux_chem_variables.find(label) == aux_chem_variables.end()) {
            sorption_isotherm_label_map[tNames[k]][str] = aux_chem_variables.size();
            aux_chem_variables[label]=aux_chem_variables.size()-1;
          }
        }
      }
    }

    ICParmPair cation_exchange_options;
    cation_exchange_options["Cation_Exchange_Capacity"] = 0;
    {
      for (ICParmPair::const_iterator it=cation_exchange_options.begin(); it!=cation_exchange_options.end(); ++it) {
        const std::string& str = it->first;
        bool found = false;
        for (int i=0; i<nrock; ++i) {
          const std::string prefix("rock."+r_names[i]);
          ParmParse pprs(prefix.c_str());
          if (pprs.countval(str.c_str())) {
            pprs.get(str.c_str(),cation_exchange_ics[r_names[i]]);
            found = true;
          }
        }
	    
        if (found) {
          using_sorption = true;
          ncation_exchange = 1;
          for (int i=0; i<nrock; ++i) {
            if (cation_exchange_ics.count(r_names[i]) == 0) {
              cation_exchange_ics[r_names[i]] = it->second; // set to default value
            }
          }

          const std::string label = str;
          if (aux_chem_variables.find(label) == aux_chem_variables.end())  {
            cation_exchange_label_map[str] = aux_chem_variables.size();
            aux_chem_variables[label]=aux_chem_variables.size()-1;
          }
          //std::cout << "****************** cation_exchange_ics[" << r_names[i] << "] = " 
          //	  << cation_exchange_ics[r_names[i]] << std::endl;
        }
      }
    }

    ICParmPair mineralogy_options;
    mineralogy_options[      "Volume_Fraction"] = 0;
    mineralogy_options["Specific_Surface_Area"] = 0;
    for (int k=0; k<minerals.size(); ++k) {
      for (ICParmPair::const_iterator it=mineralogy_options.begin(); it!=mineralogy_options.end(); ++it) {
        const std::string& str = it->first;
        bool found = false;
        for (int i=0; i<nrock; ++i) {
          const std::string prefix("rock."+r_names[i]+".mineralogy."+minerals[k]);
          ParmParse pprs(prefix.c_str());
          if (pprs.countval(str.c_str())) {
            pprs.get(str.c_str(),mineralogy_ics[r_names[i]][minerals[k]][str]);
            found = true;
          }
        }
	    
        if (found) {
          using_sorption = true;
          for (int i=0; i<nrock; ++i) {
            if (mineralogy_ics[r_names[i]][minerals[k]].count(str) == 0) {
              mineralogy_ics[r_names[i]][minerals[k]][str] = it->second; // set to default value
            }
          }
          //std::cout << "****************** mineralogy_ics[" << r_names[i] << "][" << minerals[k] 
          //	  << "][" << str << "] = " << mineralogy_ics[r_names[i]][minerals[k]][str] 
          //	  << std::endl;

          const std::string label = str+"_"+minerals[k];
          if (aux_chem_variables.find(label) == aux_chem_variables.end()) {
            mineralogy_label_map[minerals[k]][str] = aux_chem_variables.size();
            aux_chem_variables[label]=aux_chem_variables.size()-1;
          }
        }
      }
    }

    ICParmPair complexation_options;
    complexation_options["Site_Density"] = 0;
    for (int k=0; k<sorption_sites.size(); ++k) {
      for (ICParmPair::const_iterator it=complexation_options.begin(); it!=complexation_options.end(); ++it) {
        const std::string& str = it->first;
        bool found = false;
        for (int i=0; i<nrock; ++i) {
          const std::string prefix("rock."+r_names[i]+".Surface_Complexation_Sites."+sorption_sites[k]);
          ParmParse pprs(prefix.c_str());
          if (pprs.countval(str.c_str())) {
            pprs.get(str.c_str(),surface_complexation_ics[r_names[i]][sorption_sites[k]][str]);
            found = true;
          }
        }
	    
        if (found) {
          using_sorption = true;
          for (int i=0; i<nrock; ++i) {
            if (surface_complexation_ics[r_names[i]][sorption_sites[k]].count(str) == 0) {
              surface_complexation_ics[r_names[i]][sorption_sites[k]][str] = it->second; // set to default value
            }
          }
          //std::cout << "****************** surface_complexation_ics[" << r_names[i] << "][" << sorption_sites[k] 
          //	  << "][" << str << "] = " << sorption_isotherm_ics[r_names[i]][sorption_sites[k]][str] 
          //	  << std::endl;
	      
          const std::string label = str+"_"+sorption_sites[k];
          if (aux_chem_variables.find(label) == aux_chem_variables.end()) {
            surface_complexation_label_map[sorption_sites[k]][str] = aux_chem_variables.size();
            aux_chem_variables[label]=aux_chem_variables.size()-1;
          }
        }
      }
    }

    if (using_sorption) 
    {
      ICParmPair sorption_chem_options; // these are domain-wide, specified per solute
      sorption_chem_options["Total_Sorbed"] = 1.e-40;
      for (int k=0; k<tNames.size(); ++k) {
        for (ICParmPair::const_iterator it=sorption_chem_options.begin(); it!=sorption_chem_options.end(); ++it) {
          const std::string& str = it->first;
          const std::string prefix("tracer."+tNames[k]+".Initial_Condition."+str);
          ParmParse pprs(prefix.c_str());
          sorption_chem_ics[tNames[k]][str] = it->second; // set to default value
          pprs.query(str.c_str(),sorption_chem_ics[tNames[k]][str]);                      
          //std::cout << "****************** sorption_chem_ics[" << tNames[k] 
          //              << "][" << str << "] = " << sorption_chem_ics[tNames[k]][str] << std::endl;
          const std::string label = str+"_"+tNames[k];
          if (aux_chem_variables.find(label) == aux_chem_variables.end()){
            sorption_chem_label_map[tNames[k]][str] = aux_chem_variables.size();
            aux_chem_variables[label]=aux_chem_variables.size()-1;
          }
        }
      }
    }
  }
#endif

  materialFiller = new MatFiller(geomArray,refRatio,rock);

  pp.query("Use_Shifted_Kr_Eval",use_shifted_Kr_eval);
  pp.query("Saturation_Threshold_For_Kr",saturation_threshold_for_vg_Kr);

  if (use_shifted_Kr_eval && saturation_threshold_for_vg_Kr>1) {
    if (ParallelDescriptor::IOProcessor()) {
      std::cout << "WARNING: Reducing Saturation_Threshold_For_vg_Kr to 1!" << std::endl;
    }
    saturation_threshold_for_vg_Kr = 1;
  }

#if 0
  FORT_KR_INIT(&saturation_threshold_for_vg_Kr,
               &use_shifted_Kr_eval);
#endif

}

bool
RockManager::CanDerive(const std::string& property_name)
{
  return materialFiller->CanDerive(property_name);
}

bool
RockManager::GetProperty(Real               time,
                         int                level,
                         MultiFab&          mf,
                         const std::string& pname,
                         int                dComp,
                         int                nGrow,
                         void*              ctx,
                         bool               ignore_mixed)
{
  return materialFiller->SetProperty(time,level,mf,pname,dComp,nGrow,ctx,ignore_mixed);
}

bool
RockManager::GetPropertyDirect(Real               t,
                               int                level,
                               FArrayBox&         fab,
                               const Box&         box,
                               const std::string& pname,
                               int                dComp,
                               void*              ctx)
{
  return materialFiller->SetPropertyDirect(t,level,fab,box,pname,dComp,ctx);
}

void
RockManager::GetMaterialID(int level, iMultiFab& mf, int nGrow, bool ignore_mixed)
{
  return materialFiller->SetMaterialID(level,mf,nGrow,ignore_mixed);
}


  // Capillary Pressure (given Saturation)
void
RockManager::CapillaryPressure(const Real* saturation, int* matID, Real time, Real* capillaryPressure, int npts)
{
  Array<Array<int> > mat_pts = SortPtsByMaterial(matID,npts);

  // Make temp structure to interact with Property interface
  int nComp = materialFiller->NComp(CapillaryPressureName);
  static IntVect iv(D_DECL(0,0,0));
  static Box bx(iv,iv);
  FArrayBox pc_params(bx,nComp);
  int level=0; //not really used
  int dComp=0;

  for (int n=0; n<rock.size(); ++n) {
    const Property* p = rock[n].Prop(CapillaryPressureName); BL_ASSERT(p!=0);
    p->eval(time,level,bx,pc_params,dComp);
    bool is_vG = (pc_params(iv,CPL_MODEL_ID) == (Real)CP_models[CP_model_vG]);
    bool is_BC = (pc_params(iv,CPL_MODEL_ID) == (Real)CP_models[CP_model_BC]);

    if (is_vG) {
      Real m      = pc_params(iv,VG_M);
      Real alphaI = 1/pc_params(iv,VG_ALPHA);
      Real Sr     = pc_params(iv,VG_SR);
      Real b      = -1/m;
      Real omm    = 1-m;
      Real omSrI  = 1/(1-Sr);

      for (int i=0, End=mat_pts[n].size(); i<End; ++i) {
        int idx = mat_pts[n][i];
        Real s = std::min(1.0, std::max(0.0, saturation[idx]));
        
        Real seff = (s - Sr)*omSrI;
        capillaryPressure[idx] = alphaI * std::pow(std::pow(seff,b) - 1,omm);

#if 0
        if (is_initialized) {
          Real intVal = CP_s_interps[n](saturation[idx]);
          std::cout << " " << (intVal - capillaryPressure[idx])/ capillaryPressure[idx] << std::endl;
        }
#endif
      }
    }
    else if (is_BC) {
      Real mLambdaI = -1/pc_params(iv,BC_LAMBDA);
      Real alphaI   = 1/pc_params(iv,BC_ALPHA);
      Real Sr       = pc_params(iv,BC_SR);
      Real omSrI    = 1/(1-Sr);

      for (int i=0, End=mat_pts[n].size(); i<End; ++i) {
        int idx = mat_pts[n][i];
        Real s = std::min(1.0, std::max(0.0, saturation[idx]));
        
        Real seff = (s - Sr)*omSrI;
        capillaryPressure[idx] = alphaI * std::pow(seff,mLambdaI);
      }
    }
    else {
      if (ParallelDescriptor::IOProcessor()) {
        std::cerr << "Invalid Capillary Presure model " << std::endl;
      } BoxLib::Abort();
    }
  }
}

void
RockManager::CapillaryPressure(const MultiFab&  saturation,
                               const iMultiFab& matID,
                               Real             time,
                               MultiFab&        capillaryPressure,
                               int              sComp,
                               int              dComp,
                               int              nGrow)
{
  BL_ASSERT(saturation.boxArray() == capillaryPressure.boxArray());
  BL_ASSERT(sComp < saturation.nComp()  && dComp < capillaryPressure.nComp());
  BL_ASSERT(nGrow <= saturation.nGrow() && nGrow <= capillaryPressure.nGrow());
  FArrayBox s, pc;
  IArrayBox id;
  for (MFIter mfi(saturation); mfi.isValid(); ++mfi) {
    Box box = Box(mfi.validbox()).grow(nGrow);
    s.resize(box,1); s.copy(saturation[mfi],box,sComp,box,0,1);
    id.resize(box,1); id.copy(matID[mfi],box,0,box,0,1);
    pc.resize(box,1);
    CapillaryPressure(s.dataPtr(),id.dataPtr(),time,pc.dataPtr(),box.numPts());
    capillaryPressure[mfi].copy(pc,box,0,box,dComp,1);
  }
}

Array<Array<int> >
RockManager::SortPtsByMaterial(int* matID, int npts)
{
  Array<Array<int> > mat_pts(rock.size());
  for (int i=0; i<npts; ++i) {
    mat_pts[matID[i]].push_back(i);
  }
  return mat_pts;
}

// Inverse Capillary Pressure (Saturation given Capillary Pressure)
void
RockManager::InverseCapillaryPressure(const Real* capillaryPressure, int* matID, Real time, Real* saturation, int npts)
{
  Array<Array<int> > mat_pts = SortPtsByMaterial(matID,npts);

  // Make temp structure to interact with Property interface
  int nComp = materialFiller->NComp(CapillaryPressureName);
  static IntVect iv(D_DECL(0,0,0));
  static Box bx(iv,iv);
  FArrayBox pc_params(bx,nComp);
  int level=0; //not really used
  int dComp=0;

  for (int j=0; j<rock.size(); ++j) {
    const Property* p = rock[j].Prop(CapillaryPressureName); BL_ASSERT(p!=0);
    p->eval(time,level,bx,pc_params,dComp);
    bool is_vG = (pc_params(iv,CPL_MODEL_ID) == (Real)CP_models[CP_model_vG]);
    bool is_BC = (pc_params(iv,CPL_MODEL_ID) == (Real)CP_models[CP_model_BC]);

    if (is_vG) {
      Real m     = pc_params(iv,VG_M);
      Real alpha = pc_params(iv,VG_ALPHA);
      Real Sr    = pc_params(iv,VG_SR);
      Real n     = 1./(1-m);
      
      for (int i=0, End=mat_pts[j].size(); i<End; ++i) {
        int idx = mat_pts[j][i];
        Real seff = std::pow( std::pow(alpha*capillaryPressure[idx],n) + 1, -m);
        saturation[idx] = seff*(1 - Sr) + Sr;
      }
    }
    else if (is_BC) {
      Real mLambda = -pc_params(iv,BC_LAMBDA);
      Real alpha   = pc_params(iv,BC_ALPHA);
      Real Sr      = pc_params(iv,BC_SR);
      
      for (int i=0, End=mat_pts[j].size(); i<End; ++i) {
        int idx = mat_pts[j][i];
        Real seff = std::pow(alpha*capillaryPressure[idx],mLambda);
        saturation[idx] = seff*(1 - Sr) + Sr;
      }
    }
    else {
      if (ParallelDescriptor::IOProcessor()) {
        std::cerr << "Invalid Capillary Presure model " << std::endl;
      } BoxLib::Abort();
    }
  }
}

void
RockManager::InverseCapillaryPressure(const MultiFab&  capillaryPressure,
                                      const iMultiFab& matID,
                                      Real             time,
                                      MultiFab&        saturation,
                                      int              sComp,
                                      int              dComp,
                                      int              nGrow)
{
  BL_ASSERT(saturation.boxArray() == capillaryPressure.boxArray());
  BL_ASSERT(dComp < saturation.nComp()  && sComp < capillaryPressure.nComp());
  BL_ASSERT(nGrow <= saturation.nGrow() && nGrow <= capillaryPressure.nGrow());
  FArrayBox s, pc;
  IArrayBox id;
  for (MFIter mfi(saturation); mfi.isValid(); ++mfi) {
    Box box = Box(mfi.validbox()).grow(nGrow);
    pc.resize(box,1); pc.copy(capillaryPressure[mfi],box,sComp,box,0,1);
    id.resize(box,1); id.copy(matID[mfi],box,0,box,0,1);
    s.resize(box,1);
    InverseCapillaryPressure(pc.dataPtr(),id.dataPtr(),time,s.dataPtr(),box.numPts());
    saturation[mfi].copy(s,box,0,box,dComp,1);
  }
}

// D (Saturation) / D(CapillaryPressure)
void
RockManager::DInverseCapillaryPressure(const Real* saturation, int* matID, Real time, Real* DsaturationDpressure, int npts)
{
  Array<Array<int> > mat_pts = SortPtsByMaterial(matID,npts);

  // Make temp structure to interact with Property interface
  int nComp = materialFiller->NComp(CapillaryPressureName);
  static IntVect iv(D_DECL(0,0,0));
  static Box bx(iv,iv);
  FArrayBox pc_params(bx,nComp);
  int level=0; //not really used
  int dComp=0;

  for (int j=0; j<rock.size(); ++j) {
    const Property* p = rock[j].Prop(CapillaryPressureName); BL_ASSERT(p!=0);
    p->eval(time,level,bx,pc_params,dComp);
    bool is_vG = (pc_params(iv,CPL_MODEL_ID) == (Real)CP_models[CP_model_vG]);
    bool is_BC = (pc_params(iv,CPL_MODEL_ID) == (Real)CP_models[CP_model_BC]);

    if (is_vG) {
      Real m     = pc_params(iv,VG_M);
      Real alpha = pc_params(iv,VG_ALPHA);
      Real Sr    = pc_params(iv,VG_SR);
      Real n     = 1./(1-m);
      Real b     = -1/m;
      Real fac   = - (1 - Sr)*alpha*m*n;
      Real omSrI = 1/(1 - Sr);

      for (int i=0, End=mat_pts[j].size(); i<End; ++i) {
        int idx = mat_pts[j][i];
        Real seff = (saturation[idx] - Sr)*omSrI;
        Real sb = std::pow(seff,b);
        DsaturationDpressure[idx] = fac*std::pow(sb-1,m)*seff/sb;
      }
    }
    else if (is_BC) {
      Real lambda = pc_params(iv,BC_LAMBDA);
      Real alpha  = pc_params(iv,BC_ALPHA);
      Real Sr     = pc_params(iv,BC_SR);
      Real fac    = -alpha*lambda;
      Real oplI   = 1+1/lambda;
      Real omSrI = 1/(1 - Sr);

      for (int i=0, End=mat_pts[j].size(); i<End; ++i) {
        int idx = mat_pts[j][i];
        Real seff = (saturation[idx] - Sr)*omSrI;
        DsaturationDpressure[idx] = fac*std::pow(seff,oplI);
      }
    }
  }
}

void
RockManager::DInverseCapillaryPressure(const MultiFab&  saturation,
                                       const iMultiFab& matID,
                                       Real             time,
                                       MultiFab&        DsaturationDpressure,
                                       int              sComp,
                                       int              dComp,
                                       int              nGrow)
{
  BL_ASSERT(saturation.boxArray() == DsaturationDpressure.boxArray());
  BL_ASSERT(dComp < saturation.nComp()  && sComp < DsaturationDpressure.nComp());
  BL_ASSERT(nGrow <= saturation.nGrow() && nGrow <= DsaturationDpressure.nGrow());
  FArrayBox s, dsdp;
  IArrayBox id;
  for (MFIter mfi(saturation); mfi.isValid(); ++mfi) {
    Box box = Box(mfi.validbox()).grow(nGrow);
    s.resize(box,1); s.copy(saturation[mfi],box,sComp,box,0,1);
    id.resize(box,1); id.copy(matID[mfi],box,0,box,0,1);
    dsdp.resize(box,1);
    DInverseCapillaryPressure(s.dataPtr(),id.dataPtr(),time,dsdp.dataPtr(),box.numPts());
    DsaturationDpressure[mfi].copy(dsdp,box,0,box,dComp,1);
  }
}

void
RockManager::RelativePermeability(const Real* saturation, int* matID, Real time, Real* kappa, int npts)
{
  Array<Array<int> > mat_pts = SortPtsByMaterial(matID,npts);

  // Make temp structure to interact with Property interface
  int nComp = materialFiller->NComp(CapillaryPressureName);
  static IntVect iv(D_DECL(0,0,0));
  static Box bx(iv,iv);
  FArrayBox pc_params(bx,nComp);
  int level=0; //not really used
  int dComp=0;

  for (int j=0; j<rock.size(); ++j) {
    const Property* p = rock[j].Prop(CapillaryPressureName); BL_ASSERT(p!=0);
    p->eval(time,level,bx,pc_params,dComp);
    bool is_vG = (pc_params(iv,CPL_MODEL_ID) == (Real)CP_models[CP_model_vG]);
    bool is_BC = (pc_params(iv,CPL_MODEL_ID) == (Real)CP_models[CP_model_BC]);

    if (is_vG) {

      Real m     = pc_params(iv,VG_M);
      Real alpha = pc_params(iv,VG_ALPHA);
      Real Sr    = pc_params(iv,VG_SR);
      Real ell   = pc_params(iv,VG_ELL);
      Real omSrI = 1/(1 - Sr);

      bool is_Mualem  = (pc_params(iv,VG_KR_MODEL_ID) == (Real)Kr_models[Kr_model_vG_Mualem]);
      bool is_Burdine = (pc_params(iv,VG_KR_MODEL_ID) == (Real)Kr_models[Kr_model_vG_Burdine]);

      Real oom = 1./m;
      if (is_Mualem) {
        for (int i=0, End=mat_pts[j].size(); i<End; ++i) {
          int idx = mat_pts[j][i];
          Real seff = (saturation[idx] - Sr)*omSrI;
          Real tmp = 1 - std::pow(1 - std::pow(seff,oom),m);
          kappa[idx] = std::pow(seff,ell) * tmp * tmp;
        }
      }
      else if (is_Burdine) {
        for (int i=0, End=mat_pts[j].size(); i<End; ++i) {
          int idx = mat_pts[j][i];
          Real seff = (saturation[idx] - Sr)*omSrI;
          Real tmp = 1 - std::pow(1 - std::pow(seff,oom),m);
          kappa[idx] = std::pow(seff,ell) * tmp;
        }
      }

    } else if (is_BC) {

      Real lambda = pc_params(iv,BC_LAMBDA);
      Real alpha  = pc_params(iv,BC_ALPHA);
      Real Sr     = pc_params(iv,BC_SR);
      Real ell    = pc_params(iv,BC_ELL);
      Real omSrI  = 1/(1 - Sr);

      bool is_Mualem  = (pc_params(iv,VG_KR_MODEL_ID) == (Real)Kr_models[Kr_model_vG_Mualem]);
      bool is_Burdine = (pc_params(iv,VG_KR_MODEL_ID) == (Real)Kr_models[Kr_model_vG_Burdine]);

      BL_ASSERT(is_Mualem || is_Burdine);
      Real f = (is_Mualem ? ell + 2 + 2/lambda : ell + 1 + 2/lambda);

      for (int i=0, End=mat_pts[j].size(); i<End; ++i) {
        int idx = mat_pts[j][i];
        Real seff = (saturation[idx] - Sr)*omSrI;
        kappa[idx] = std::pow(seff,f);
      }

    }
  }
}

void
RockManager::RelativePermeability(const MultiFab&  saturation,
                                  const iMultiFab& matID,
                                  Real             time,
                                  MultiFab&        kappa,
                                  int              sComp,
                                  int              dComp,
                                  int              nGrow)
{
  BL_ASSERT(saturation.boxArray() == kappa.boxArray());
  BL_ASSERT(dComp < saturation.nComp()  && sComp < kappa.nComp());
  BL_ASSERT(nGrow <= saturation.nGrow() && nGrow <= kappa.nGrow());
  FArrayBox s, k;
  IArrayBox id;
  for (MFIter mfi(saturation); mfi.isValid(); ++mfi) {
    Box box = Box(mfi.validbox()).grow(nGrow);
    s.resize(box,1); s.copy(saturation[mfi],box,sComp,box,0,1);
    id.resize(box,1); id.copy(matID[mfi],box,0,box,0,1);
    k.resize(box,1);
    RelativePermeability(s.dataPtr(),id.dataPtr(),time,k.dataPtr(),box.numPts());
    kappa[mfi].copy(k,box,0,box,dComp,1);
  }
}


