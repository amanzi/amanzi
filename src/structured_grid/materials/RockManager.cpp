
#include <RockManager.H>

#include <ParmParse.H>

RockManager::RockManager(const RegionManager*   _region_manager,
                         const Array<Geometry>& geomArray,
                         const Array<IntVect>&  refRatio)
  : region_manager(_region_manager)
{
  Initialize(geomArray,refRatio);
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

  ParmParse pp("rock");
  int nrock = pp.countval("rock");
  if (nrock <= 0) {
    BoxLib::Abort("At least one rock type must be defined.");
  }
  Array<std::string> r_names;  pp.getarr("rock",r_names,0,nrock);

  materials.clear();
  materials.resize(nrock,PArrayManage);
  Array<std::string> material_regions;

  // Scan materials for properties that must be defined for all
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
    std::string phi_str = "porosity";
    Array<Real> rpvals(1), rptimes;
    Array<std::string> rpforms;
    if (ppr.countval("porosity.vals")) {
      ppr.getarr("porosity.vals",rpvals,0,ppr.countval("porosity.vals"));
      int nrpvals = rpvals.size();
      if (nrpvals>1) {
        ppr.getarr("porosity.times",rptimes,0,nrpvals);
        ppr.getarr("porosity.forms",rpforms,0,nrpvals-1);
        TabularFunction pft(rptimes,rpvals,rpforms);
        phi_func = new TabularInTimeProperty(phi_str,pft,arith_crsn,pc_refine);
      }
      else {
        phi_func = new ConstantProperty(phi_str,rpvals[0],arith_crsn,pc_refine);
      }
    } else if (ppr.countval("porosity")) {
      ppr.get("porosity",rpvals[0]); // FIXME: For backward compatibility
      phi_func = new ConstantProperty(phi_str,rpvals[0],arith_crsn,pc_refine);
    } else {
      BoxLib::Abort(std::string("No porosity function specified for rock: \""+rname).c_str());
    }

    Property* kappa_func = 0;
    std::string kappa_str = "permeability";
    Array<Real> rvpvals(1), rhpvals(1), rvptimes(1), rhptimes(1);
    Array<std::string> rvpforms, rhpforms;

    Array<Real> rperm_in(2);
    if (ppr.countval("permeability")) {
      ppr.getarr("permeability",rperm_in,0,2);
      rhpvals[0] = rperm_in[0];
      rvpvals[0] = rperm_in[1];
    }
    else {

      int nrvpvals = ppr.countval("permeability.vertical.vals");
      int nrhpvals = ppr.countval("permeability.horizontal.vals");
      if (nrvpvals>0 && nrhpvals>0) {
        ppr.getarr("permeability.vertical.vals",rvpvals,0,nrvpvals);
        if (nrvpvals>1) {
          ppr.getarr("permeability.vertical.times",rvptimes,0,nrvpvals);
          ppr.getarr("permeability.vertical.forms",rvpforms,0,nrvpvals-1);
        }

        ppr.getarr("permeability.horizontal.vals",rhpvals,0,nrhpvals);
        if (nrhpvals>1) {
          ppr.getarr("permeability.horizontal.times",rhptimes,0,nrhpvals);
          ppr.getarr("permeability.horizontal.forms",rhpforms,0,nrhpvals-1);
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
      kappa_func = new TabularInTimeProperty(kappa_str,pft,harm_crsn,pc_refine);
    }
    else {
      Array<Real> vals(2); vals[0] = rhpvals[0]; vals[1] = rvpvals[0];
      kappa_func = new ConstantProperty(kappa_str,vals,harm_crsn,pc_refine);
    }

    // Set old-style values
    Array<Real> rpermeability(BL_SPACEDIM,rvpvals[0]);
    for (int j=0;j<BL_SPACEDIM-1;j++) rpermeability[j] = rhpvals[0];

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

    // capillary pressure: include cpl_coef, residual_saturation, sigma
    int rcplType = 0;  ppr.query("cpl_type", rcplType);
    Array<Real> rcplParam;
    if (rcplType > 0) {
      ppr.getarr("cpl_param",rcplParam,0,ppr.countval("cpl_param"));
    }
    Array<Real> cplPt(rcplParam.size()+1);
    cplPt[0] = Real(rcplType);
    for (int j=0; j<rcplParam.size(); ++j) {
      cplPt[j+1] = rcplParam[j];
    }
    std::string cpl_str = "capillary_pressure";
    Property* cpl_func = new ConstantProperty(cpl_str,cplPt,arith_crsn,pc_refine);

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
    materials.set(i,new Material(rname,rregions,properties));
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

  materialFiller = new MatFiller(geomArray,refRatio,materials);

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
