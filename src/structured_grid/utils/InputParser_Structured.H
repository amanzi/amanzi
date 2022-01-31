#ifndef Input_Parser_Structured_H_
#define Input_Parser_Structured_H_
#include "Teuchos_ParameterList.hpp"

namespace Amanzi {
  namespace AmanziInput {
    
    using Teuchos::ParameterList;
    using Teuchos::Array;

    static int ndim;
    static double simulation_time;
    static Array<double> domlo,domhi;      
    static double geometry_eps;

    static Array<std::string> user_derive_list;

    ParameterList convert_to_structured(const ParameterList& parameter_list);
    std::string underscore(const std::string& instring);

      struct ICBCFunc
      {
          ICBCFunc() {};
          ICBCFunc(const ParameterList& plist,
                   const std::string& Amanzi_type,
                   const std::string& label)
              : Amanzi_type(Amanzi_type),
                label(label),
                plist(plist) {} 
          
          const ParameterList& PList() const {return plist;}
          const std::string& Amanzi_Type() const {return Amanzi_type;}
          const std::string& Label() const {return label;}
          
      protected:
          std::string Amanzi_type, label;
          ParameterList plist;
      };

      struct CompFunc
      {
          typedef std::map<std::string,ICBCFunc> ICBCFuncMap;
          ICBCFuncMap& getICBCFuncMap() {return soluteFunc_map;}
          ICBCFunc& operator [] (const std::string& str) {return soluteFunc_map[str];}
          int size() const {return soluteFunc_map.size();}

      protected:
          ICBCFuncMap soluteFunc_map;
      };

      struct PhaseFunc
      {
          typedef std::map<std::string,CompFunc> CompFuncMap;

          CompFuncMap& getCompFuncMap() {return compFunc_map;}
          CompFunc& operator [] (const std::string& compName) {return compFunc_map[compName];}
          int size() const {return compFunc_map.size();}
          
          CompFuncMap compFunc_map;
      };
      
      struct StateFunc
      {
          StateFunc() {}
          StateFunc(const std::string&        label,
                    const std::string&        Amanzi_type,
                    const ParameterList&      func_plist,
                    const Array<std::string>& regions)
              : func_plist(func_plist),
                Amanzi_type(Amanzi_type),
                label(label),
                regions(regions) {}


          typedef std::map<std::string,PhaseFunc> PhaseFuncMap;

          PhaseFuncMap& getPhaseFuncMap() {return phaseFunc_map;}
          PhaseFunc& operator [] (const std::string& phaseName) {return phaseFunc_map[phaseName];}
          int size() const {return phaseFunc_map.size();}

          const Array<std::string>& Regions() const {return regions;}
          const std::string& Amanzi_Type() const {return Amanzi_type;}
          const ParameterList& FuncPList() const {return func_plist;}

      protected:
          ParameterList func_plist;
          std::string Amanzi_type, label;
          PhaseFuncMap phaseFunc_map;
          Array<std::string> regions;
      };
      
      struct TRACER
      {
        TRACER(const std::string& _name)
          : name(_name),
            molecularDiffusivity(0),
            firstOrderDecayConstant(0) {}
        std::string name;
        double molecularDiffusivity;
        double firstOrderDecayConstant;
      };

      struct COMP
      {
          Array<TRACER>& getTracerArray() {return tracers;}
          const Array<TRACER>& getTracerArray() const {return tracers;}
          bool HasTracer(const std::string& name) {
              for (int i=0; i<tracers.size(); ++i) {
                  if (tracers[i].name==name) {
                      return true;
                  }
              }
              return false;
          }
          std::string& operator [] (int i) {return tracers[i].name;}
          int size() const {return tracers.size();}
          void push_back(const std::string& str)
              {
                tracers.push_back(TRACER(str));
              }
          
      protected:
          Array<TRACER> tracers;
      };
      
      struct PHASE
      {
          PHASE() {}
          PHASE(double density, double viscosity, double diffusivity)
              : density(density), viscosity(viscosity), diffusivity(diffusivity) {}

          typedef std::map<std::string,COMP> CompMap;

          const CompMap& getCompMap() const {return comp_map;}
          COMP& operator [] (const std::string& compName) {return comp_map[compName];}
          int size() const {return comp_map.size();}
          double Density() const {return density;}
          double Viscosity() const {return viscosity;}
          double Diffusivity() const {return diffusivity;}

      protected:
          CompMap comp_map;
          double density;
          double viscosity;
          double diffusivity;
      };
      

      struct SolidChem
      {
          SolidChem() {}
          struct MineralData
          {
              MineralData()
                  : volume_frac(0), specific_surface_area(0) {}
              double volume_frac, specific_surface_area;
              ParameterList BuildPL() const {
                  ParameterList pl;
                  pl.set("Volume_Fraction",volume_frac);
                  pl.set("Specific_Surface_Area",specific_surface_area);
                  return pl;
              }
          };
          struct SorptionSiteData
          {
              SorptionSiteData() 
                  : site_density(0) {}
              double site_density;
              ParameterList BuildPL() const {
                  ParameterList pl;
                  pl.set("Site_Density",site_density);
                  return pl;
              }
          };
          struct SorptionIsothermData
          {
              SorptionIsothermData() 
                  : Kd(0), Langmuir_b(-1.0), Freundlich_n(-1.0) {}
              bool IsFreundlich() const {return Freundlich_n>0;}
              bool IsLangmuir() const {return Langmuir_b>0;}
              double Kd, Langmuir_b, Freundlich_n;
              ParameterList BuildPL() const {
                  ParameterList pl;
                  pl.set("Kd",Kd);
                  if (IsFreundlich()) {
                      pl.set("Freundlich_n",Freundlich_n);
                  }
                  else if (IsLangmuir()) {
                      pl.set("Langmuir_b",Langmuir_b);
                  }
                  return pl;
              }
          };
          void InitializeMinerals(const Array<std::string>& names) {
              minerals.clear();
              for (int i=0; i<names.size(); ++i) {
                  minerals[names[i]] = MineralData();
              }
          }
          SolidChem::MineralData& Mineral(const std::string& name) {return minerals[name];}
          const std::map<std::string,MineralData>& Minerals() const {return minerals;}
          bool IsAMineral(const std::string& name) const {
              return !(minerals.find(name)==minerals.end());
          }

          void InitializeSorptionSites(const Array<std::string>& names) {
              sorption_sites.clear();
              for (int i=0; i<names.size(); ++i) {
                  sorption_sites[names[i]] = SorptionSiteData();
              }
          }
          SolidChem::SorptionSiteData& SorptionSite(const std::string& name) {return sorption_sites[name];}
          const std::map<std::string,SorptionSiteData>& SorptionSites() const {return sorption_sites;}
          bool IsASorptionSite(const std::string& name) const {
              return !(sorption_sites.find(name)==sorption_sites.end());
          }

	  void InitializeSorptionIsotherm(const Array<std::string>& names) {
              sorption_isotherms.clear();
              for (int i=0; i<names.size(); ++i) {
                  sorption_isotherms[names[i]] = SorptionIsothermData();
              }
          }
          SolidChem::SorptionIsothermData& SorptionIsotherm(const std::string& solute_name) {return sorption_isotherms[solute_name];}
          const std::map<std::string,SorptionIsothermData>& SorptionIsotherm() const {return sorption_isotherms;}          
          bool HasSorptionIsotherm(const std::string& solute)  const
	      {
		return sorption_isotherms.find(solute) != sorption_isotherms.end();
              }
	bool InsertSorptionIsotherm(const std::string& soluteName, const SorptionIsothermData& iso) {
              if (HasSorptionIsotherm(soluteName)) {
                  return false;
              }
              else {
                  sorption_isotherms[soluteName] = iso;
              }
              return true;
          }
	
      protected:
          std::map<std::string,MineralData> minerals;
          std::map<std::string,SorptionSiteData> sorption_sites;
          std::map<std::string,SorptionIsothermData> sorption_isotherms; 
      };


      struct Solid
      {
          Solid() : has_cation_exchange(false) {}

          bool IsAMineral(const std::string& name) const {
              for (int i=0; i<mineral_names.size(); ++i) {
                  if (name==mineral_names[i]) return true;
              }
              return false;
          }
          bool IsASorptionSite(const std::string& name) const {
              for (int i=0; i<sorption_site_names.size(); ++i) {
                  if (name==sorption_site_names[i]) return true;
              }
              return false;
          }
          bool HasSorptionIsotherm(const std::string& name) const {
              for (int i=0; i<sorption_isotherm_names.size(); ++i) {
                  if (name==sorption_isotherm_names[i]) return true;
              }
              return false;
          }
          bool HasStuff() const {
              return mineral_names.size() || sorption_site_names.size() || has_cation_exchange || sorption_isotherm_names.size();
          }

          bool UsingSorption() const {
              return sorption_site_names.size() || has_cation_exchange || sorption_isotherm_names.size();

          }

          Array<std::string> mineral_names;
          Array<std::string> sorption_site_names;
          Array<std::string> sorption_isotherm_names;
          bool has_cation_exchange;
      };


      typedef std::map<std::string,StateFunc> StateFuncMap;
      struct StateDef
      {
          StateDef(const ParameterList& parameter_list);

          typedef PHASE::CompMap CompMap;

          typedef std::map<std::string,PHASE> Phases; // phase name to PHASE
          typedef std::map<std::string,PHASE::CompMap> PhaseCompMap; // phase name to CompMap

          Phases& getPhases() {return phases;}
          Solid& getSolid() {return solid;}
          const Solid& getSolid() const {return solid;}
          PhaseCompMap& getPhaseCompMap() {return phase_comp_map;}
          CompMap& operator [] (const std::string& phase) {return phase_comp_map[phase];}
          StateFuncMap& IC() {return state_ics;}
          StateFuncMap& BC() {return state_bcs;}
          const ParameterList& PL() const {return parameter_list;}
          void clear();
          bool HasSolidChem() const {return getSolid().HasStuff();}

      protected:
          void build_state_def();
          void build_solute_funcs(StateFuncMap& s,
                                  const std::string& top_level_label,
                                  const std::string& solute_section_label);


          const ParameterList& parameter_list;
          Phases phases;
          PhaseCompMap phase_comp_map;
          StateFuncMap state_ics, state_bcs;
          Solid solid;
          bool has_solid_chem;
      };

      class PLoptions
      {
      public:
          PLoptions(const ParameterList& p,
                    const Array<std::string>& reqLists,
                    const Array<std::string>& reqParms,
                    bool reqListsOnly,
                    bool reqParmsOnly)
              : reqLists(reqLists),
                reqParms(reqParms),
                reqListsOnly(reqListsOnly),
                reqParmsOnly(reqParmsOnly)
              {
                  init(p);
              }

          const Array<std::string>& OptLists() {return optLists;}
          const Array<std::string>& ReqLists() {return reqLists;}
          const Array<std::string>& OptParms() {return optParms;}
          const Array<std::string>& ReqParms() {return reqParms;}

      protected:
          void init(const ParameterList& p)
              {
                  optLists.clear();
                  optParms.clear();
                  for (ParameterList::ConstIterator i=p.begin(); i!=p.end(); ++i) {
                      const std::string& name = p.name(i);
                      if (p.isSublist(name))
                      {
                          bool found = false;
                          for (int j=0; j<reqLists.size() && !found; ++j) {
                              if (name==reqLists[j]) {
                                  found = true;
                              }
                          }

                          if (!found) {
                              optLists.push_back(name);
                          }
                      }
                      else {
                          bool found = false;
                          for (int j=0; j<reqParms.size() && !found; ++j) {
                              if (name==reqParms[j]) {
                                  found = true;
                              }
                          }
                          
                          if (!found) {
                              optParms.push_back(name);
                          }
                      }
                  }
                  if (reqListsOnly && optLists.size()) {
                      bailLists(p);
                  }
                  if (reqParmsOnly && optParms.size()) {
                      bailParms(p);
                  }
              }   

          void bailParms(const ParameterList& p) 
              {
		  std::cerr << "Unsupported option parameters(s) in ParameterList: " << p.name() << std::endl;
                  for (int i=0; i<optParms.size(); ++i) {
                      std::cerr << "\"" << optParms[i] << "\" ";
                  }
                  std::cerr << '\n';
                  std::cerr << "Required parameters: ";
                  if (reqParms.size()==0) {
                      std::cerr << "[None]";
                  } else {
                      for (int i=0; i<reqParms.size(); ++i) {
                          std::cerr << "\"" << reqParms[i] << "\" ";
                      }
                  }
                  std::cerr << '\n';
                  throw std::exception();
              }
          void bailLists(const ParameterList& p) 
              {
		  std::cerr << "Unsupported option list(s) in ParameterList: " << p.name() << std::endl;
                  for (int i=0; i<optLists.size(); ++i) {
                      std::cerr << "\"" << optLists[i] << "\" ";
                  }
                  std::cerr << '\n';
                  std::cerr << "Required list(s): ";
                  if (reqLists.size()==0) {
                      std::cerr << "[None]";
                  } else {
                      for (int i=0; i<reqLists.size(); ++i) {
                          std::cerr << "\"" << reqLists[i] << "\" ";
                      }
                  }
                  std::cerr << '\n';
                  throw std::exception();
              }

          Array<std::string> reqLists, optLists, reqParms, optParms;
          bool reqListsOnly, reqParmsOnly;
      };
  }
}
#endif // Input_Parser_Structured_H_
