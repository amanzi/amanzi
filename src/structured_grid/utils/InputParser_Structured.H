#ifndef Input_Parser_Structured_H_
#define Input_Parser_Structured_H_
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_MPISession.hpp"

namespace Amanzi {
  namespace AmanziInput {
    
    using Teuchos::ParameterList;
    using Teuchos::Array;

    static int ndim;
    static double simulation_time;
    static Array<double> domlo,domhi;      
    static double geometry_eps;

    static Array<std::string> user_derive_list;

    typedef struct {
	std::string phase;
	std::string component;
	std::string solute;
      } tridata;

    static std::map<std::string,std::string> AMR_to_Amanzi_label_map;
    std::string underscore(const std::string& instring);

    ParameterList convert_to_structured(const ParameterList& parameter_list);

    ParameterList setup_structured();
    
    void convert_to_structured_mesh    (const ParameterList& parameter_list, 
					ParameterList&       struc_list);
    void convert_to_structured_control (const ParameterList& parameter_list, 
					ParameterList&       struc_list);
    void convert_to_structured_region  (const ParameterList& parameter_list, 
					ParameterList&       struc_list);
    void convert_to_structured_material(const ParameterList& parameter_list, 
					ParameterList&       struc_list);
    void convert_to_structured_output  (const ParameterList& parameter_list, 
					ParameterList&       struc_list);

      struct ICBCFunc
      {
          ICBCFunc() {};
          ICBCFunc(const ParameterList& plist,
                   const std::string& Amanzi_type,
                   const std::string& units,
                   const std::string& label)
              : plist(plist), Amanzi_type(Amanzi_type), units(units), label(label) {}
          
          const ParameterList& PList() const {return plist;}
          const std::string& Amanzi_Type() const {return Amanzi_type;}
          const std::string& Label() const {return label;}
          const std::string& Units() const {return units;}
          
      protected:
          std::string Amanzi_type, label, units;
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
              : label(label), Amanzi_type(Amanzi_type), func_plist(func_plist), regions(regions) {}


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
      


      struct COMP
      {
          const Array<std::string>& getTracerArray() const {return tracers;}
          std::string& operator [] (int i) {return tracers[i];}
          int size() const {return tracers.size();}
          void push_back(const std::string& str)
              {
                  tracers.push_back(str);
              }

      protected:
          Array<std::string> tracers;
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
      



      typedef std::map<std::string,StateFunc> StateFuncMap;
      struct StateDef
      {
          StateDef(const ParameterList& parameter_list);

          typedef PHASE::CompMap CompMap;

          typedef std::map<std::string,PHASE> Phases; // phase name to PHASE
          typedef std::map<std::string,PHASE::CompMap> PhaseCompMap; // phase name to CompMap

          Phases& getPhases() {return phases;}
          PhaseCompMap& getPhaseCompMap() {return phase_comp_map;}
          CompMap& operator [] (const std::string& phase) {return phase_comp_map[phase];}
          StateFuncMap& IC() {return state_ics;}
          StateFuncMap& BC() {return state_bcs;}
          const ParameterList& PL() const {return parameter_list;}
          void clear();

      protected:
          void build_state_def();
          void build_solute_funcs(StateFuncMap& s,
                                  const std::string& top_level_label,
                                  const std::string& solute_section_label);


          const ParameterList& parameter_list;
          Phases phases;
          PhaseCompMap phase_comp_map;
          StateFuncMap state_ics, state_bcs;          
      };

      class PLoptions
      {
      public:
          PLoptions(const ParameterList& p,
                    const Array<std::string>& reqLists,
                    const Array<std::string>& reqParms,
                    bool reqListsOnly,
                    bool reqParmsOnly)
              : reqListsOnly(reqListsOnly), reqParmsOnly(reqParmsOnly),
                reqLists(reqLists), reqParms(reqParms)
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
                      bailLists();
                  }
                  if (reqParmsOnly && optParms.size()) {
                      bailParms();
                  }
              }   

          void bailParms() 
              {
                  std::cerr << "Unsupported option parameters(s) list in ParameterList: " << std::endl;
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
          void bailLists() 
              {
                  std::cerr << "Unsupported option list(s) list in ParameterList: " << std::endl;
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
