#include "InputParser.H"

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_MPISession.hpp"
#include "Teuchos_StrUtils.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

namespace Amanzi {
  namespace AmanziInput {

    static Rock INVALID_ROCK;
    
    const Rock&
    find_rock_for_region(const RockList& rock_list,
                         const std::string& region_name)
    {
      for (RockList::const_iterator it=rock_list.begin(); it!=rock_list.end(); ++it) {
        const Rock& rock = *it;
        for (int i=0; i<rock.regions.size(); ++i) {
          if (rock.regions[i] == region_name) {
            return rock;
          }
        }
      }
      return INVALID_ROCK;
    }
    
    
    
    RockList
    build_rock_list(const Teuchos::ParameterList& parameter_list)
    {
      const ParameterList& rock_sublist = parameter_list.sublist("Rock");
      
      RockList rock_list;
      for (ParameterList::ConstIterator i=rock_sublist.begin(); i!=rock_sublist.end(); ++i) {
        
        Rock rock;
        rock.rock_label = rock_sublist.name(i);
        
        const Teuchos::ParameterEntry& rock_entry = rock_sublist.getEntry(rock.rock_label);
        
        if ( !rock_entry.isList() ) {
          if (Teuchos::MPISession::getRank() == 0) {
            std::cerr << "Rock section must define only rocks.  " << rock.rock_label << " is not a rock" << std::endl;
          }
          throw std::exception();
        }
        
        const ParameterList& rock_def = rock_sublist.sublist(rock.rock_label);
        
        for (ParameterList::ConstIterator j=rock_def.begin(); j!=rock_def.end(); ++j) {
          
          const std::string& rlabel = rock_def.name(j);
          const Teuchos::ParameterEntry& rentry = rock_def.getEntry(rlabel);
          
          if (rentry.isList()) {
            
            if (rlabel=="porosity: uniform") {
              
              rock.porosity_model = rlabel;
              const ParameterList& porosity_model_parameter_list = rock_def.sublist(rock.porosity_model);
              rock.porosity_params.resize(1);
              rock.porosity_params[0] = porosity_model_parameter_list.get<double>("porosity");
            }
            
            else if (rlabel=="permeability: uniform") {
              
              rock.permeability_model = rlabel;
              const ParameterList& permeability_model_parameter_list = rock_def.sublist(rock.permeability_model);
              rock.permeability_params.resize(1);
              rock.permeability_params[0] = permeability_model_parameter_list.get<double>("permeability");
            }
            
            if (rlabel=="water retention: vG") {
              
              rock.retention_model = rlabel;
              const ParameterList& retention_model_parameter_list = rock_def.sublist(rock.retention_model);
              rock.retention_params = retention_model_parameter_list.get<Array<double> >("m alpha slr");
            }
            
          } else {
            
            rock.regions = rock_def.get<Array<std::string> >("regions");
            
          }
        }
        
        rock_list.push_back(rock);
      }
      
      return rock_list;
    }

    void
    define_regions_rock_new_mesh(const Teuchos::ParameterList& parameter_list,
                                 RegionMap&                    region_map,
                                 RockList&                     rock_list,
                                 ParameterList&                new_generate_sublist,
                                 UserToMeshLabelMap&           user_to_mesh_labelID_map)
    {
      rock_list = build_rock_list(parameter_list);
      
      // Build set of all regions used by rock list
      std::set<std::string> region_set;
      for (RockList::const_iterator rit=rock_list.begin(); rit!=rock_list.end(); ++rit) {
        const Rock& rock = *rit;
        for (int i=0; i<rock.regions.size(); ++i) {
          region_set.insert(rock.regions[i]);
        }
      }
      
      const ParameterList& mesh_sublist = parameter_list.sublist("Mesh");
      const ParameterList& generate_sublist = mesh_sublist.sublist("Generate");
      
      std::set<Layer,LayerLT> layers;
      
      // As we build the regions, check for ones used in setting the rocks, and pull additional information
      const ParameterList& region_sublist = parameter_list.sublist("Regions");
      
      for (ParameterList::ConstIterator i=region_sublist.begin(); i!=region_sublist.end(); ++i) {
        
        const std::string& region_name = region_sublist.name(i);
        
        const Teuchos::ParameterEntry& region_entry = region_sublist.getEntry(region_name);
        
        if ( !region_entry.isList() ) {
          if (Teuchos::MPISession::getRank() == 0) {
            std::cerr << "Region section must define only regions.  " << region_name << " is not a region" << std::endl;
          }
          throw std::exception();
        }
        
        const ParameterList& region_def = region_sublist.sublist(region_name);
        ParameterList::ConstIterator region_params_it = region_def.begin();
        
        const std::string& region_type = region_def.name(region_params_it);
        const ParameterEntry& region_params_entry = region_def.getEntry(region_type);
        
        if ( !region_params_entry.isList() ) {
          if (Teuchos::MPISession::getRank() == 0) {
            std::cerr << "Region definition must be a list.  " << region_type << " is not a list" << std::endl;
          }
          throw std::exception();
        }
        
        if (region_type != "box") {
          if (Teuchos::MPISession::getRank() == 0) {
            std::cerr << "Region must be a box. " << region_type << " is not a supported" << std::endl;
          }
          throw std::exception();
        }
        
        const ParameterList& params_sublist = region_def.sublist(region_type);
        
        const Array<double>& lo = params_sublist.get<Array<double> >("lo");
        const Array<double>& hi = params_sublist.get<Array<double> >("hi");
        
        Point p0(lo[0],lo[1],lo[2]);
        Point p1(hi[0],hi[1],hi[2]);
        region_map[region_name] = RegionPtr(new RectangularRegion(p0,p1));
        
        layers.insert(Layer(p0,p1,region_name));
      }
      
      // Add default regions
      const Array<double> lo = generate_sublist.get<Array<double> >("Domain Low Corner");
      const Array<double> hi = generate_sublist.get<Array<double> >("Domain High Corner");
      Point p0(lo[0],lo[1],lo[2]);
      Point p1(hi[0],hi[1],hi[2]);
      region_map["all"] = RegionPtr(new RectangularRegion(p0,p1));
      
      
      // Make a new Mesh sublist
      ParameterList new_mesh_sublist = ParameterList(parameter_list.sublist("Mesh"));
      new_generate_sublist = new_mesh_sublist.sublist("Generate");
      const Array<int>& n = new_generate_sublist.get<Array<int> >("Number of Cells");
      new_generate_sublist.set<int>("Number of Cells in X",n[0]);
      new_generate_sublist.set<int>("Number of Cells in Y",n[1]);
      new_generate_sublist.set<int>("Number of Cells in Z",n[2]);
      new_generate_sublist.remove("Number of Cells");
      
      Array<double>& dlo = new_generate_sublist.get<Array<double> >("Domain Low Corner");
      new_generate_sublist.set<double>("X_Min",dlo[0]);
      new_generate_sublist.set<double>("Y_Min",dlo[1]);
      new_generate_sublist.set<double>("Z_Min",dlo[2]);
      new_generate_sublist.remove("Domain Low Corner");
      
      Array<double>& dhi = new_generate_sublist.get<Array<double> >("Domain High Corner");
      new_generate_sublist.set<double>("X_Max",dhi[0]);
      new_generate_sublist.set<double>("Y_Max",dhi[1]);
      new_generate_sublist.set<double>("Z_Max",dhi[2]);
      new_generate_sublist.remove("Domain High Corner");
      
      new_generate_sublist.set<int>("Number of mesh blocks",layers.size());
      int cnt = 1;
      for (std::set<Layer>::const_iterator it=layers.begin(); it!=layers.end(); ++it) {
        
        ParameterList sublist;
        sublist.set<double>("Z0", it->p0[2]);
        sublist.set<double>("Z1", it->p1[2]);
        
        // FIXME: This assumes that we know the scheme by which the mesher derives the labels 
        //   and can reproduce it here.
        std::stringstream mesh_label;
        mesh_label << "Mesh block " << cnt;
        new_generate_sublist.set(mesh_label.str(), sublist);
        
        user_to_mesh_labelID_map[it->label] = std::pair<std::string,int>(mesh_label.str(),cnt);
        cnt++;
      }
    }

    bool version_at_or_below(const std::string& input_version, const std::string version_id) {

      Array<std::string> version_tokens
        = Teuchos::StrUtils::stringTokenizer(Teuchos::StrUtils::varSubstitute(input_version,"."," "));
    
      Array<std::string> version_id_tokens
        = Teuchos::StrUtils::stringTokenizer(Teuchos::StrUtils::varSubstitute(version_id,"."," "));
    
      if (version_tokens.size() >= 1) {
          
        if (version_id_tokens.size() == 0) {
            
          // input_version is a number, id is not -> version is newer
          return false;
        }
        else {
            
          if (version_id_tokens.size() >= 1) {
              
            if (version_tokens[0] > version_id_tokens[0]) {
                
              // If first number is higher, is newer
              return false;
                
            } else {
                
              if (version_tokens[0] == version_id_tokens[0]) {
                  
                if (version_tokens.size() >= 2) {
                    
                  if (version_id_tokens.size() >= 2) {
                      
                    // Only allowed two digits, if first agrees, then use second to decide
                    return version_tokens[1] <= version_id_tokens[1];
                      
                  }
                    
                } else {
                    
                  // If first number agrees, and version has less digits, is older
                  return version_id_tokens.size() >= 2;
                }
              }
            }
          }
        }
      }
    }

    
    ParameterList translate_state_sublist(const ParameterList& orig_params)
    {
      ParameterList new_params(orig_params);

      std::string input_version = orig_params.get<std::string>("Amanzi input format version");

      // Bail if format is older than 0.1
      if (version_at_or_below(input_version,"0.1"))
          return new_params;

      // Build a rock list
      RegionMap region_map;
      RockList rock_list;
      ParameterList new_mesh_sublist;
      UserToMeshLabelMap user_to_mesh_label_map;

      define_regions_rock_new_mesh(orig_params, region_map, rock_list, new_mesh_sublist, user_to_mesh_label_map);

      ParameterList& mesh_sublist = new_params.sublist("Mesh");
      mesh_sublist.remove("Generate");
      mesh_sublist.set("Generate",new_mesh_sublist);

      PhaseMap phase_map;
      const ParameterList& state_sublist = orig_params.sublist("State");
      for (ParameterList::ConstIterator i=state_sublist.begin(); i!=state_sublist.end(); ++i) {

        const std::string& name = state_sublist.name(i);
        const ParameterEntry& entry = state_sublist.getEntry(name);

        // If this entry is a list, it is either a component spec or a BC
        if ( entry.isList() ) {

          // If list is named after a region, this is a bc list
          RegionMap::const_iterator it=region_map.find(name);
          if (it != region_map.end()) {
        
            const ParameterList& sublist = state_sublist.sublist(name);
            for (ParameterList::ConstIterator j=sublist.begin(); j!=sublist.end(); ++j) {
          
              const std::string& bc_entry_name = sublist.name(j);

              const ParameterEntry& bc_entry = sublist.getEntry(bc_entry_name);
              if (!bc_entry.isList()) {
                if (Teuchos::MPISession::getRank() == 0) {
                  std::cerr << "bc functional must be a list.  " << bc_entry_name << " is not a list" << std::endl;
                }
                throw std::exception();
              }

              if (name == "bc: inflow" 
                  || name == "bc: outflow"
                  || name == "bc: seepage"
                  || name == "bc: noflow") {

                BoundaryCondition bc;
                bc.functional = name;
                const ParameterList& bc_param_list = sublist.sublist(bc_entry_name);
                if (name == "bc: inflow") {
                  bc.value = bc_param_list.get<double>("value");
                } else if (name == "bc: outflow") {
                  bc.value = bc_param_list.get<double>("value");
                } else if (name == "bc: seepage") {
                  bc.value = bc_param_list.get<double>("water table height");
                } else if (name == "bc: noflow") {
                }

              } else {
                if (Teuchos::MPISession::getRank() == 0) {
                  std::cerr << "unknonwn bc functional: " << bc_entry_name << std::endl;
                }
                throw std::exception();
              } // recognized bc functional
            } // iter over bc list
          } // bc list

          // Otherwise, this list is a component
          else {
            Component comp;
            comp.name = name;
            const ParameterList& comp_list = state_sublist.sublist(comp.name);
        
            for (ParameterList::ConstIterator k=comp_list.begin(); k!=comp_list.end(); ++k) {
              const std::string& cname = comp_list.name(k);
          
              const ParameterEntry& centry = comp_list.getEntry(cname);
          
              // If this is a list, it is either an IC or a source.  If it is named after a region (or "default"), it is an IC, 
              // otherwise it defines a named source term
              if (centry.isList()) {
            
                RegionMap::const_iterator it=region_map.find(cname);
                if (it != region_map.end() || cname=="default") {
              
                  // If this is named after a region, it must be an IC
                  InitialCondition ic;
                  const ParameterList& ic_list_for_region = comp_list.sublist(cname);
              
                  for (ParameterList::ConstIterator L=ic_list_for_region.begin(); L!=ic_list_for_region.end(); ++L) {
                
                    const std::string& ic_list_element_name = ic_list_for_region.name(L);
                
                    const ParameterEntry& ic_entry = ic_list_for_region.getEntry(ic_list_element_name);
                
                    if (ic_entry.isList()) {
                  
                      if (ic_list_element_name=="ic: uniform" 
                          || ic_list_element_name=="ic: coordinate-aligned-linear" 
                          || ic_list_element_name=="ic: quadratic" 
                          || ic_list_element_name=="ic: exponential") {
                    
                        ic.functional = ic_list_element_name;
                        ic.value.resize(1);
                        const ParameterList& ic_param_list = ic_list_for_region.sublist(ic.functional);
                        if (ic.functional == "ic: uniform") {
                          ic.value[0] = ic_param_list.get<double>("value");
                        } else if (ic.functional == "ic: coordinate-aligned-linear") {
                          ic.dir = ic_param_list.get<std::string>("dir");
                          ic.value = ic_param_list.get<Array<double> >("x0 y0 slope");
                        } else if (ic.functional == "ic: quadratic") {
                          ic.dir = ic_param_list.get<std::string>("dir");
                          ic.value = ic_param_list.get<Array<double> >("x0 y0 slope curvature");
                        } else if (ic.functional == "bc: exponential") {
                          ic.dir = ic_param_list.get<std::string>("dir");
                          ic.value = ic_param_list.get<Array<double> >("x0 y0 offset growth");
                        }
                    
                        comp.ic[cname] = ic;
                    
                      } else {
                        // This is not an IC or a source, it is an error
                        if (Teuchos::MPISession::getRank() == 0) {
                          std::cerr << "Unknown list element " << ic_list_element_name << " in component list " << comp.name << std::endl;
                        }
                        throw std::exception();
                      }
                    } // is functional or source
                  } // iterator over ic list
                } // is an ic
                else {
              
                  // Otherwise, not an IC, so must be a source                    
                  const std::string& source_name = cname; 
                  if (cname=="source") {
                    if (Teuchos::MPISession::getRank() == 0) {
                      std::cerr << "Sources not yet supported: " << source_name <<  std::endl;
                    }
                    throw std::exception();
                  }
                } // an IC or a source
              } //is a list under comp
              else {

                // is not a list under comp
                if (cname == "mass density") {
                  comp.mass_density = comp_list.get<double>(cname);
                } else if (cname == "viscosity") {
                  comp.viscosity = comp_list.get<double>(cname);
                } else if (cname == "diffusivity") {
                  comp.diffusivity = comp_list.get<double>(cname);
                } else {
                  comp.phase = comp_list.get<std::string>(cname);
                }
              }
            } // list/no-list under comp

            // Insert this component into the phase map
            phase_map[comp.phase][comp.name] = comp;

          } // iter over comp
        }
      }
  
      std::string aqueous_name = "Aqueous";
      if (! phase_map.count(aqueous_name) ) {
        if (Teuchos::MPISession::getRank() == 0) {
          std::cerr << "Inputs translator must define at least one component in "
                    << aqueous_name << " phase " <<  std::endl;
        }
        throw std::exception();
      }
      const ComponentMap& aqueous = phase_map[aqueous_name];

      std::string water_name = "Water";
      if (!aqueous.count(water_name)) {
        if (Teuchos::MPISession::getRank() == 0) {
          std::cerr << "Inputs translator must define " << water_name << " component in " 
                    << aqueous_name << " phase " <<  std::endl;
        }
        throw std::exception();
      }
      const Component& water = phase_map[aqueous_name][water_name];


      int number_of_phases = phase_map.size();
      if (number_of_phases > 1) {
        if (Teuchos::MPISession::getRank() == 0) {
          std::cerr << "Inputs translater does not support more than one phase" <<  std::endl;
        }
        throw std::exception();
      }


      int number_of_aqueous_components = aqueous.size();
      if (number_of_aqueous_components > 1) {
        if (Teuchos::MPISession::getRank() == 0) {
          std::cerr << "Inputs translater does not support more than one aqueous component" <<  std::endl;
        }
        throw std::exception();
      }

      // Build new State parameter list

      // Write gravity array as 3 entries
      new_params.remove("State");
      const ParameterList& orig_state = orig_params.sublist("State");
      ParameterList new_state;
      const Array<double>& g = orig_state.get<Array<double> >("Gravity");
      new_state.set<double>("Gravity x",g[0]);
      new_state.set<double>("Gravity y",g[1]);
      new_state.set<double>("Gravity z",g[2]);

      new_state.set<int>("Number of mesh blocks",water.ic.size());
      new_state.set<int>("Number of component concentrations",number_of_aqueous_components);
      new_state.set<double>("Constant water saturation",water.ic.begin()->second.value[0]);
      new_state.set<double>("Constant water density",water.mass_density);
      new_state.set<double>("Constant viscosity",water.viscosity);

      ParameterList retention_list;
      for (std::map<std::string,InitialCondition>::const_iterator it=water.ic.begin(); it!=water.ic.end(); ++it) {
        ParameterList reg_block;
        std::stringstream reg_block_name;
        reg_block_name << user_to_mesh_label_map[it->first].first;
        reg_block.set<int>("Mesh block ID",user_to_mesh_label_map[it->first].second);


        const Rock& rock = find_rock_for_region(rock_list,it->first);


        // Retention model
        std::string this_retention_label = rock.retention_model;
        const Array<double>& this_retention_values = rock.retention_params;

        ParameterList this_new_retention_list;
        if (this_retention_label=="water retention: vG") {
          this_new_retention_list.set<std::string>("Water retention model", "van Genuchten");
          if (this_retention_values.size() != 3) {
            if (Teuchos::MPISession::getRank() == 0) {
              std::cerr << "Wrong number of vG parameters in water retention model for rock: " << rock.rock_label <<  std::endl;
            }
            throw std::exception();
          }
          this_new_retention_list.set<double>("van Genuchten m",this_retention_values[0]);
          this_new_retention_list.set<double>("van Genuchten alpha",this_retention_values[1]);
          this_new_retention_list.set<double>("van Genuchten residual saturation",this_retention_values[2]);
        }
        else {
          if (Teuchos::MPISession::getRank() == 0) {
            std::cerr << "Unsupported water retention model for rock: " << rock.rock_label <<  std::endl;
          }
          throw std::exception();
        }
        this_new_retention_list.set<int>("Region ID",user_to_mesh_label_map[it->first].second);
        std::stringstream retention_label;
        retention_label << "Water retention model " << user_to_mesh_label_map[it->first].second - 1; // FIXME: Numbering weird!
        retention_list.set(retention_label.str(), this_new_retention_list);

        // Porosity
        if (rock.porosity_model!="porosity: uniform") {
          if (Teuchos::MPISession::getRank() == 0) {
            std::cerr << "Inputs translater does not support nonuniform porosity" <<  std::endl;
          }
          throw std::exception();
        } else {
          reg_block.set<double>("Constant porosity",rock.porosity_params[0]);
        }

        if (rock.permeability_model!="permeability: uniform") {
          if (Teuchos::MPISession::getRank() == 0) {
            std::cerr << "Inputs translater does not support nonuniform permeability **" << rock.permeability_model << "*" <<  std::endl;
          }
          throw std::exception();
        } else {
          reg_block.set<double>("Constant permeability",rock.permeability_params[0]);
        }

        std::map<std::string,InitialCondition>::const_iterator wit = water.ic.find(it->first);
        if (wit != water.ic.end()) {

          std::stringstream label;; 
          label << "Constant component concentration " << phase_map[aqueous_name].size() - 1;
          reg_block.set<double>(label.str(),wit->second.value[0]);
        }
        else {
          if (Teuchos::MPISession::getRank() == 0) {
            std::cerr << "Inputs translater does not support nonuniform component initial condition" <<  std::endl;
          }
          throw std::exception();
        }

        new_state.set(reg_block_name.str(), reg_block);
      }

      new_params.sublist("Flow").sublist("Richards Problem").set("Water retention models",retention_list);
      new_params.set("State", new_state);

      new_params.remove("Regions");
      new_params.remove("Rock");


      return new_params;
    }
  }
}
