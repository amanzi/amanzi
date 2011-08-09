#include "Teuchos_ParameterList.hpp"

#include "Region.hh"
#include "RectangularRegion.hh"
#include "MeshDefs.hh"
#include "Mesh.hh"

#include "Teuchos_MPISession.hpp"

using Amanzi::AmanziGeometry::Region;
using Amanzi::AmanziGeometry::RegionPtr;
using Amanzi::AmanziGeometry::RectangularRegion;
using Amanzi::AmanziGeometry::Point;
using Amanzi::AmanziMesh::Mesh;

using Teuchos::Array;
using Teuchos::ParameterList;
using Teuchos::ParameterEntry;

namespace Amanzi {
  namespace AmanziInput {
    struct Rock
    {
      std::string rock_label;
      std::string porosity_model;
      Array<double> porosity_params;
      std::string permeability_model;
      Array<double> permeability_params;
      std::string retention_model;
      Array<double> retention_params;
      Array<double> component_concentration;
      Array<std::string> regions;
    };
    
    struct BoundaryCondition
    {
      double value;
      std::string functional;
    };
    
    struct InitialCondition
    {
      std::string dir;
      Array<double> value;
      std::string functional;
    };
    
    struct Component
    {
      std::string name;
      std::map<std::string,InitialCondition> ic;
      std::vector<std::string> tracers;
      std::string phase;
      double mass_density;
      double viscosity;
      double diffusivity;
    };
    
    struct Layer
    {
      Layer(const Point& p0_,
            const Point& p1_,
            const std::string& label_)
        : p0(p0_), p1(p1_), label(label_)  {}
      Point p0, p1;
      std::string label;
    };
    
    struct LayerLT
    {
      bool operator() (const Layer& lhs, const Layer& rhs) {
        return lhs.p0[2] < rhs.p0[2];
      }
    };
    
    typedef std::map<std::string,RegionPtr> RegionMap;
    typedef std::list<Rock> RockList;
    typedef std::map<std::string,Component> ComponentMap;
    typedef std::map<std::string,ComponentMap> PhaseMap;
    typedef std::map<std::string,std::string> StringStringMap;
    
    const Rock&
    find_rock_for_region(const RockList& rock_list,
                         const std::string& region_name);
    
    RockList
    build_rock_list(const Teuchos::ParameterList& parameter_list);
    
    
    void
    define_regions_rock_new_mesh(const Teuchos::ParameterList& parameter_list,
                                 RegionMap&                    region_map,
                                 RockList&                     rock_list,
                                 ParameterList&                new_generate_sublist,
                                 StringStringMap&              user_to_mesh_label_map);
      
    ParameterList translate_state_sublist(const ParameterList& orig_params);
  }
}

      
