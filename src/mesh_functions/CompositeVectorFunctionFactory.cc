/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
  Mesh Functions

  Factory for a CV function on a mesh.
*/

#include "errors.hh"

#include "MeshFunction.hh"
#include "MultiFunction.hh"

#include "CompositeVectorFunctionFactory.hh"

namespace Amanzi {
namespace Functions {

Teuchos::RCP<CompositeVectorFunction>
CreateCompositeVectorFunction(Teuchos::ParameterList& plist,
        CompositeVectorSpace& sample,
        std::vector<std::string>& componentname_list,
        bool dot_with_normal,
        const std::string& spatial_dist_method)
{
  Teuchos::RCP<MeshFunction> mesh_func = Teuchos::rcp(new MeshFunction(sample.Mesh()));
  componentname_list.clear();

  // top level plist contains sublists containing the entry
  for (const auto& lcv : plist) {
    std::string name = lcv.first;

    if (plist.isSublist(name)) {
      Teuchos::ParameterList& sublist = plist.sublist(name);

      // grab regions from the sublist
      std::vector<std::string> regions;
      if (sublist.isParameter("region")) {
        if (sublist.isType<std::string>("region")) {
          regions.emplace_back(sublist.get<std::string>("region"));
        } else {
          Errors::Message msg;
          msg << "CompositeVectorFunctionFactory \"" << plist.name() << "(" << name << ")"
              << "\": parameter \"region\" should be a string.";
          Exceptions::amanzi_throw(msg);
        }
      } else if (sublist.isParameter("regions")) {
        if (sublist.isType<Teuchos::Array<std::string>>("regions")) {
          regions = sublist.get<Teuchos::Array<std::string>>("regions").toVector();
        } else {
          Errors::Message msg;
          msg << "CompositeVectorFunctionFactory \"" << plist.name() << "(" << name << ")"
              << "\": parameter \"regions\" should be an Array(string).";
          Exceptions::amanzi_throw(msg);
        }
      } else {
        Errors::Message msg;
        msg << "CompositeVectorFunctionFactory \"" << plist.name() << "(" << name << ")"
            << "\": parameter \"region\" or \"regions\" must exist.";
        Exceptions::amanzi_throw(msg);
      }

      // grab the name of the components from the list
      std::vector<std::string> components;
      if (sublist.isParameter("component")) {
        if (sublist.isType<std::string>("component")) {
          components.emplace_back(sublist.get<std::string>("component"));
        } else {
          Errors::Message msg;
          msg << "CompositeVectorFunctionFactory \"" << plist.name() << "(" << name << ")"
              << "\": parameter \"component\" should be a string.";
          Exceptions::amanzi_throw(msg);
        }
      } else if (sublist.isParameter("components")) {
        if (sublist.isType<Teuchos::Array<std::string>>("components")) {
          components = sublist.get<Teuchos::Array<std::string>>("components").toVector();
        } else {
          Errors::Message msg;
          msg << "CompositeVectorFunctionFactory \"" << plist.name() << "(" << name << ")"
              << "\": parameter \"components\" should be an Array(string).";
          Exceptions::amanzi_throw(msg);
        }
      } else {
        Errors::Message msg;
        msg << "CompositeVectorFunctionFactory \"" << plist.name() << "(" << name << ")"
            << "\": parameter \"component\" or \"components\" must exist.";
        Exceptions::amanzi_throw(msg);
      }

      // parse special case: initialize all existing components
      if (components.size() == 1 && components[0] == "*") {
        components.clear();
        for (auto it = sample.begin(); it != sample.end(); ++it) components.emplace_back(*it);
      }

      // grab the name of the locations from the list
      std::vector<AmanziMesh::Entity_kind> locations;
      if (sublist.isParameter("location")) {
        if (sublist.isType<std::string>("location")) {
          locations.emplace_back(AmanziMesh::createEntityKind(sublist.get<std::string>("location")));
        } else {
          Errors::Message msg;
          msg << "CompositeVectorFunctionFactory \"" << plist.name() << "(" << name << ")"
              << "\": parameter \"location\" should be a string.";
          Exceptions::amanzi_throw(msg);
        }
      } else if (sublist.isParameter("locations")) {
        if (sublist.isType<Teuchos::Array<std::string>>("locations")) {
          auto location_strings = sublist.get<Teuchos::Array<std::string>>("locations").toVector();
          for (const auto& loc : location_strings) {
            locations.emplace_back(AmanziMesh::createEntityKind(loc));
          }
        } else {
          Errors::Message msg;
          msg << "CompositeVectorFunctionFactory \"" << plist.name() << "(" << name << ")"
              << "\": parameter \"locations\" should be an Array(string).";
          Exceptions::amanzi_throw(msg);
        }
      } else {
        for (const auto& comp : components) {
          locations.emplace_back(AmanziMesh::createEntityKind(comp));
        }
      }

      // create the function
      Teuchos::RCP<MultiFunction> func;
      if (sublist.isSublist("function")) {
        Teuchos::ParameterList& func_plist = sublist.sublist("function");
        func = Teuchos::rcp(new MultiFunction(func_plist));
      } else {
        Errors::Message msg;
        msg << "CompositeVectorFunctionFactory \"" << plist.name() << "(" << name << ")"
            << "\": missing \"function\" sublist.";
        Exceptions::amanzi_throw(msg);
      }

      // From the above data, add to the cv function.
      // Loop through components, adding a spec/component name for each.
      for (int i = 0; i != components.size(); ++i) {
        // get the entity kind based upon the sample vector
        const std::string& component = components[i];
        AmanziMesh::Entity_kind location = locations[i];

        // this call will either Add the component if it doesn't exist, or
        // confirm the structure if it does.
        if (dot_with_normal) {
          sample.AddComponent(component, location, 1);
        } else {
          sample.AddComponent(component, location, func->size());
        }

        // -- Create the domain,
        auto domain = Teuchos::rcp(new MeshFunction::Domain(regions, location));

        // -- and the spec,
        auto spec = Teuchos::rcp(new MeshFunction::Spec(domain, func));

        mesh_func->AddSpec(spec);
        componentname_list.emplace_back(component);
      }

    } else {
      Errors::Message msg;
      msg << "CompositeVectorFunctionFactory \"" << plist.name() << "(" << name << ")"
          << "\": is not a sublist.";
      Exceptions::amanzi_throw(msg);
    }
  }

  // create the function
  return Teuchos::rcp(new CompositeVectorFunction(mesh_func, componentname_list, dot_with_normal, spatial_dist_method));
};



Teuchos::RCP<CompositeVectorFunction>
CreateCompositeVectorFunction(Teuchos::ParameterList& plist,
        const CompositeVectorSpace& sample,
        std::vector<std::string>& componentname_list,
        bool dot_with_normal,
        const std::string& spatial_dist_method)
{
  Teuchos::RCP<MeshFunction> mesh_func = Teuchos::rcp(new MeshFunction(sample.Mesh()));
  componentname_list.clear();

  // top level plist contains sublists containing the entry
  for (const auto& lcv : plist) {
    std::string name = lcv.first;

    if (plist.isSublist(name)) {
      Teuchos::ParameterList& sublist = plist.sublist(name);

      // grab regions from the sublist
      std::vector<std::string> regions;
      if (sublist.isParameter("region")) {
        if (sublist.isType<std::string>("region")) {
          regions.emplace_back(sublist.get<std::string>("region"));
        } else {
          Errors::Message msg;
          msg << "CompositeVectorFunctionFactory \"" << plist.name() << "(" << name << ")"
              << "\": parameter \"region\" should be a string.";
          Exceptions::amanzi_throw(msg);
        }
      } else if (sublist.isParameter("regions")) {
        if (sublist.isType<Teuchos::Array<std::string>>("regions")) {
          regions = sublist.get<Teuchos::Array<std::string>>("regions").toVector();
        } else {
          Errors::Message msg;
          msg << "CompositeVectorFunctionFactory \"" << plist.name() << "(" << name << ")"
              << "\": parameter \"regions\" should be an Array(string).";
          Exceptions::amanzi_throw(msg);
        }
      } else {
        Errors::Message msg;
        msg << "CompositeVectorFunctionFactory \"" << plist.name() << "(" << name << ")"
            << "\": parameter \"region\" or \"regions\" must exist.";
        Exceptions::amanzi_throw(msg);
      }

      // grab the name of the components from the list
      std::vector<std::string> components;
      if (sublist.isParameter("component")) {
        if (sublist.isType<std::string>("component")) {
          components.emplace_back(sublist.get<std::string>("component"));
        } else {
          Errors::Message msg;
          msg << "CompositeVectorFunctionFactory \"" << plist.name() << "(" << name << ")"
              << "\": parameter \"component\" should be a string.";
          Exceptions::amanzi_throw(msg);
        }
      } else if (sublist.isParameter("components")) {
        if (sublist.isType<Teuchos::Array<std::string>>("components")) {
          components = sublist.get<Teuchos::Array<std::string>>("components").toVector();
        } else {
          Errors::Message msg;
          msg << "CompositeVectorFunctionFactory \"" << plist.name() << "(" << name << ")"
              << "\": parameter \"components\" should be an Array(string).";
          Exceptions::amanzi_throw(msg);
        }
      } else {
        Errors::Message msg;
        msg << "CompositeVectorFunctionFactory \"" << plist.name() << "(" << name << ")"
            << "\": parameter \"component\" or \"components\" must exist.";
        Exceptions::amanzi_throw(msg);
      }

      // parse special case: initialize all existing components
      if (components.size() == 1 && components[0] == "*") {
        components.clear();
        for (auto it = sample.begin(); it != sample.end(); ++it) components.emplace_back(*it);
      }

      // grab the name of the locations from the list
      std::vector<AmanziMesh::Entity_kind> locations;
      if (sublist.isParameter("location")) {
        if (sublist.isType<std::string>("location")) {
          locations.emplace_back(AmanziMesh::createEntityKind(sublist.get<std::string>("location")));
        } else {
          Errors::Message msg;
          msg << "CompositeVectorFunctionFactory \"" << plist.name() << "(" << name << ")"
              << "\": parameter \"location\" should be a string.";
          Exceptions::amanzi_throw(msg);
        }
      } else if (sublist.isParameter("locations")) {
        if (sublist.isType<Teuchos::Array<std::string>>("locations")) {
          auto location_strings = sublist.get<Teuchos::Array<std::string>>("locations").toVector();
          for (const auto& loc : location_strings) {
            locations.emplace_back(AmanziMesh::createEntityKind(loc));
          }
        } else {
          Errors::Message msg;
          msg << "CompositeVectorFunctionFactory \"" << plist.name() << "(" << name << ")"
              << "\": parameter \"locations\" should be an Array(string).";
          Exceptions::amanzi_throw(msg);
        }
      } else {
        for (const auto& comp : components) {
          locations.emplace_back(AmanziMesh::createEntityKind(comp));
        }
      }

      // create the function
      Teuchos::RCP<MultiFunction> func;
      if (sublist.isSublist("function")) {
        Teuchos::ParameterList& func_plist = sublist.sublist("function");
        func = Teuchos::rcp(new MultiFunction(func_plist));
      } else {
        Errors::Message msg;
        msg << "CompositeVectorFunctionFactory \"" << plist.name() << "(" << name << ")"
            << "\": missing \"function\" sublist.";
        Exceptions::amanzi_throw(msg);
      }

      // From the above data, COMPARE WITH the CV function
      for (int i = 0; i != components.size(); ++i) {
        // get the entity kind based upon the sample vector
        const std::string& component = components[i];
        if (!sample.HasComponent(component)) {
          Errors::Message msg;
          msg << "CompositeVectorFunctionFactory \"" << plist.name() << "(" << name << ")"
              << "\": specified component \"" << component
              << "\" is either not valid or this vector does not include such a component.";
          Exceptions::amanzi_throw(msg);
        }

        // check n_dofs matches
        if (!dot_with_normal) {
          int n_dofs = sample.NumVectors(component);
          if (n_dofs != func->size()) {
            Errors::Message msg;
            msg << "CompositeVectorFunctionFactory \"" << plist.name() << "(" << name << ")"
                << "\": component \"" << component << "\" provided a function with " << func->size()
                << " degrees of freedom, but the vector expects " << n_dofs << " dofs.";
            Exceptions::amanzi_throw(msg);
          }
        }

        // this call will either Add the component if it doesn't exist, or
        // confirm the structure if it does.
        AmanziMesh::Entity_kind location = sample.Location(component);

        // -- Create the domain,
        auto domain = Teuchos::rcp(new MeshFunction::Domain(regions, location));

        // -- and the spec,
        auto spec = Teuchos::rcp(new MeshFunction::Spec(domain, func));

        mesh_func->AddSpec(spec);
        componentname_list.emplace_back(component);
      }

    } else {
      Errors::Message msg;
      msg << "CompositeVectorFunctionFactory \"" << plist.name() << "(" << name << ")"
          << "\": is not a sublist.";
      Exceptions::amanzi_throw(msg);
    }
  }

  // create the function
  return Teuchos::rcp(new CompositeVectorFunction(mesh_func, componentname_list, dot_with_normal, spatial_dist_method));
};


} // namespace Functions
} // namespace Amanzi
