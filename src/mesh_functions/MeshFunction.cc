/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! <MISSING_ONELINE_DOCSTRING>


#include "MeshFunction.hh"

namespace Amanzi {
namespace Functions {


//
// process a list for regions, components, and functions
//
Spec
processSpecWithFunction(const Teuchos::ParameterList& list,
                        std::string function_name);
{
  Teuchos::Array<std::string> regions;
  if (list.isParameter("region")) {
    regions.push_back(list.get<std::string>("region"));
  } else {
    regions = list.get<Teuchos::Array<std::string> >("regions");
  }

  Teuchos::Array<std::string> components;
  if (list.isParameter("component")) {
    components.push_back(list.get<std::string>("component"));
  } else {
    Teuchos::Array<std::string> def{"cell", "face", "node"};
    components = list.get<Teuchos::Array<std::string> >("components", def);
  }

  Teuchos::ParameterList f_list = list.sublist(function_name, true);

  // Make the function.
  Teuchos::RCP<Function> f;
  FunctionFactory f_fact;
  try {
    f = Teuchos::rcp(f_fact.Create(f_list));
  } catch (std::exception& msg) {
    Errors::Message m;
    m << "error in sublist " << function_name << ": " << msg.what();
    throw(m);
  }

  // TODO: Currently this does not support multiple-DoF functions.   --etc
  Teuchos::RCP<MultiFunction> func = Teuchos::rcp(new MultiFunction(f));

  return std::make_tuple(std::move(regions), std::move(components), func);
}

//
// parse a list of region-based specs
std::pair<MultiPatchSpace,
          std::vector<Teuchos::RCP<const MultiFunction>>>
processListWithFunction(const Teuchos::ParameterList& list,
                        std::string function_name)
{
  MultiPatchSpace space;
  std::vector<Teuchos::RCP<MultiFunction>> functions;

  // All are expected to be sublists of identical structure.
  for (auto sublist : list) {
    std::string name = sublist.first;
    if (list.isSublist(name)) {
      Teuchos::ParameterList spec_plist = list.sublist(name);


      Spec spec;
      try {
        spec = processSpecWithFunction(spec_plist, function_name);
      } catch (Errors::Message& msg) {
        Errors::Message m;
        m << "in sublist " << name << ": " << msg.what();
        throw(m);
      }

      for (const auto& region : std::get<0>(spec)) {
        for (const auto& comp : std::get<1>(spec)) {
          try {
            auto entity_kind = AmanziMesh::entity_kind(comp);
          } catch (Errors::Message& msg) {
            Errors::Message m;
            m << "component \"" << comp << "\" provided in sublist \""
              << function_name << "\" must be one of \"cell\", \"face\", or \"node\".";
            throw(m);
          }
          space.AddPatch(region, entity_kind, 1);
          functions.emplace_back(std::get<2>(spec.second));
        }
      }
    } else { // ERROR -- parameter is not a sublist
      Errors::Message m;
      m << "parameter " << name << " is not a sublist";
      throw(m);
    }
  }
  return std::make_pair(std::move(space), std::move(functions));
}



// void ProcessListWithoutFunction_(const Teuchos::ParameterList&,
//         const Teuchos::RCP<Functions::BoundaryFunction>&) const;



} // namespace Functions
} // namesapce Amanzi
