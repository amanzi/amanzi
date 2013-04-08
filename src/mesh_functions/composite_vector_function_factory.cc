/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Factory for a CV function on a mesh.
------------------------------------------------------------------------- */

#include "mesh_function.hh"
#include "vector_function_factory.hh"

#include "composite_vector_function_factory.hh"


namespace Amanzi {
namespace Functions {

Teuchos::RCP<CompositeVectorFunction>
CreateCompositeVectorFunction(Teuchos::ParameterList& plist,
        const CompositeVector& sample) {

  Teuchos::RCP<MeshFunction> mesh_func =
    Teuchos::rcp(new MeshFunction(sample.mesh()));
  std::vector<std::string> componentname_list;
  VectorFunctionFactory vfunc_factory;

  // top level plist contains sublists containing the entry
  for (Teuchos::ParameterList::ConstIterator lcv=plist.begin();
       lcv!=plist.end(); ++lcv) {
    std::string name = lcv->first;
    if (plist.isSublist(name)) {
      Teuchos::ParameterList& sublist = plist.sublist(name);

      // grab regions from the sublist
      std::vector<std::string> regions;
      if (sublist.isParameter("region")) {
        if (sublist.isType<std::string>("region")) {
          regions.push_back(sublist.get<std::string>("region"));
        } else {
          // ERROR -- invalid Region parameter
          ASSERT(0);
        }
      } else if (sublist.isParameter("regions")) {
        if (sublist.isType<Teuchos::Array<std::string> >("regions")) {
          regions = sublist.get<Teuchos::Array<std::string> >("regions").toVector();
        } else {
          // ERROR -- invalid Regions parameter
          ASSERT(0);
        }
      } else {
        // ERROR -- missing Region/Regions parameter
        ASSERT(0);
      }

      // grab the name of the components from the list
      std::vector<std::string> components;
      if (sublist.isParameter("component")) {
        if (sublist.isType<std::string>("component")) {
          components.push_back(sublist.get<std::string>("component"));
        } else {
          // ERROR -- invalid component parameter
          ASSERT(0);
        }
      } else if (sublist.isParameter("components")) {
        if (sublist.isType<Teuchos::Array<std::string> >("components")) {
          components = sublist.get<Teuchos::Array<std::string> >("components").toVector();
        } else {
          // ERROR -- invalid Regions parameter
          ASSERT(0);
        }
      } else {
        // ERROR -- missing component parameter
        ASSERT(0);
      }

      // get the function
      Teuchos::RCP<VectorFunction> func;
      if (sublist.isSublist("function")) {
        Teuchos::ParameterList& func_plist = sublist.sublist("function");
        func = vfunc_factory.Create(func_plist);
      } else {
        // ERROR -- missing function plist
        ASSERT(0);
      }

      // From the above data, add to the cv function.
      // Loop through components, adding a spec/component name for each.
      for (std::vector<std::string>::const_iterator component=components.begin();
           component!=components.end(); ++component) {

        // get the entity kind based upon the sample vector
        if (!sample.has_component(*component)) {
          // ERROR -- invalid component name
          ASSERT(0);
        }
        AmanziMesh::Entity_kind kind = sample.location(*component);

        // -- Create the domain,
        Teuchos::RCP<MeshFunction::Domain> domain =
          Teuchos::rcp(new MeshFunction::Domain(regions, kind));

        // -- and the spec,
        Teuchos::RCP<MeshFunction::Spec> spec =
          Teuchos::rcp(new MeshFunction::Spec(domain, func));

        mesh_func->AddSpec(spec);
        componentname_list.push_back(*component);
      }
    } else {
      // ERROR -- not a sublist
      ASSERT(0);
    }
  }

  // create the function
  return Teuchos::rcp(new CompositeVectorFunction(mesh_func,
          componentname_list));
};


Teuchos::RCP<CompositeVectorFunction>
CreateCompositeVectorFunction(Teuchos::ParameterList& plist,
        const CompositeVectorFactory& factory) {
  Teuchos::RCP<CompositeVector> sample = factory.CreateVector(true);
  return CreateCompositeVectorFunction(plist, *sample);
};

} // namespace
} // namespace


