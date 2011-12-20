/*
This is the flow component of the Amanzi code. 
License: BSD
Authors: Neil Carlson (versions 1 & 2)  (nnc@lanl.gov)
*/

#include "boundary-function.hh"
#include "function-factory.hh"
#include "errors.hh"

#include "Flow_BC_Factory.hpp"

namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* TBA
****************************************************************** */
BoundaryFunction* FlowBCFactory::CreatePressure() const
{
  BoundaryFunction *bc = new BoundaryFunction(mesh_);
  try {
    process_pressure_list(params_->sublist("pressure"), bc);
  } catch (Errors::Message &msg) {
    Errors::Message m;
    m << "FlowBCFactory: \"pressure\" sublist error: " << msg.what();
    Exceptions::amanzi_throw(m);
  } catch (Teuchos::Exceptions::InvalidParameterType &msg) {
    Errors::Message m;
    m << "FlowBCFactory: \"pressure\" sublist error: not a sublist: " << msg.what();
    Exceptions::amanzi_throw(m);
  }
  return bc;
}


/* ******************************************************************
* TBA
****************************************************************** */
BoundaryFunction* FlowBCFactory::CreateMassFlux() const
{
  BoundaryFunction *bc = new BoundaryFunction(mesh_);
  try {
    process_mass_flux_list(params_->sublist("mass flux"), bc);
  } catch (Errors::Message &msg) {
    Errors::Message m;
    m << "FlowBCFactory: \"mass flux\" sublist error: " << msg.what();
    Exceptions::amanzi_throw(m);
  } catch (Teuchos::Exceptions::InvalidParameterType &msg) {
    Errors::Message m;
    m << "FlowBCFactory: \"mass flux\" sublist error: not a sublist: " << msg.what();
    Exceptions::amanzi_throw(m);
  }
  return bc;
}


/* ******************************************************************
* TBA
****************************************************************** */
BoundaryFunction* FlowBCFactory::CreateStaticHead(double p0, double density, double gravity) const
{
  BoundaryFunction *bc = new BoundaryFunction(mesh_);
  try {
    process_static_head_list(p0, density, gravity, params_->sublist("static head"), bc);
  } catch (Errors::Message &msg) {
    Errors::Message m;
    m << "FlowBCFactory: \"static head\" sublist error: " << msg.what();
    Exceptions::amanzi_throw(m);
  } catch (Teuchos::Exceptions::InvalidParameterType &msg) {
    Errors::Message m;
    m << "FlowBCFactory: \"static head\" sublist error: not a sublist: " << msg.what();
    Exceptions::amanzi_throw(m);
  }
  return bc;
}


/* ******************************************************************
* TBA
****************************************************************** */
void FlowBCFactory::process_pressure_list(Teuchos::ParameterList &list,
                                           BoundaryFunction *bc) const
{
  // Iterate through the BC specification sublists in the list.
  // All are expected to be sublists of identical structure.
  for (Teuchos::ParameterList::ConstIterator i = list.begin(); i != list.end(); ++i) {
    std::string name = i->first;
    if (list.isSublist(name)) {
      Teuchos::ParameterList &spec = list.sublist(name);
      try {
        process_pressure_spec(spec, bc);
      } catch (Errors::Message &msg) {
        Errors::Message m;
        m << "in sublist \"" << spec.name().c_str() << "\": " << msg.what();
        Exceptions::amanzi_throw(m);
      }
    } else { // ERROR -- parameter is not a sublist
      Errors::Message m;
      m << "parameter \"" << name.c_str() << "\" is not a sublist";
      Exceptions::amanzi_throw(m);
    }
  }
}


/* ******************************************************************
* TBA
****************************************************************** */
void FlowBCFactory::process_pressure_spec(Teuchos::ParameterList &list, BoundaryFunction *bc) const
{
  Errors::Message m;
  // Get the regions parameter value.
  std::vector<std::string> regions;
  if (list.isParameter("regions")) {
    if (list.isType<Teuchos::Array<std::string> >("regions")) {
      regions = list.get<Teuchos::Array<std::string> >("regions").toVector();
    } else { // ERROR -- "regions" is not of type Array<int>
      m << "parameter \"regions\" is not of type \"Array string\"";
      Exceptions::amanzi_throw(m);
    }
  } else { // ERROR -- "regions" is missing
    m << "parameter \"regions\" is missing";
    Exceptions::amanzi_throw(m);
  }
  // Get the boundary pressure function sublist.
  Teuchos::ParameterList *f_list;
  if (list.isParameter("boundary pressure")) {
    if (list.isSublist("boundary pressure")) {
      // validate function-factory sublist
      f_list = &list.sublist("boundary pressure");
    } else { // ERROR -- "boundary pressure" is not a sublist
      m << "parameter \"boundary pressure\" is not a sublist";
      Exceptions::amanzi_throw(m);
    }
  } else { // ERROR -- "boundary pressure" sublist is missing
    m << "sublist \"boundary pressure\" is missing";
    Exceptions::amanzi_throw(m);
  }
  // Make the boundary pressure function.
  Teuchos::RCP<Function> f;
  FunctionFactory f_fact;
  try {
    f = Teuchos::rcp(f_fact.Create(*f_list));
  } catch (Errors::Message &msg) {
    m << "error in sublist \"boundary pressure\": " << msg.what();
    Exceptions::amanzi_throw(m);
  }
  // Add this BC specification to the boundary function.
  bc->Define(regions, f);
}


/* ******************************************************************
* TBA
****************************************************************** */
void FlowBCFactory::process_mass_flux_list(Teuchos::ParameterList &list,
                                           BoundaryFunction *bc) const
{
  // Iterate through the BC specification sublists in the list.
  // All are expected to be sublists of identical structure.
  for (Teuchos::ParameterList::ConstIterator i = list.begin(); i != list.end(); ++i) {
    std::string name = i->first;
    if (list.isSublist(name)) {
      Teuchos::ParameterList &spec = list.sublist(name);
      try {
        process_mass_flux_spec(spec, bc);
      } catch (Errors::Message &msg) {
        Errors::Message m;
        m << "in sublist \"" << spec.name().c_str() << "\": " << msg.what();
        Exceptions::amanzi_throw(m);
      }
    } else { // ERROR -- parameter is not a sublist
      Errors::Message m;
      m << "parameter \"" << name.c_str() << "\" is not a sublist";
      Exceptions::amanzi_throw(m);
    }
  }
}


/* ******************************************************************
* TBA
****************************************************************** */
void FlowBCFactory::process_mass_flux_spec(Teuchos::ParameterList &list,
                                           BoundaryFunction *bc) const
{
  Errors::Message m;
  // Get the regions parameter value.
  std::vector<std::string> regions;
  if (list.isParameter("regions")) {
    if (list.isType<Teuchos::Array<std::string> >("regions")) {
      regions = list.get<Teuchos::Array<std::string> >("regions").toVector();
    } else { // ERROR -- "regions" is not of type Array<int>
      m << "parameter \"regions\" is not of type \"Array string\"";
      Exceptions::amanzi_throw(m);
    }
  } else { // ERROR -- "regions" is missing
    m << "parameter \"regions\" is missing";
    Exceptions::amanzi_throw(m);
  }
  // Get the outward mass flux function sublist.
  Teuchos::ParameterList *f_list;
  if (list.isParameter("outward mass flux")) {
    if (list.isSublist("outward mass flux")) {
      // validate function-factory sublist
      f_list = &list.sublist("outward mass flux");
    } else { // ERROR -- "outward mass flux" is not a sublist
      m << "parameter \"outward mass flux\" is not a sublist";
      Exceptions::amanzi_throw(m);
    }
  } else { // ERROR -- "outward mass flux" sublist is missing
    m << "sublist \"outward mass flux\" is missing";
    Exceptions::amanzi_throw(m);
  }
  // Make the outward mass flux function.
  Teuchos::RCP<Function> f;
  FunctionFactory f_fact;
  try {
    f = Teuchos::rcp(f_fact.Create(*f_list));
  } catch (Errors::Message &msg) {
    m << "error in sublist \"outward mass flux\": " << msg.what();
    Exceptions::amanzi_throw(m);
  }
  // Add this BC specification to the boundary function.
  bc->Define(regions, f);
}


/* ******************************************************************
* TBA
****************************************************************** */
void FlowBCFactory::process_static_head_list(
    double p0, double density, double gravity,
    Teuchos::ParameterList &list, BoundaryFunction *bc) const
{
  // Iterate through the BC specification sublists in the list.
  // All are expected to be sublists of identical structure.
  for (Teuchos::ParameterList::ConstIterator i = list.begin(); i != list.end(); ++i) {
    std::string name = i->first;
    if (list.isSublist(name)) {
      Teuchos::ParameterList &spec = list.sublist(name);
      try {
        process_static_head_spec(p0, density, gravity, spec, bc);
      } catch (Errors::Message &msg) {
        Errors::Message m;
        m << "in sublist \"" << spec.name().c_str() << "\": " << msg.what();
        Exceptions::amanzi_throw(m);
      }
    } else { // ERROR -- parameter is not a sublist
      Errors::Message m;
      m << "parameter \"" << name.c_str() << "\" is not a sublist";
      Exceptions::amanzi_throw(m);
    }
  }
}


/* ******************************************************************
* TBA
****************************************************************** */
void FlowBCFactory::process_static_head_spec(
    double p0, double density, double gravity,
    Teuchos::ParameterList &list, BoundaryFunction *bc) const
{
  Errors::Message m;
  // Get the regions parameter value.
  std::vector<std::string> regions;
  if (list.isParameter("regions")) {
    if (list.isType<Teuchos::Array<std::string> >("regions")) {
      regions = list.get<Teuchos::Array<std::string> >("regions").toVector();
    } else { // ERROR -- "regions" is not of type Array<int>
      m << "parameter \"regions\" is not of type \"Array string\"";
      Exceptions::amanzi_throw(m);
    }
  } else { // ERROR -- "regions" is missing
    m << "parameter \"regions\" is missing";
    Exceptions::amanzi_throw(m);
  }
  // Get the water table elevation function sublist.
  Teuchos::ParameterList *water_table_list;
  if (list.isParameter("water table elevation")) {
    if (list.isSublist("water table elevation")) {
      // validate function-factory sublist
      water_table_list = &list.sublist("water table elevation");
    } else { // ERROR -- "water table elevation" is not a sublist
      m << "parameter \"water table elevation\" is not a sublist";
      Exceptions::amanzi_throw(m);
    }
  } else { // ERROR -- "water table elevation" sublist is missing
    m << "sublist \"water table elevation\" is missing";
    Exceptions::amanzi_throw(m);
  }
  // Form the parameter list to create static head function.
  Teuchos::ParameterList f_list;
  Teuchos::ParameterList &static_head_list = f_list.sublist("function-static-head");
  static_head_list.set("p0", p0);
  static_head_list.set("density", density);
  static_head_list.set("gravity", gravity);
  static_head_list.set("water table elevation", *water_table_list);
  // Make the static head function.
  Teuchos::RCP<Function> f;
  FunctionFactory f_fact;
  try {
    f = Teuchos::rcp(f_fact.Create(f_list));
  } catch (Errors::Message &msg) {
    m << "error in sublist \"water table elevation\": " << msg.what();
    Exceptions::amanzi_throw(m);
  }
  // Add this BC specification to the boundary function.
  bc->Define(regions, f);
}

}  // namespace AmanziFlow
}  // namespace Amanzi
