/*
This is the flow component of the Amanzi code. 
License: BSD
Authors: Neil Carlson (versions 1 & 2)  (nnc@lanl.gov)
         Ethan Coon (ATS version)
*/

#include "constant-function.hh"
#include "function-factory.hh"
#include "vector_function.hh"
#include "composite_function.hh"
#include "errors.hh"

#include "flow_bc_factory.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* Process Dirichet BC (pressure), step 1.
****************************************************************** */
Teuchos::RCP<Functions::BoundaryFunction> FlowBCFactory::CreatePressure() const {
  Teuchos::RCP<Functions::BoundaryFunction> bc = Teuchos::rcp(new Functions::BoundaryFunction(mesh_));
  try {
    ProcessPressureList(plist_.sublist("pressure"), bc);
  } catch (Errors::Message& msg) {
    Errors::Message m;
    m << "FlowBCFactory: \"pressure\" sublist error: " << msg.what();
    Exceptions::amanzi_throw(m);
  } catch (Teuchos::Exceptions::InvalidParameterType& msg) {
    Errors::Message m;
    m << "FlowBCFactory: \"pressure\" sublist error: not a sublist: " << msg.what();
    Exceptions::amanzi_throw(m);
  }
  return bc;
}


/* ******************************************************************
* Process Neumann BC (mass flux), step 1.
****************************************************************** */
Teuchos::RCP<Functions::BoundaryFunction> FlowBCFactory::CreateMassFlux() const {
  Teuchos::RCP<Functions::BoundaryFunction> bc = Teuchos::rcp(new Functions::BoundaryFunction(mesh_));
  try {
    ProcessMassFluxList(plist_.sublist("mass flux"), bc);
  } catch (Errors::Message& msg) {
    Errors::Message m;
    m << "FlowBCFactory: \"mass flux\" sublist error: " << msg.what();
    Exceptions::amanzi_throw(m);
  } catch (Teuchos::Exceptions::InvalidParameterType& msg) {
    Errors::Message m;
    m << "FlowBCFactory: \"mass flux\" sublist error: not a sublist: " << msg.what();
    Exceptions::amanzi_throw(m);
  }
  return bc;
}


/* ******************************************************************
* Process Dirichet BC (static head), step 1.
****************************************************************** */
Teuchos::RCP<Functions::BoundaryFunction> FlowBCFactory::CreateStaticHead(double p0,
        double density, AmanziGeometry::Point& gravity) const {
  Teuchos::RCP<Functions::BoundaryFunction> bc = Teuchos::rcp(new Functions::BoundaryFunction(mesh_));
  try {
    ProcessStaticHeadList(p0, density, gravity, plist_.sublist("static head"), bc);
  } catch (Errors::Message& msg) {
    Errors::Message m;
    m << "FlowBCFactory: \"static head\" sublist error: " << msg.what();
    Exceptions::amanzi_throw(m);
  } catch (Teuchos::Exceptions::InvalidParameterType& msg) {
    Errors::Message m;
    m << "FlowBCFactory: \"static head\" sublist error: not a sublist: " << msg.what();
    Exceptions::amanzi_throw(m);
  }
  return bc;
}


/* ******************************************************************
* Process Zero Gradient BC (pressure), step 1.
****************************************************************** */
Teuchos::RCP<Functions::BoundaryFunction> FlowBCFactory::CreateZeroGradient() const {
  Teuchos::RCP<Functions::BoundaryFunction> bc = Teuchos::rcp(new Functions::BoundaryFunction(mesh_));
  try {
    ProcessZeroGradientList(plist_.sublist("zero-gradient"), bc);
  } catch (Errors::Message& msg) {
    Errors::Message m;
    m << "FlowBCFactory: \"zero-gradient\" sublist error: " << msg.what();
    Exceptions::amanzi_throw(m);
  } catch (Teuchos::Exceptions::InvalidParameterType& msg) {
    Errors::Message m;
    m << "FlowBCFactory: \"zero-gradient\" sublist error: not a sublist: " << msg.what();
    Exceptions::amanzi_throw(m);
  }
  return bc;
}


/* ******************************************************************
* Process Dirichet BC (pressure), step 2.
****************************************************************** */
void FlowBCFactory::ProcessPressureList(const Teuchos::ParameterList& list,
                                        const Teuchos::RCP<Functions::BoundaryFunction>& bc) const {
  // Iterate through the BC specification sublists in the list.
  // All are expected to be sublists of identical structure.
  for (Teuchos::ParameterList::ConstIterator i = list.begin(); i != list.end(); ++i) {
    std::string name = i->first;
    if (list.isSublist(name)) {
      Teuchos::ParameterList spec = list.sublist(name);
      try {
        ProcessPressureSpec(spec, bc);
      } catch (Errors::Message& msg) {
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
* Process Dirichet BC (pressure), step 3.
****************************************************************** */
void FlowBCFactory::ProcessPressureSpec(const Teuchos::ParameterList& list,
        const Teuchos::RCP<Functions::BoundaryFunction>& bc) const {
  Errors::Message m;
  std::vector<std::string> regions;

  if (list.isParameter("regions")) {
    if (list.isType<Teuchos::Array<std::string> >("regions")) {
      regions = list.get<Teuchos::Array<std::string> >("regions").toVector();
    } else {
      m << "parameter \"regions\" is not of type \"Array string\"";
      Exceptions::amanzi_throw(m);
    }
  } else {  // Parameter "regions" is missing.
    m << "parameter \"regions\" is missing";
    Exceptions::amanzi_throw(m);
  }

  Teuchos::ParameterList f_list;
  if (list.isParameter("boundary pressure")) {
    if (list.isSublist("boundary pressure")) {
      // validate function-factory sublist
      f_list = list.sublist("boundary pressure");
    } else {  // "boundary pressure" is not a sublist
      m << "parameter \"boundary pressure\" is not a sublist";
      Exceptions::amanzi_throw(m);
    }
  } else {  // "boundary pressure" sublist is missing.
    m << "sublist \"boundary pressure\" is missing";
    Exceptions::amanzi_throw(m);
  }

  // Make the boundary pressure function.
  Teuchos::RCP<Function> f;
  FunctionFactory f_fact;
  try {
    f = Teuchos::rcp(f_fact.Create(f_list));
  } catch (Errors::Message& msg) {
    m << "error in sublist \"boundary pressure\": " << msg.what();
    Exceptions::amanzi_throw(m);
  }

  // A bit hacky -- this entire code needs to be revisited in light of the new
  // options for mesh_functions.
  Teuchos::RCP<VectorFunction> func = Teuchos::rcp(new CompositeFunction(f));

  // Add this BC specification to the boundary function.
  bc->Define(regions, func);
}


/* ******************************************************************
* Process Neumann BC (mass flux), step 2.
****************************************************************** */
void FlowBCFactory::ProcessMassFluxList(const Teuchos::ParameterList& list,
                                        const Teuchos::RCP<Functions::BoundaryFunction>& bc) const {
  // Iterate through the BC specification sublists in the list.
  // All are expected to be sublists of identical structure.
  for (Teuchos::ParameterList::ConstIterator i = list.begin(); i != list.end(); ++i) {
    std::string name = i->first;
    if (list.isSublist(name)) {
      Teuchos::ParameterList spec = list.sublist(name);
      try {
        ProcessMassFluxSpec(spec, bc);
      } catch (Errors::Message& msg) {
        Errors::Message m;
        m << "in sublist \"" << spec.name().c_str() << "\": " << msg.what();
        Exceptions::amanzi_throw(m);
      }
    } else {  // Parameter is not a sublist
      Errors::Message m;
      m << "parameter \"" << name.c_str() << "\" is not a sublist";
      Exceptions::amanzi_throw(m);
    }
  }
}


/* ******************************************************************
* Process Neumann BC (mass flux), step 3.
****************************************************************** */
void FlowBCFactory::ProcessMassFluxSpec(const Teuchos::ParameterList& list,
                                        const Teuchos::RCP<Functions::BoundaryFunction>& bc) const {
  Errors::Message m;
  // Get the regions parameter value.
  std::vector<std::string> regions;
  if (list.isParameter("regions")) {
    if (list.isType<Teuchos::Array<std::string> >("regions")) {
      regions = list.get<Teuchos::Array<std::string> >("regions").toVector();
    } else {
      m << "parameter \"regions\" is not of type \"Array string\"";
      Exceptions::amanzi_throw(m);
    }
  } else {  // Parameter "regions" is missing.
    m << "parameter \"regions\" is missing";
    Exceptions::amanzi_throw(m);
  }

  Teuchos::ParameterList f_list;
  if (list.isParameter("outward mass flux")) {
    if (list.isSublist("outward mass flux")) {
      // validate function-factory sublist
      f_list = list.sublist("outward mass flux");
    } else {  // "outward mass flux" is not a sublist.
      m << "parameter \"outward mass flux\" is not a sublist";
      Exceptions::amanzi_throw(m);
    }
  } else {  // "outward mass flux" sublist is missing
    m << "sublist \"outward mass flux\" is missing";
    Exceptions::amanzi_throw(m);
  }

  // Make the outward mass flux function.
  Teuchos::RCP<Function> f;
  FunctionFactory f_fact;
  try {
    f = Teuchos::rcp(f_fact.Create(f_list));
  } catch (Errors::Message& msg) {
    m << "error in sublist \"outward mass flux\": " << msg.what();
    Exceptions::amanzi_throw(m);
  }

  // A bit hacky -- this entire code needs to be revisited in light of the new
  // options for mesh_functions.
  Teuchos::RCP<VectorFunction> func = Teuchos::rcp(new CompositeFunction(f));

  // Add this BC specification to the boundary function.
  bc->Define(regions, func);
}


/* ******************************************************************
* Process Dirichet BC (static head), step 2.
****************************************************************** */
void FlowBCFactory::ProcessStaticHeadList(double p0, double density,
        AmanziGeometry::Point& gravity, const Teuchos::ParameterList& list,
        const Teuchos::RCP<Functions::BoundaryFunction>& bc) const {
  // Iterate through the BC specification sublists in the list.
  // All are expected to be sublists of identical structure.
  for (Teuchos::ParameterList::ConstIterator i = list.begin(); i != list.end(); ++i) {
    std::string name = i->first;
    if (list.isSublist(name)) {
      Teuchos::ParameterList spec = list.sublist(name);
      try {
        ProcessStaticHeadSpec(p0, density, gravity, spec, bc);
      } catch (Errors::Message& msg) {
        Errors::Message m;
        m << "in sublist \"" << spec.name().c_str() << "\": " << msg.what();
        Exceptions::amanzi_throw(m);
      }
    } else {  // Parameter is not a sublist.
      Errors::Message m;
      m << "parameter \"" << name.c_str() << "\" is not a sublist";
      Exceptions::amanzi_throw(m);
    }
  }
}


/* ******************************************************************
* Process Dirichet BC (static head), step 3.
****************************************************************** */
void FlowBCFactory::ProcessStaticHeadSpec(double p0, double density,
        AmanziGeometry::Point& gravity, const Teuchos::ParameterList& list,
        const Teuchos::RCP<Functions::BoundaryFunction>& bc) const {
  Errors::Message m;
  std::vector<std::string> regions;

  if (list.isParameter("regions")) {
    if (list.isType<Teuchos::Array<std::string> >("regions")) {
      regions = list.get<Teuchos::Array<std::string> >("regions").toVector();
    } else {
      m << "parameter \"regions\" is not of type \"Array string\"";
      Exceptions::amanzi_throw(m);
    }
  } else {  // Parameter "regions" is missing.
    m << "parameter \"regions\" is missing";
    Exceptions::amanzi_throw(m);
  }

  // Get the water table elevation function sublist.
  Teuchos::ParameterList water_table_list;
  if (list.isParameter("water table elevation")) {
    if (list.isSublist("water table elevation")) {
      // validate function-factory sublist
      water_table_list = list.sublist("water table elevation");
    } else {  // "water table elevation" is not a sublist
      m << "parameter \"water table elevation\" is not a sublist";
      Exceptions::amanzi_throw(m);
    }
  } else {  // "water table elevation" sublist is missing.
    m << "sublist \"water table elevation\" is missing";
    Exceptions::amanzi_throw(m);
  }

  // Form the parameter list to create static head function.
  Teuchos::ParameterList f_list;
  Teuchos::ParameterList& static_head_list = f_list.sublist("function-static-head");
  int dim = gravity.dim();

  static_head_list.set("p0", p0);
  static_head_list.set("density", density);
  static_head_list.set("gravity", -gravity[dim-1]);
  static_head_list.set("space dimension", dim);
  static_head_list.set("water table elevation", water_table_list);

  Teuchos::RCP<Function> f;  // Make the static head function.
  FunctionFactory f_fact;
  try {
    f = Teuchos::rcp(f_fact.Create(f_list));
  } catch (Errors::Message& msg) {
    m << "error in sublist \"water table elevation\": " << msg.what();
    Exceptions::amanzi_throw(m);
  }

  // A bit hacky -- this entire code needs to be revisited in light of the new
  // options for mesh_functions.
  Teuchos::RCP<VectorFunction> func = Teuchos::rcp(new CompositeFunction(f));

  // Add this BC specification to the boundary function.
  bc->Define(regions, func);
}


/* ******************************************************************
* Process Dirichet BC (pressure), step 2.
****************************************************************** */
void FlowBCFactory::ProcessZeroGradientList(const Teuchos::ParameterList& list,
                                        const Teuchos::RCP<Functions::BoundaryFunction>& bc) const {
  // Iterate through the BC specification sublists in the list.
  // All are expected to be sublists of identical structure.
  for (Teuchos::ParameterList::ConstIterator i = list.begin(); i != list.end(); ++i) {
    std::string name = i->first;
    if (list.isSublist(name)) {
      Teuchos::ParameterList spec = list.sublist(name);
      try {
        ProcessZeroGradientSpec(spec, bc);
      } catch (Errors::Message& msg) {
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
* Process Dirichet BC (pressure), step 3.
****************************************************************** */
void FlowBCFactory::ProcessZeroGradientSpec(const Teuchos::ParameterList& list,
        const Teuchos::RCP<Functions::BoundaryFunction>& bc) const {
  Errors::Message m;
  std::vector<std::string> regions;

  if (list.isParameter("regions")) {
    if (list.isType<Teuchos::Array<std::string> >("regions")) {
      regions = list.get<Teuchos::Array<std::string> >("regions").toVector();
    } else {
      m << "parameter \"regions\" is not of type \"Array string\"";
      Exceptions::amanzi_throw(m);
    }
  } else {  // Parameter "regions" is missing.
    m << "parameter \"regions\" is missing";
    Exceptions::amanzi_throw(m);
  }

  // Make the boundary function.  This simply puts in zero -- the PK
  // is expected to handle this condition correctly.
  Teuchos::RCP<Function> f = Teuchos::rcp(new ConstantFunction(0.0));

  // A bit hacky -- this entire code needs to be revisited in light of the new
  // options for mesh_functions.
  Teuchos::RCP<VectorFunction> func = Teuchos::rcp(new CompositeFunction(f));

  // Add this BC specification to the boundary function.
  bc->Define(regions, func);
}
}  // namespace AmanziFlow
}  // namespace Amanzi
