/*
This is the flow component of the Amanzi code. 
License: BSD
Authors: Neil Carlson (version 1)  (nnc@lanl.gov)
         Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

#include <vector>
#include <string>

#include "function-factory.hh"
#include "errors.hh"

#include "Flow_BC_Factory.hpp"

namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* Process Dirichet BC (pressure), step 1.
****************************************************************** */
BoundaryFunction* FlowBCFactory::CreatePressure() const
{
  BoundaryFunction* bc = new BoundaryFunction(mesh_);
  try {
    ProcessPressureList(params_->sublist("pressure"), bc);
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
BoundaryFunction* FlowBCFactory::CreateMassFlux(std::vector<double>& rainfall_factor) const
{
  int ncells = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  rainfall_factor.resize(ncells, 1.0);

  BoundaryFunction* bc = new BoundaryFunction(mesh_);
  try {
    ProcessMassFluxList(params_->sublist("mass flux"), rainfall_factor, bc);
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
BoundaryFunction* FlowBCFactory::CreateStaticHead(
    double p0, double density, AmanziGeometry::Point& gravity) const
{
  BoundaryFunction* bc = new BoundaryFunction(mesh_);
  try {
    ProcessStaticHeadList(p0, density, gravity, params_->sublist("static head"), bc);
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
* Seepage Face BC, step 1.
****************************************************************** */
BoundaryFunction* FlowBCFactory::CreateSeepageFace() const
{
  BoundaryFunction* bc = new BoundaryFunction(mesh_);
  try {
    ProcessSeepageFaceList(params_->sublist("seepage face"), bc);
  } catch (Errors::Message& msg) {
    Errors::Message m;
    m << "FlowBCFactory: \"seepage face\" sublist error: " << msg.what();
    Exceptions::amanzi_throw(m);
  } catch (Teuchos::Exceptions::InvalidParameterType& msg) {
    Errors::Message m;
    m << "FlowBCFactory: \"seepage face\" sublist error: not a sublist: " << msg.what();
    Exceptions::amanzi_throw(m);
  }
  return bc;
}


/* ******************************************************************
* Process Dirichet BC (pressure), step 2.
* Loop over sublists with typical names "BC 0", "BC 1", etc.
****************************************************************** */
void FlowBCFactory::ProcessPressureList(Teuchos::ParameterList& list,
                                        BoundaryFunction* bc) const
{
  // Iterate through the BC specification sublists in the list.
  // All are expected to be sublists of identical structure.
  for (Teuchos::ParameterList::ConstIterator i = list.begin(); i != list.end(); ++i) {
    std::string name = i->first;
    if (list.isSublist(name)) {
      Teuchos::ParameterList& spec = list.sublist(name);
      try {
        ProcessPressureSpec(spec, bc);
      } catch (Errors::Message& msg) {
        Errors::Message m;
        m << "in sublist \"" << spec.name().c_str() << "\": " << msg.what();
        Exceptions::amanzi_throw(m);
      }
    } else {
      Errors::Message m;
      m << "parameter \"" << name.c_str() << "\" is not a sublist";
      Exceptions::amanzi_throw(m);
    }
  }
}


/* ******************************************************************
* Process Dirichet BC (pressure), step 3.
****************************************************************** */
void FlowBCFactory::ProcessPressureSpec(Teuchos::ParameterList& list, BoundaryFunction* bc) const
{
  Errors::Message m;
  std::vector<std::string> regions;

  if (list.isParameter("regions")) {
    if (list.isType<Teuchos::Array<std::string> >("regions")) {
      regions = list.get<Teuchos::Array<std::string> >("regions").toVector();
    } else {
      m << "parameter \"regions\" is not of type \"Array string\"";
      Exceptions::amanzi_throw(m);
    }
  } else {
    m << "parameter \"regions\" is missing";
    Exceptions::amanzi_throw(m);
  }

  Teuchos::ParameterList* f_list;
  if (list.isSublist("boundary pressure")) {
    f_list = &list.sublist("boundary pressure");
  } else {  // "boundary pressure" is not a sublist
    m << "parameter \"boundary pressure\" is not a sublist";
    Exceptions::amanzi_throw(m);
  }

  // Make the boundary pressure function.
  Teuchos::RCP<Function> f;
  FunctionFactory f_fact;
  try {
    f = Teuchos::rcp(f_fact.Create(*f_list));
  } catch (Errors::Message& msg) {
    m << "error in sublist \"boundary pressure\": " << msg.what();
    Exceptions::amanzi_throw(m);
  }

  // Add this BC specification to the boundary function.
  bc->Define(regions, f);
}


/* ******************************************************************
* Process Neumann BC (mass flux), step 2.
* Iterate through the BC specification sublists in the list with 
* typical names "BC 0", "BC 1", etc. All are expected to be sublists 
* of identical structure.
****************************************************************** */
void FlowBCFactory::ProcessMassFluxList(
    Teuchos::ParameterList& list, 
    std::vector<double>& rainfall_factor, BoundaryFunction* bc) const
{
  for (Teuchos::ParameterList::ConstIterator i = list.begin(); i != list.end(); ++i) {
    std::string name = i->first;
    if (list.isSublist(name)) {
      Teuchos::ParameterList& spec = list.sublist(name);
      try {
        ProcessMassFluxSpec(spec, rainfall_factor, bc);
      } catch (Errors::Message& msg) {
        Errors::Message m;
        m << "in sublist \"" << spec.name().c_str() << "\": " << msg.what();
        Exceptions::amanzi_throw(m);
      }
    } else {
      Errors::Message m;
      m << "parameter \"" << name.c_str() << "\" is not a sublist";
      Exceptions::amanzi_throw(m);
    }
  }
}


/* ******************************************************************
* Process Neumann BC (mass flux), step 3.
****************************************************************** */
void FlowBCFactory::ProcessMassFluxSpec(
    Teuchos::ParameterList& list,
    std::vector<double>& rainfall_factor, BoundaryFunction* bc) const
{
  Errors::Message m;
  std::vector<std::string> regions;

  if (list.isParameter("regions")) {
    if (list.isType<Teuchos::Array<std::string> >("regions")) {
      regions = list.get<Teuchos::Array<std::string> >("regions").toVector();
    } else {
      m << "parameter \"regions\" is not of type \"Array string\"";
      Exceptions::amanzi_throw(m);
    }
  } else {
    m << "parameter \"regions\" is missing";
    Exceptions::amanzi_throw(m);
  }

  Teuchos::ParameterList* f_list;
  if (list.isSublist("outward mass flux")) {
    f_list = &list.sublist("outward mass flux");
  } else {
    m << "parameter \"outward mass flux\" is not a sublist";
    Exceptions::amanzi_throw(m);
  }

  // Make the outward mass flux function.
  Teuchos::RCP<Function> f;
  FunctionFactory f_fact;
  try {
    f = Teuchos::rcp(f_fact.Create(*f_list));
  } catch (Errors::Message& msg) {
    m << "error in sublist \"outward mass flux\": " << msg.what();
    Exceptions::amanzi_throw(m);
  }

  // Add this BC specification to the boundary function.
  bc->Define(regions, f);

  // Calculate rainfall factor
  if (list.get<bool>("rainfall", false)) {
    int dim = mesh_->space_dimension();

    int nregions = regions.size();
    for (int n = 0; n < nregions; n++) {
      AmanziMesh::Entity_ID_List faces;
      mesh_->get_set_entities(regions[n], AmanziMesh::FACE, AmanziMesh::OWNED, &faces);

      int nfaces = faces.size();
      for (int m = 0; m < nfaces; m++) {
        int f = faces[m];
        const AmanziGeometry::Point& normal = mesh_->face_normal(f);
        rainfall_factor[f] = fabs(normal[dim - 1]) / mesh_->face_area(f);
      }
    }
  }
}


/* ******************************************************************
* Process Dirichet BC (static head), step 2.
* Iterate through the BC specification sublists with typical names 
* "BC 0", "BC 1", etc. All are expected to be sublists of identical 
* structure.
****************************************************************** */
void FlowBCFactory::ProcessStaticHeadList(
    double p0, double density, AmanziGeometry::Point& gravity,
    Teuchos::ParameterList& list, BoundaryFunction* bc) const
{
  for (Teuchos::ParameterList::ConstIterator i = list.begin(); i != list.end(); ++i) {
    std::string name = i->first;
    if (list.isSublist(name)) {
      Teuchos::ParameterList& spec = list.sublist(name);
      try {
        ProcessStaticHeadSpec(p0, density, gravity, spec, bc);
      } catch (Errors::Message& msg) {
        Errors::Message m;
        m << "in sublist \"" << spec.name().c_str() << "\": " << msg.what();
        Exceptions::amanzi_throw(m);
      }
    } else {
      Errors::Message m;
      m << "parameter \"" << name.c_str() << "\" is not a sublist";
      Exceptions::amanzi_throw(m);
    }
  }
}


/* ******************************************************************
* Process Dirichet BC (static head), step 3.
****************************************************************** */
void FlowBCFactory::ProcessStaticHeadSpec(
    double p0, double density, AmanziGeometry::Point& gravity,
    Teuchos::ParameterList& list, BoundaryFunction* bc) const
{
  Errors::Message m;
  std::vector<std::string> regions;

  if (list.isParameter("regions")) {
    if (list.isType<Teuchos::Array<std::string> >("regions")) {
      regions = list.get<Teuchos::Array<std::string> >("regions").toVector();
    } else {
      m << "parameter \"regions\" is not of type \"Array string\"";
      Exceptions::amanzi_throw(m);
    }
  } else {
    m << "parameter \"regions\" is missing";
    Exceptions::amanzi_throw(m);
  }

  // Get the water table elevation function sublist.
  Teuchos::ParameterList* water_table_list;
  if (list.isSublist("water table elevation")) {
    water_table_list = &list.sublist("water table elevation");
  } else {
    m << "parameter \"water table elevation\" is not a sublist";
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
  static_head_list.set("water table elevation", *water_table_list);

  Teuchos::RCP<Function> f;  // Make the static head function.
  FunctionFactory f_fact;
  try {
    f = Teuchos::rcp(f_fact.Create(f_list));
  } catch (Errors::Message& msg) {
    m << "error in sublist \"water table elevation\": " << msg.what();
    Exceptions::amanzi_throw(m);
  }

  // Add this BC specification to the boundary function.
  bc->Define(regions, f);
}


/* ******************************************************************
* Process Seepage Face BC, step 2.
****************************************************************** */
void FlowBCFactory::ProcessSeepageFaceList(Teuchos::ParameterList& list,
                                           BoundaryFunction* bc) const
{
  // Iterate through the BC specification sublists in the list.
  // All are expected to be sublists of identical structure.
  for (Teuchos::ParameterList::ConstIterator i = list.begin(); i != list.end(); ++i) {
    std::string name = i->first;
    if (list.isSublist(name)) {
      Teuchos::ParameterList& spec = list.sublist(name);
      try {
        ProcessSeepageFaceSpec(spec, bc);
      } catch (Errors::Message& msg) {
        Errors::Message m;
        m << "in sublist \"" << spec.name().c_str() << "\": " << msg.what();
        Exceptions::amanzi_throw(m);
      }
    } else {  // ERROR -- parameter is not a sublist
      Errors::Message m;
      m << "parameter \"" << name.c_str() << "\" is not a sublist";
      Exceptions::amanzi_throw(m);
    }
  }
}


/* ******************************************************************
* Process Seepage Face BC, step 3.
****************************************************************** */
void FlowBCFactory::ProcessSeepageFaceSpec(Teuchos::ParameterList& list,
                                           BoundaryFunction* bc) const
{
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

  Teuchos::ParameterList* f_list;
  if (list.isSublist("outward mass flux")) {
    f_list = &list.sublist("outward mass flux");
  } else {
    m << "parameter \"outward mass flux\" is not a sublist";
    Exceptions::amanzi_throw(m);
  }

  // Make the seepage face function.
  Teuchos::RCP<Function> f;
  FunctionFactory f_fact;
  try {
    f = Teuchos::rcp(f_fact.Create(*f_list));
  } catch (Errors::Message& msg) {
    m << "error in sublist \"outward mass flux\": " << msg.what();
    Exceptions::amanzi_throw(m);
  }

  // Add this BC specification to the boundary function.
  bc->Define(regions, f);
}

}  // namespace AmanziFlow
}  // namespace Amanzi
