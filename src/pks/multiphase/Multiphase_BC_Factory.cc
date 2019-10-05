/*
  This is the multiphase component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Quan Bui (mquanbui@math.umd.edu)
  Adaptation of Flow_BC_Factory for multiphase. We need to prescribe boundary
  conditions for both pressure and saturation. Flow_BC_Factory only has option
  for pressure.
*/

#include <string>
#include <sstream>
#include <vector>

#include "errors.hh"
#include "MultiFunction.hh"

#include "Multiphase_BC_Factory.hh"
#include "MultiphaseDefs.hh"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* Constructor
****************************************************************** */
MultiphaseBCFactory::MultiphaseBCFactory(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                             const Teuchos::RCP<Teuchos::ParameterList>& plist)
   : mesh_(mesh), plist_(plist)
{
  // create verbosity object
  Teuchos::ParameterList vlist;
  // vlist.set("Verbosity Level", "medium");
  vo_ = new VerboseObject("MultiphaseBCFactory", vlist); 
}


/* ******************************************************************
* Process Dirichet BC (pressure), step 1.
****************************************************************** */
Flow::FlowBoundaryFunction* MultiphaseBCFactory::CreatePressure(std::vector<int>& submodel) const
{
  Flow::FlowBoundaryFunction* bc = new Flow::FlowBoundaryFunction(mesh_);
  try {
    ProcessPressureList(plist_->sublist("pressure"), submodel, bc);
  } catch (Errors::Message& msg) {
    Errors::Message m;
    m << "MultiphaseBCFactory: \"pressure\" sublist error: " << msg.what();
    Exceptions::amanzi_throw(m);
  } catch (Teuchos::Exceptions::InvalidParameterType& msg) {
    Errors::Message m;
    m << "MultiphaseBCFactory: \"pressure\" sublist error: not a sublist: " << msg.what();
    Exceptions::amanzi_throw(m);
  }
  return bc;
}


/* ******************************************************************
* Process Dirichet BC (pressure), step 2.
* Loop over sublists with typical names "BC 0", "BC 1", etc.
****************************************************************** */
void MultiphaseBCFactory::ProcessPressureList(
    Teuchos::ParameterList& list, std::vector<int>& submodel, Flow::FlowBoundaryFunction* bc) const
{
  // Iterate through the BC specification sublists in the list.
  // All are expected to be sublists of identical structure.
  for (Teuchos::ParameterList::ConstIterator i = list.begin(); i != list.end(); ++i) {
    std::string name = i->first;
    if (list.isSublist(name)) {
      Teuchos::ParameterList& spec = list.sublist(name);
      try {
        ProcessPressureSpec(spec, submodel, bc);
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
void MultiphaseBCFactory::ProcessPressureSpec(
    Teuchos::ParameterList& list, std::vector<int>& submodel, Flow::FlowBoundaryFunction* bc) const
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
  Teuchos::RCP<Amanzi::MultiFunction> f;
  try {
    f = Teuchos::rcp(new Amanzi::MultiFunction(*f_list));
  } catch (Errors::Message& msg) {
    m << "error in sublist \"boundary pressure\": " << msg.what();
    Exceptions::amanzi_throw(m);
  }

  // Add this BC specification to the boundary function.
  bc->Define(regions, f, Amanzi::CommonDefs::BOUNDARY_FUNCTION_ACTION_NONE);
}


/* ******************************************************************
* Process Dirichet BC (pressure), step 1.
****************************************************************** */
Flow::FlowBoundaryFunction* MultiphaseBCFactory::CreateHydrogenDensity(std::vector<int>& submodel) const
{
  Flow::FlowBoundaryFunction* bc = new Flow::FlowBoundaryFunction(mesh_);
  try {
    ProcessHydrogenDensityList(plist_->sublist("hydrogen density"), submodel, bc);
  } catch (Errors::Message& msg) {
    Errors::Message m;
    m << "MultiphaseBCFactory: \"hydrogen density\" sublist error: " << msg.what();
    Exceptions::amanzi_throw(m);
  } catch (Teuchos::Exceptions::InvalidParameterType& msg) {
    Errors::Message m;
    m << "MultiphaseBCFactory: \"hydrogen density\" sublist error: not a sublist: " << msg.what();
    Exceptions::amanzi_throw(m);
  }
  return bc;
}


/* ******************************************************************
* Process Dirichet BC (pressure), step 2.
* Loop over sublists with typical names "BC 0", "BC 1", etc.
****************************************************************** */
void MultiphaseBCFactory::ProcessHydrogenDensityList(
    Teuchos::ParameterList& list, std::vector<int>& submodel, Flow::FlowBoundaryFunction* bc) const
{
  // Iterate through the BC specification sublists in the list.
  // All are expected to be sublists of identical structure.
  for (Teuchos::ParameterList::ConstIterator i = list.begin(); i != list.end(); ++i) {
    std::string name = i->first;
    if (list.isSublist(name)) {
      Teuchos::ParameterList& spec = list.sublist(name);
      try {
        ProcessHydrogenDensitySpec(spec, submodel, bc);
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
void MultiphaseBCFactory::ProcessHydrogenDensitySpec(
    Teuchos::ParameterList& list, std::vector<int>& submodel, Flow::FlowBoundaryFunction* bc) const
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
  if (list.isSublist("boundary hydrogen density")) {
    f_list = &list.sublist("boundary hydrogen density");
  } else {  // "hydrogen density" is not a sublist
    m << "parameter \"boundary hydrogen density\" is not a sublist";
    Exceptions::amanzi_throw(m);
  }

  // Make the boundary pressure function.
  Teuchos::RCP<Amanzi::MultiFunction> f;
  try {
    f = Teuchos::rcp(new Amanzi::MultiFunction(*f_list));
  } catch (Errors::Message& msg) {
    m << "error in sublist \"boundary hydrogen density\": " << msg.what();
    Exceptions::amanzi_throw(m);
  }

  // Add this BC specification to the boundary function.
  bc->Define(regions, f, Amanzi::CommonDefs::BOUNDARY_FUNCTION_ACTION_NONE);
}


/* ******************************************************************
* Process Dirichet BC (saturation), step 1.
****************************************************************** */
Flow::FlowBoundaryFunction* MultiphaseBCFactory::CreateSaturation(std::vector<int>& submodel) const
{
  Flow::FlowBoundaryFunction* bc = new Flow::FlowBoundaryFunction(mesh_);
  try {
    ProcessSaturationList(plist_->sublist("saturation"), submodel, bc);
  } catch (Errors::Message& msg) {
    Errors::Message m;
    m << "MultiphaseBCFactory: \"saturation\" sublist error: " << msg.what();
    Exceptions::amanzi_throw(m);
  } catch (Teuchos::Exceptions::InvalidParameterType& msg) {
    Errors::Message m;
    m << "MultiphaseBCFactory: \"saturation\" sublist error: not a sublist: " << msg.what();
    Exceptions::amanzi_throw(m);
  }
  return bc;
}


/* ******************************************************************
* Process Dirichet BC (saturation), step 2.
* Loop over sublists with typical names "BC 0", "BC 1", etc.
****************************************************************** */
void MultiphaseBCFactory::ProcessSaturationList(
    Teuchos::ParameterList& list, std::vector<int>& submodel, Flow::FlowBoundaryFunction* bc) const
{
  // Iterate through the BC specification sublists in the list.
  // All are expected to be sublists of identical structure.
  for (Teuchos::ParameterList::ConstIterator i = list.begin(); i != list.end(); ++i) {
    std::string name = i->first;
    if (list.isSublist(name)) {
      Teuchos::ParameterList& spec = list.sublist(name);
      try {
        ProcessSaturationSpec(spec, submodel, bc);
      } catch (Errors::Message& msg) {Errors::Message m; m << "in sublist \"" << spec.name().c_str() << "\": " << msg.what();
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
* Process Dirichet BC (Saturation), step 3.
****************************************************************** */
void MultiphaseBCFactory::ProcessSaturationSpec(
    Teuchos::ParameterList& list, std::vector<int>& submodel, Flow::FlowBoundaryFunction* bc) const
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
  if (list.isSublist("boundary saturation")) {
    f_list = &list.sublist("boundary saturation");
  } else {  // "boundary saturation" is not a sublist
    m << "parameter \"boundary saturation\" is not a sublist";
    Exceptions::amanzi_throw(m);
  }

  // Make the boundary saturation function.
  Teuchos::RCP<Amanzi::MultiFunction> f;
  try {
    f = Teuchos::rcp(new Amanzi::MultiFunction(*f_list));
  } catch (Errors::Message& msg) {
    m << "error in sublist \"boundary saturation\": " << msg.what();
    Exceptions::amanzi_throw(m);
  }

  // Add this BC specification to the boundary function.
  bc->Define(regions, f, Amanzi::CommonDefs::BOUNDARY_FUNCTION_ACTION_NONE);
}


/* ******************************************************************
* Process Dirichet BC (molar fraction), step 1.
****************************************************************** */
Flow::FlowBoundaryFunction* MultiphaseBCFactory::CreateMolarFraction(std::vector<int>& submodel, int phase) const
{
  std::string sublist_name("molar fraction phase");
  std::stringstream ss;
  ss << phase;
  sublist_name = sublist_name + ss.str(); 
  //std::cout << "Processing sulbist: " << sublist_name << "\n";
  Flow::FlowBoundaryFunction* bc = new Flow::FlowBoundaryFunction(mesh_);
  try {
    ProcessMolarFractionList(plist_->sublist(sublist_name), submodel, bc);
  } catch (Errors::Message& msg) {
    Errors::Message m;
    m << "MultiphaseBCFactory: \"" << sublist_name << "\" sublist error: " << msg.what();
    Exceptions::amanzi_throw(m);
  } catch (Teuchos::Exceptions::InvalidParameterType& msg) {
    Errors::Message m;
    m << "MultiphaseBCFactory: \"" << sublist_name << "\" sublist error: not a sublist: " << msg.what();
    Exceptions::amanzi_throw(m);
  }
  return bc;
}


/* ******************************************************************
* Process Dirichet BC (molar fraction), step 2.
* Loop over sublists with typical names "BC 0", "BC 1", etc.
****************************************************************** */
void MultiphaseBCFactory::ProcessMolarFractionList(
    Teuchos::ParameterList& list, std::vector<int>& submodel, Flow::FlowBoundaryFunction* bc) const
{
  // Iterate through the BC specification sublists in the list.
  // All are expected to be sublists of identical structure.
  for (Teuchos::ParameterList::ConstIterator i = list.begin(); i != list.end(); ++i) {
    std::string name = i->first;
    if (list.isSublist(name)) {
      Teuchos::ParameterList& spec = list.sublist(name);
      try {
        ProcessMolarFractionSpec(spec, submodel, bc);
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
* Process Dirichet BC (molar fraction), step 3.
****************************************************************** */
void MultiphaseBCFactory::ProcessMolarFractionSpec(
    Teuchos::ParameterList& list, std::vector<int>& submodel, Flow::FlowBoundaryFunction* bc) const
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
  if (list.isSublist("boundary molar fraction")) {
    f_list = &list.sublist("boundary molar fraction");
  } else {  // "boundary molar fraction" is not a sublist
    m << "parameter \"boundary molar fraction\" is not a sublist";
    Exceptions::amanzi_throw(m);
  }

  // Make the boundary molar fraction function.
  Teuchos::RCP<Amanzi::MultiFunction> f;
  try {
    f = Teuchos::rcp(new Amanzi::MultiFunction(*f_list));
  } catch (Errors::Message& msg) {
    m << "error in sublist \"boundary molar fraction\": " << msg.what();
    Exceptions::amanzi_throw(m);
  }

  // Add this BC specification to the boundary function.
  bc->Define(regions, f, Amanzi::CommonDefs::BOUNDARY_FUNCTION_ACTION_NONE);
}


/* ******************************************************************
* Process Neumann BC (mass flux), step 1.
****************************************************************** */
Flow::FlowBoundaryFunction* MultiphaseBCFactory::CreateMassFlux(std::vector<int>& submodel, int phase) const
{
  std::string sublist_base("mass flux");
  std::string sublist_name;
  std::stringstream ss;
  ss << phase;
  if (phase == 0) {
    sublist_name = sublist_base;
  } else {
    sublist_name = sublist_base + " phase" + ss.str(); 
  }

  Flow::FlowBoundaryFunction* bc = new Flow::FlowBoundaryFunction(mesh_);
  try {
    ProcessMassFluxList(plist_->sublist(sublist_name), submodel, bc);
  } catch (Errors::Message& msg) {
    Errors::Message m;
    m << "MultiphaseBCFactory: \"" << sublist_name << "\" sublist error: " << msg.what();
    Exceptions::amanzi_throw(m);
  } catch (Teuchos::Exceptions::InvalidParameterType& msg) {
    Errors::Message m;
    m << "MultiphaseBCFactory: \"" << sublist_name << "\" sublist error: not a sublist: " << msg.what();
    Exceptions::amanzi_throw(m);
  }
  return bc;
}


/* ******************************************************************
* Process Neumann BC (mass flux), step 1.
****************************************************************** */
Flow::FlowBoundaryFunction* MultiphaseBCFactory::CreateMassFlux(std::vector<int>& submodel) const
{
  std::string sublist_name("mass flux total");

  Flow::FlowBoundaryFunction* bc = new Flow::FlowBoundaryFunction(mesh_);
  try {
    ProcessMassFluxList(plist_->sublist(sublist_name), submodel, bc);
  } catch (Errors::Message& msg) {
    Errors::Message m;
    m << "MultiphaseBCFactory: \"" << sublist_name << "\" sublist error: " << msg.what();
    Exceptions::amanzi_throw(m);
  } catch (Teuchos::Exceptions::InvalidParameterType& msg) {
    Errors::Message m;
    m << "MultiphaseBCFactory: \"" << sublist_name << "\" sublist error: not a sublist: " << msg.what();
    Exceptions::amanzi_throw(m);
  }
  return bc;
}


/* ******************************************************************
* Process Neumann BC (mass flux), step 2.
* Iterate through the BC specification sublists in the list with 
* typical names "BC 0", "BC 1", etc. All are expected to be sublists 
* of identical structure.
****************************************************************** */
void MultiphaseBCFactory::ProcessMassFluxList(
    Teuchos::ParameterList& list, std::vector<int>& submodel, Flow::FlowBoundaryFunction* bc) const
{
  for (Teuchos::ParameterList::ConstIterator i = list.begin(); i != list.end(); ++i) {
    std::string name = i->first;
    if (list.isSublist(name)) {
      Teuchos::ParameterList& spec = list.sublist(name);
      try {
        ProcessMassFluxSpec(spec, submodel, bc);
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
void MultiphaseBCFactory::ProcessMassFluxSpec(
    Teuchos::ParameterList& list, std::vector<int>& submodel, Flow::FlowBoundaryFunction* bc) const
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
  Teuchos::RCP<Amanzi::MultiFunction> f;
  try {
    f = Teuchos::rcp(new Amanzi::MultiFunction(*f_list));
  } catch (Errors::Message& msg) {
    m << "error in sublist \"outward mass flux\": " << msg.what();
    Exceptions::amanzi_throw(m);
  }

  // Add this BC specification to the boundary function.
  bc->Define(regions, f, Amanzi::CommonDefs::BOUNDARY_FUNCTION_ACTION_NONE);

  // Populate submodel flags.
  if (list.get<bool>("rainfall", false))
      PopulateSubmodelFlag(regions, FLOW_BC_SUBMODEL_RAINFALL, submodel);
}


/* ******************************************************************
* Populate submodel flags.
****************************************************************** */
void MultiphaseBCFactory::PopulateSubmodelFlag(
    const std::vector<std::string>& regions, int flag, std::vector<int>& submodel) const
{
  int nregions = regions.size();
  for (int n = 0; n < nregions; n++) {
    AmanziMesh::Entity_ID_List faces;
    mesh_->get_set_entities(regions[n], AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED, &faces);

    int nfaces = faces.size();
    for (int m = 0; m < nfaces; m++) {
      int f = faces[m];
      submodel[f] += flag;
    }
  }
}

}
}
