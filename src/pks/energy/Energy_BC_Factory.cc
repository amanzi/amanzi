/*
  Energy

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Neil Carlson (nnc@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <string>
#include <vector>

#include "errors.hh"
#include "MultiFunction.hh"

#include "Energy_BC_Factory.hh"
#include "EnergyDefs.hh"

namespace Amanzi {
namespace Energy {

/* ******************************************************************
* Constructor
****************************************************************** */
EnergyBCFactory::EnergyBCFactory(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                                 const Teuchos::RCP<Teuchos::ParameterList>& plist) :
   mesh_(mesh),
   plist_(plist) {};


/* ******************************************************************
* Process Dirichet BC (temperature), step 1.
****************************************************************** */
EnergyBoundaryFunction* EnergyBCFactory::CreateTemperature(std::vector<int>& submodel) const
{
  EnergyBoundaryFunction* bc = new EnergyBoundaryFunction(mesh_);
  try {
    ProcessTemperatureList(plist_->sublist("temperature"), submodel, bc);
  } catch (Errors::Message& msg) {
    Errors::Message m;
    m << "EnergyBCFactory: \"temperature\" sublist error: " << msg.what();
    Exceptions::amanzi_throw(m);
  } catch (Teuchos::Exceptions::InvalidParameterType& msg) {
    Errors::Message m;
    m << "EnergyBCFactory: \"temperature\" sublist error: not a sublist: " << msg.what();
    Exceptions::amanzi_throw(m);
  }
  return bc;
}


/* ******************************************************************
* Process Neumann BC (energy flux), step 1.
****************************************************************** */
EnergyBoundaryFunction* EnergyBCFactory::CreateEnergyFlux(std::vector<int>& submodel) const
{
  int ncells = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);

  EnergyBoundaryFunction* bc = new EnergyBoundaryFunction(mesh_);
  try {
    ProcessEnergyFluxList(plist_->sublist("energy flux"), submodel, bc);
  } catch (Errors::Message& msg) {
    Errors::Message m;
    m << "EnergyBCFactory: \"energy flux\" sublist error: " << msg.what();
    Exceptions::amanzi_throw(m);
  } catch (Teuchos::Exceptions::InvalidParameterType& msg) {
    Errors::Message m;
    m << "EnergyBCFactory: \"energy flux\" sublist error: not a sublist: " << msg.what();
    Exceptions::amanzi_throw(m);
  }
  return bc;
}


/* ******************************************************************
* Process Dirichet BC (temperature), step 2.
* Loop over sublists with typical names "BC 0", "BC 1", etc.
****************************************************************** */
void EnergyBCFactory::ProcessTemperatureList(
    Teuchos::ParameterList& list, std::vector<int>& submodel, EnergyBoundaryFunction* bc) const
{
  // Iterate through the BC specification sublists in the list.
  // All are expected to be sublists of identical structure.
  for (Teuchos::ParameterList::ConstIterator i = list.begin(); i != list.end(); ++i) {
    std::string name = i->first;
    if (list.isSublist(name)) {
      Teuchos::ParameterList& spec = list.sublist(name);
      try {
        ProcessTemperatureSpec(spec, submodel, bc);
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
* Process Dirichet BC (temperature), step 3.
****************************************************************** */
void EnergyBCFactory::ProcessTemperatureSpec(
    Teuchos::ParameterList& list, std::vector<int>& submodel, EnergyBoundaryFunction* bc) const
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
  if (list.isSublist("boundary temperature")) {
    f_list = &list.sublist("boundary temperature");
  } else {  // "boundary temperature" is not a sublist
    m << "parameter \"boundary temperature\" is not a sublist";
    Exceptions::amanzi_throw(m);
  }

  // Make the boundary temperature function.
  Teuchos::RCP<Amanzi::MultiFunction> f;
  try {
    f = Teuchos::rcp(new Amanzi::MultiFunction(*f_list));
  } catch (Errors::Message& msg) {
    m << "error in sublist \"boundary temperature\": " << msg.what();
    Exceptions::amanzi_throw(m);
  }

  // Add this BC specification to the boundary function.
  bc->Define(regions, f, CommonDefs::BOUNDARY_FUNCTION_ACTION_NONE);
}


/* ******************************************************************
* Process Neumann BC (energy flux), step 2.
* Iterate through the BC specification sublists in the list with 
* typical names "BC 0", "BC 1", etc. All are expected to be sublists 
* of identical structure.
****************************************************************** */
void EnergyBCFactory::ProcessEnergyFluxList(
    Teuchos::ParameterList& list, std::vector<int>& submodel, EnergyBoundaryFunction* bc) const
{
  for (Teuchos::ParameterList::ConstIterator i = list.begin(); i != list.end(); ++i) {
    std::string name = i->first;
    if (list.isSublist(name)) {
      Teuchos::ParameterList& spec = list.sublist(name);
      try {
        ProcessEnergyFluxSpec(spec, submodel, bc);
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
* Process Neumann BC (energy flux), step 3.
****************************************************************** */
void EnergyBCFactory::ProcessEnergyFluxSpec(
    Teuchos::ParameterList& list, std::vector<int>& submodel, EnergyBoundaryFunction* bc) const
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
  if (list.isSublist("outward energy flux")) {
    f_list = &list.sublist("outward energy flux");
  } else {
    m << "parameter \"outward energy flux\" is not a sublist";
    Exceptions::amanzi_throw(m);
  }

  // Make the outward energy flux function.
  Teuchos::RCP<Amanzi::MultiFunction> f;
  try {
    f = Teuchos::rcp(new Amanzi::MultiFunction(*f_list));
  } catch (Errors::Message& msg) {
    m << "error in sublist \"outward energy flux\": " << msg.what();
    Exceptions::amanzi_throw(m);
  }

  // Add this BC specification to the boundary function.
  bc->Define(regions, f, CommonDefs::BOUNDARY_FUNCTION_ACTION_NONE);
}


/* ******************************************************************
* Populate submodel flags.
****************************************************************** */
void EnergyBCFactory::PopulateSubmodelFlag(
    const std::vector<std::string>& regions, int flag, std::vector<int>& submodel) const
{
  int nregions = regions.size();
  for (int n = 0; n < nregions; n++) {
    AmanziMesh::Entity_ID_List faces;
    std::vector<double> vofs;
    mesh_->get_set_entities(regions[n], AmanziMesh::FACE, AmanziMesh::OWNED, &faces, &vofs);

    int nfaces = faces.size();
    for (int m = 0; m < nfaces; m++) {
      int f = faces[m];
      submodel[f] += flag;
    }
  }
}

}  // namespace Energy
}  // namespace Amanzi
