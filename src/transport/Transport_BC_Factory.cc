/*
This is the transport component of the Amanzi code. 

Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Authors:  Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <vector>

#include "errors.hh"

#include "Transport_BC_Factory.hh"
#include "TransportDefs.hh"

namespace Amanzi {
namespace AmanziTransport {

/* ******************************************************************
* Process Dirichet BC (concentration), step 1.
****************************************************************** */
void TransportBCFactory::CreateConcentration(
    std::vector<Functions::TransportBoundaryFunction*>& bcs, std::vector<std::string>& bcs_tcc_name) const
{
  Errors::Message msg;

#if ALQUIMIA_ENABLED
  if (chem_engine_ != Teuchos::null)
  {
    // We create a set of boundary Functions using the Chemistry Engine.
    if (list_->isSublist("geochemical condition")) {
      Teuchos::ParameterList& clist = list_->get<Teuchos::ParameterList>("geochemical condition");
      for (Teuchos::ParameterList::ConstIterator it = clist.begin(); it != clist.end(); ++it) {
        std::string cond_name = it->first; // name of geochemical condition
        if (clist.isSublist(cond_name)) {
          Teuchos::ParameterList& cond_list = clist.sublist(cond_name);

          // Set up a transport boundary function for this species that uses the Chemistry Engine 
          // to satisfy the given geochemical condition.
          Functions::TransportBoundaryFunction* bc = new Functions::TransportBoundaryFunction(mesh_);
          ProcessGeochemicalCondition(cond_name, cond_list, bc);
          bcs.push_back(bc);
//          bcs_tcc_name.push_back(name); // FIXME: This is not correct!
        }
      }
    }
  }
  else
#endif
  {
    if (list_->isSublist("concentration")) {
      Teuchos::ParameterList& clist = list_->get<Teuchos::ParameterList>("concentration");

      for (Teuchos::ParameterList::ConstIterator it = clist.begin(); it != clist.end(); ++it) {
        std::string name = it->first;
        if (clist.isSublist(name)) {
          Teuchos::ParameterList& bclist = clist.sublist(name);
          for (Teuchos::ParameterList::ConstIterator it1 = bclist.begin(); it1 != bclist.end(); ++it1) {
            std::string specname = it1->first;

            if (bclist.isSublist(specname)) {
              Teuchos::ParameterList& spec = bclist.sublist(specname);
              try {
                Functions::TransportBoundaryFunction* bc = new Functions::TransportBoundaryFunction(mesh_);
                ProcessConcentrationSpec(spec, bc);
                bcs.push_back(bc);
                bcs_tcc_name.push_back(name);
              } catch (Errors::Message& m) {
                msg << "in sublist \"" << specname.c_str() << "\": " << m.what();
                Exceptions::amanzi_throw(msg);
              }
            } else {
              msg << "parameter \"" << specname.c_str() << "\" is not a sublist.\n";
              Exceptions::amanzi_throw(msg);
            }
          }
        } else {
          msg << "parameter \"" << name.c_str() << "\" is not a sublist.\n";
          Exceptions::amanzi_throw(msg);
        }
      }
    } else {
      msg << "Transport PK: BC has no sublist \"concentration\".\n";
      Exceptions::amanzi_throw(msg);  
    }
  }
}


/* ******************************************************************
* Process Dirichet BC (concentration), step 3.
****************************************************************** */
void TransportBCFactory::ProcessConcentrationSpec(
    Teuchos::ParameterList& spec, Functions::TransportBoundaryFunction* bc) const
{
  Errors::Message msg;
  std::vector<std::string> regions;

  if (spec.isParameter("regions")) {
    if (spec.isType<Teuchos::Array<std::string> >("regions")) {
      regions = spec.get<Teuchos::Array<std::string> >("regions").toVector();
    } else {
      msg << "parameter \"regions\" is not of type \"Array string\"";
      Exceptions::amanzi_throw(msg);
    }
  } else {
    msg << "parameter \"regions\" is missing";
    Exceptions::amanzi_throw(msg);
  }

  Teuchos::ParameterList* f_list;
  if (spec.isSublist("boundary concentration")) {
    f_list = &spec.sublist("boundary concentration");
  } else {  // "boundary concentration" is not a sublist
    msg << "parameter \"boundary concentration\" is not a sublist";
    Exceptions::amanzi_throw(msg);
  }

  // Make the boundary pressure function.
  Teuchos::RCP<Amanzi::MultiFunction> f;
  try {
    f = Teuchos::rcp(new Amanzi::MultiFunction(*f_list));
  } catch (Errors::Message& m) {
    msg << "error in sublist \"boundary concentration\": " << m.what();
    Exceptions::amanzi_throw(msg);
  }
  
  // Add this BC specification to the boundary function.
  bc->Define(regions, f); // , Amanzi::BOUNDARY_FUNCTION_ACTION_NONE); // needs to be fixed
}

#ifdef ALQUIMIA_ENABLED
// Process a geochemical condition.
void TransportBCFactory::ProcessGeochemicalCondition(const std::string& cond_name,
                                                     const Teuchos::ParameterList& cond_list, 
                                                     Functions::TransportBoundaryFunction* bc) const
{
  Errors::Message msg;

  // Figure out which regions this geochemical condition applies to.
  std::vector<std::string> regions;
  if (cond_list.isParameter("regions")) {
    if (cond_list.isType<Teuchos::Array<std::string> >("regions")) {
      regions = cond_list.get<Teuchos::Array<std::string> >("regions").toVector();
    } else {
      msg << "parameter \"regions\" is not of type \"Array string\"";
      Exceptions::amanzi_throw(msg);
    }
  } else {
    msg << "parameter \"regions\" is missing";
    Exceptions::amanzi_throw(msg);
  }

  // Construct functions for each of the species and create multi-functions for each. Here, we use a 
  // Geochemical condition "context" that enforces a geochemical condition and broadcasts the resulting 
  // species concentrations to various GeochemicalConcentrationFunctions.
  // FIXME: I don't think we need to do this anymore.
#if 0
  Teuchos::RCP<AmanziChemistry::GeochemicalConditionContext> geochem_context = 
    Teuchos::rcp(new AmanziChemistry::GeochemicalConditionContext(chem_engine_, cond_name));
  std::vector<std::string> speciesNames;
  chem_engine_->GetPrimarySpeciesNames(speciesNames);
  for (size_t i = 0; i < speciesNames.size(); ++i)
  {
    // Create a function representing the geochemical condition acting upon the species.
    Teuchos::RCP<Function> geochem_func = geochem_context->speciesFunction(speciesNames[i]);
    Teuchos::RCP<Amanzi::MultiFunction> f = Teuchos::rcp(new Amanzi::MultiFunction(geochem_func));

    // Add this BC specification to the boundary function.
    bc->Define(regions, f); // , Amanzi::BOUNDARY_FUNCTION_ACTION_NONE); // needs to be fixed
  }
#endif
}
#endif

}  // namespace AmanziTransport
}  // namespace Amanzi
