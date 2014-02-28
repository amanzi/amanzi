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

#include "TransportBCFactory.hh"
#include "TransportDefs.hh"

namespace Amanzi {
namespace AmanziTransport {

/* ******************************************************************
* Process Dirichet BCs.
****************************************************************** */
void TransportBCFactory::CreateConcentration(std::vector<TransportBoundaryFunction*>& bcs) const
{
  Errors::Message msg;

  if (list_->isSublist("concentration")) {
    ProcessTracerList(bcs);
  } else if (list_->isSublist("geochemical conditions")) {
    ProcessGeochemicalConditionList(bcs);
  } else {
    msg << "Transport PK: BC sublist has not been recognized\n";
    Exceptions::amanzi_throw(msg);  
  }
}


/* ******************************************************************
* Process Dirichet BC (concentration), step 1.
* **************************************************************** */
void TransportBCFactory::ProcessTracerList(std::vector<TransportBoundaryFunction*>& bcs) const
{
  Errors::Message msg;
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
            TransportBoundaryFunction_Tracer* bc = new TransportBoundaryFunction_Tracer(mesh_);
            ProcessTracerSpec(spec, bc);
            bc->tcc_names().push_back(name);

            TransportBoundaryFunction* bc_base = bc;
            bcs.push_back(bc_base);
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
}


/* ******************************************************************
* Process Dirichet BC (concentration), step 3.
****************************************************************** */
void TransportBCFactory::ProcessTracerSpec(
    Teuchos::ParameterList& spec, TransportBoundaryFunction_Tracer* bc) const
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


/* ******************************************************************
* Process Dirichet BC (concentration), step 1.
* **************************************************************** */
void TransportBCFactory::ProcessGeochemicalConditionList(std::vector<TransportBoundaryFunction*>& bcs) const
{
#ifdef ALQUIMIA_ENABLED
  Errors::Message msg;
  if (!list_->isSublist("geochemical conditions"))
  {
    msg << "  No 'geochemical conditions' list was found in 'Transport->boundary conditions'!\n";
    Exceptions::amanzi_throw(msg);
  }
  Teuchos::ParameterList& clist = list_->get<Teuchos::ParameterList>("geochemical conditions");

  for (Teuchos::ParameterList::ConstIterator it = clist.begin(); it != clist.end(); ++it) {
    std::string cond_name = it->first;
    if (clist.isSublist(cond_name)) {
      Teuchos::ParameterList& spec = clist.sublist(cond_name);
      try {
        TransportBoundaryFunction_Alquimia* bc = 
          new TransportBoundaryFunction_Alquimia(cond_name, mesh_, chem_state_, chem_engine_);
        ProcessGeochemicalConditionSpec(spec, bc);

        TransportBoundaryFunction* bc_base = bc;
        bcs.push_back(bc_base);
      } catch (Errors::Message& m) {
        msg << "in sublist \"" << cond_name << "\": " << m.what();
        Exceptions::amanzi_throw(msg);
      }
    } else {
      msg << "parameter \"" << cond_name << "\" is not a sublist.\n";
      Exceptions::amanzi_throw(msg);
    }
  }
#endif
}

// Process a geochemical condition.
void TransportBCFactory::ProcessGeochemicalConditionSpec(Teuchos::ParameterList& spec, 
                                                         TransportBoundaryFunction_Alquimia* bc) const
{
#ifdef ALQUIMIA_ENABLED
  Errors::Message msg;
  std::vector<std::string> regions;

  if (spec.isParameter("regions")) {
    if (spec.isType<Teuchos::Array<std::string> >("regions")) {
      regions = spec.get<Teuchos::Array<std::string> >("regions").toVector();
      bc->Define(regions);
    } else {
      msg << "parameter \"regions\" is not of type \"Array string\"";
      Exceptions::amanzi_throw(msg);
    }
  } else {
    msg << "parameter \"regions\" is missing";
    Exceptions::amanzi_throw(msg);
  }
#else
  msg << "Internal error: Alquimia is not enabled!";
  Exceptions::amanzi_throw(msg);
#endif
}

}  // namespace AmanziTransport
}  // namespace Amanzi
