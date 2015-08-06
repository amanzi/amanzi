/*
  This is the input component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <sstream>
#include <string>

//TPLs
#include <xercesc/dom/DOM.hpp>
#include <xercesc/util/XMLString.hpp>

// Amanzi's
#include "errors.hh"
#include "exceptions.hh"
#include "dbc.hh"

#include "InputConverterU.hh"
#include "InputConverterU_Defs.hh"

namespace Amanzi {
namespace AmanziInput {

XERCES_CPP_NAMESPACE_USE

/* ******************************************************************
* Create MPC list, version 2, dubbed cycle driver.
****************************************************************** */
Teuchos::ParameterList InputConverterU::TranslateCycleDriver_()
{
  Teuchos::ParameterList out_list;

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    *vo_->os() << "Translating cycle driver" << std::endl;
  }

  // available PKs
  bool transport_pk(false), chemistry_pk(false), flow_pk(false), energy_pk(false);
  std::string flow_model, chemistry_model, energy_model;

  char* tagname;

  DOMNodeList* node_list = doc_->getElementsByTagName(XMLString::transcode("process_kernels"));
  DOMNodeList* children = node_list->item(0)->getChildNodes();

  for (int i = 0; i < children->getLength(); ++i) {
    DOMNode* inode = children->item(i);
    if (inode->getNodeType() != DOMNode::ELEMENT_NODE) continue;
    tagname = XMLString::transcode(inode->getNodeName());

    if (strcmp(tagname, "transport") == 0) {
      char* state = GetAttributeValueC_(static_cast<DOMElement*>(inode), "state");
      transport_pk = (strcmp(state, "on") == 0);
      XMLString::release(&state);
    } else if (strcmp(tagname, "flow") == 0) {
      char* state = GetAttributeValueC_(static_cast<DOMElement*>(inode), "state");
      flow_pk = (strcmp(state, "on") == 0);
      XMLString::release(&state);

      char* model = GetAttributeValueC_(static_cast<DOMElement*>(inode), "model");
      flow_model = TrimString_(model);
      XMLString::release(&model);
    } else if (strcmp(tagname, "chemistry") == 0) {
      char* state = GetAttributeValueC_(static_cast<DOMElement*>(inode), "state");
      flow_pk = (strcmp(state, "on") == 0);
      XMLString::release(&state);

      char* model = GetAttributeValueC_(static_cast<DOMElement*>(inode), "engine");
      chemistry_model = TrimString_(model);
      XMLString::release(&model);
    }
    XMLString::release(&tagname);
  }

  // ???
  /*
  Teuchos::ParameterList pk_tree_list, pk_tree_list_pre;
  double start_time(0.0), switch_time(0.0), end_time(0.0);
  double dt_steady(1.0), dt_tran(1.0);
  int max_cycle_number(-1);

  int model(0), model_pre(0);
  model = chemistry_pk + 2*transport_pk + 4*flow_pk;

  Teuchos::RCP<Teuchos::ParameterList> tim_list = Teuchos::sublist(exe_list, "Time Integration Mode", true);

  if (tim_list->isSublist("Initialize To Steady")) {
    start_time = tim_list->sublist("Initialize To Steady").get<double>("Start");
    switch_time = tim_list->sublist("Initialize To Steady").get<double>("Switch");
    end_time = tim_list->sublist("Initialize To Steady").get<double>("End");
      
    dt_steady = tim_list->sublist("Initialize To Steady").get<double>("Steady Initial Time Step", 1);
    dt_tran = tim_list->sublist("Initialize To Steady").get<double>("Transient Initial Time Step", 1);
    flow_name = "Flow";

    model_pre = 4*flow_pk;
  }
  if (tim_list->isSublist("Steady")) {
    start_time = tim_list->sublist("Steady").get<double>("Start");
    end_time = tim_list->sublist("Steady").get<double>("End");
    dt_steady = tim_list->sublist("Steady").get<double>("Initial Time Step",  1);
    switch_time = start_time;
    dt_tran = dt_steady;
    flow_name = "Flow Steady";
  }
  if (tim_list->isSublist("Transient")) {
    start_time = tim_list->sublist("Transient").get<double>("Start");
    end_time = tim_list->sublist("Transient").get<double>("End");
    switch_time = start_time;
    dt_tran = tim_list->sublist("Transient").get<double>("Initial Time Step", 1);
    flow_name = "Flow";
        
    if (tim_list->sublist("Transient").isParameter("Maximum Cycle Number"))
      max_cycle_number = tim_list->sublist("Transient").get<double>("Maximum Cycle Number");
  }
  if (tim_list->isSublist("Transient with Static Flow")) {
    start_time = tim_list->sublist("Transient with Static Flow").get<double>("Start");
    end_time = tim_list->sublist("Transient with Static Flow").get<double>("End");
    switch_time = start_time;
    dt_tran = tim_list->sublist("Transient with Static Flow").get<double>("Initial Time Step", 1);
    flow_name = "Flow";

    if (tim_list->sublist("Transient with Static Flow").isParameter("Maximum Cycle Number"))
      max_cycle_number = tim_list->sublist("Transient with Static Flow").get<double>("Maximum Cycle Number");
  }      

  int time_pr_id = 0;
  std::string tp_list_name;

  if (model_pre > 0) {
    std::ostringstream ss; ss << time_pr_id;
    tp_list_name = "TP "+ ss.str();

    pk_tree_list_pre.sublist("Flow Steady").set<std::string>("PK type", flow_pk);
    out_list.sublist("time periods").sublist(tp_list_name.data()).sublist("PK Tree") = pk_tree_list_pre;

    out_list.sublist("time periods").sublist(tp_list_name.data()).set<double>("start period time", start_time);
    out_list.sublist("time periods").sublist(tp_list_name.data()).set<double>("end period time", switch_time);
    out_list.sublist("time periods").sublist(tp_list_name.data()).set<int>("maximum cycle number", max_cycle_number);
    out_list.sublist("time periods").sublist(tp_list_name.data()).set<double>("initial time step", dt_steady);
    time_pr_id++;
  }

  switch (model) {
  case 1:
    // Chemistry
    pk_tree_list.sublist("Chemistry").set<std::string>("PK type", chemistry_pk);
    break;
  case 2:
    // Transport
    pk_tree_list.sublist("Transport").set<std::string>("PK type", "transport");;
    break;
  case 3:
    // Reactive Transport
    pk_tree_list.sublist("Reactive Transport").set<std::string>("PK type", "reactive transport");
    pk_tree_list.sublist("Reactive Transport").sublist("Transport").set<std::string>("PK type", "transport");
    pk_tree_list.sublist("Reactive Transport").sublist("Chemistry").set<std::string>("PK type", chemistry_pk);  
    break;
  case 4:
    // Flow
    pk_tree_list.sublist(flow_name).set<std::string>("PK type", flow_pk);    
    break;
  case 5:
    // Flow + Chemistry
    pk_tree_list.sublist("Flow and Chemistry").set<std::string>("PK type", "flow reactive transport");
    pk_tree_list.sublist("Flow and Chemistry").sublist("Chemistry").set<std::string>("PK type", chemistry_pk);
    pk_tree_list.sublist("Flow and Chemistry").sublist("Flow").set<std::string>("PK type", flow_pk);   
    break;
  case 6:
    // Flow + Transport
    pk_tree_list.sublist("Flow and Transport").set<std::string>("PK type", "flow reactive transport");
    pk_tree_list.sublist("Flow and Transport").sublist("Transport").set<std::string>("PK type", "transport");
    pk_tree_list.sublist("Flow and Transport").sublist("Flow").set<std::string>("PK type", flow_pk);
    break;
  case 7:
    // Flow + Reactive Transport
    pk_tree_list.sublist("Flow and Reactive Transport").set<std::string>("PK type", "flow reactive transport");
    pk_tree_list.sublist("Flow and Reactive Transport").sublist("Reactive Transport").set<std::string>("PK type", "reactive transport");
    pk_tree_list.sublist("Flow and Reactive Transport").sublist("Reactive Transport").sublist("Transport").set<std::string>("PK type", "transport");
    pk_tree_list.sublist("Flow and Reactive Transport").sublist("Reactive Transport").sublist("Chemistry").set<std::string>("PK type", chemistry_pk);
    pk_tree_list.sublist("Flow and Reactive Transport").sublist("Flow").set<std::string>("PK type", flow_pk);
    break;
  default:
    Exceptions::amanzi_throw(Errors::Message("This model does not supported by new MPC driver."));
  }

  std::ostringstream ss; ss << time_pr_id;
  tp_list_name = "TP "+ ss.str();

  out_list.sublist("time periods").sublist(tp_list_name.data()).sublist("PK Tree") = pk_tree_list;
  if (transport_pk || chemistry_pk) {
    out_list.set<Teuchos::Array<std::string> >("component names", comp_names_all_);
  }

  Teuchos::ParameterList tpc_list = CreateTimePeriodControlList_(plist);

  out_list.sublist("time periods").sublist(tp_list_name.data()).set<double>("start period time", switch_time);
  out_list.sublist("time periods").sublist(tp_list_name.data()).set<double>("end period time", end_time);
  out_list.sublist("time periods").sublist(tp_list_name.data()).set<int>("maximum cycle number", max_cycle_number);
  out_list.sublist("time periods").sublist(tp_list_name.data()).set<double>("initial time step", dt_tran);
  out_list.sublist("Time Period Control") = tpc_list;

  out_list.sublist("VerboseObject") = CreateVerbosityList_(verbosity_level);

  if (plist->sublist("Execution Control").isSublist("Restart") &&
    plist->sublist("Initial Conditions").isParameter("Init from Checkpoint File")) {
    // this is an error, you can either restart or re-init, but not both
    Exceptions::amanzi_throw(Errors::Message("You can either restart from a checkpoint or initialize from a checkpoint, but not both."));
  }
    
  if (plist->sublist("Execution Control").isSublist("Restart")) {
    out_list.sublist("Restart") = plist->sublist("Execution Control").sublist("Restart");
  }
    
  if (plist->sublist("Initial Conditions").isParameter("Init from Checkpoint File")) {
    Teuchos::ParameterList& rest_list = out_list.sublist("Restart");
    std::string file = plist->sublist("Initial Conditions").get<std::string>("Init from Checkpoint File");
    rest_list.set<std::string>("File Name", file);
    rest_list.set<bool>("initialize from checkpoint data file and do not restart", true);
  }
  */

  return out_list;
}

}  // namespace AmanziInput
}  // namespace Amanzi
