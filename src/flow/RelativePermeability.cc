/*
This is the flow component of the Amanzi code. 

Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "errors.hh"

#include "WaterRetentionModel.hh"
#include "WRM_vanGenuchten.hh"
#include "WRM_BrooksCorey.hh"
#include "WRM_fake.hh"

#include "RelativePermeability.hh"
#include "Flow_constants.hh"

namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* List to WRM models has to be provided.
****************************************************************** */
void RelativePermeability::Init()
{
  Errors::Message msg;

  int nblocks = 0;  // Find out how many WRM entries there are.
  for (Teuchos::ParameterList::ConstIterator i = list_.begin(); i != list_.end(); i++) {
    if (list_.isSublist(list_.name(i))) {
      nblocks++;
    } else {
      msg << "Flow PK: WRM contains an entry that is not a sublist.\n";
      Exceptions::amanzi_throw(msg);
    }
  }

  WRM_.resize(nblocks);

  int iblock = 0;
  for (Teuchos::ParameterList::ConstIterator i = list_.begin(); i != list_.end(); i++) {
    if (list_.isSublist(list_.name(i))) {
      Teuchos::ParameterList& wrm_list = list_.sublist(list_.name(i));

      std::string region;
      if (wrm_list.isParameter("region")) {
        region = wrm_list.get<std::string>("region");  // associated mesh block
      } else {
        msg << "Flow PK: WMR sublist \"" << list_.name(i).c_str() << "\" has no parameter \"region\".\n";
        Exceptions::amanzi_throw(msg);
      }

      if (wrm_list.get<string>("water retention model") == "van Genuchten") {
        double m = wrm_list.get<double>("van Genuchten m", FLOW_WRM_EXCEPTION);
        double alpha = wrm_list.get<double>("van Genuchten alpha", FLOW_WRM_EXCEPTION);
        double l = wrm_list.get<double>("van Genuchten l", FLOW_WRM_VANGENUCHTEN_L);
        double sr = wrm_list.get<double>("residual saturation", FLOW_WRM_EXCEPTION);
        double pc0 = wrm_list.get<double>("regularization interval", FLOW_WRM_REGULARIZATION_INTERVAL);
        std::string krel_function = wrm_list.get<std::string>("relative permeability model", "Mualem");

        VerifyWRMparameters(m, alpha, sr, pc0);
        VerifyStringMualemBurdine(krel_function);

        WRM_[iblock] = Teuchos::rcp(new WRM_vanGenuchten(region, m, l, alpha, sr, krel_function, pc0));

      } else if (wrm_list.get<string>("water retention model") == "Brooks Corey") {
        double lambda = wrm_list.get<double>("Brooks Corey lambda", FLOW_WRM_EXCEPTION);
        double alpha = wrm_list.get<double>("Brooks Corey alpha", FLOW_WRM_EXCEPTION);
        double l = wrm_list.get<double>("Brooks Corey l", FLOW_WRM_BROOKS_COREY_L);
        double sr = wrm_list.get<double>("residual saturation", FLOW_WRM_EXCEPTION);
        double pc0 = wrm_list.get<double>("regularization interval", FLOW_WRM_REGULARIZATION_INTERVAL);
        std::string krel_function = wrm_list.get<std::string>("relative permeability model", "Mualem");

        VerifyWRMparameters(lambda, alpha, sr, pc0);
        VerifyStringMualemBurdine(krel_function);

        WRM_[iblock] = Teuchos::rcp(new WRM_BrooksCorey(region, lambda, l, alpha, sr, krel_function, pc0));

      } else if (wrm_list.get<string>("water retention model") == "fake") {
        WRM_[iblock] = Teuchos::rcp(new WRM_fake(region));

      } else {
        msg << "Flow PK: unknown water retention model.";
        Exceptions::amanzi_throw(msg);
      }
      iblock++;
    }
  }
}


/* ****************************************************************
* Verify string for the relative permeability model.
**************************************************************** */
void RelativePermeability::VerifyWRMparameters(double m, double alpha, double sr, double pc0)
{
  Errors::Message msg;
  if (m < 0.0 || alpha < 0.0 || sr < 0.0 || pc0 < 0.0) {
    msg << "Flow PK: Negative/undefined parameter in a water retention model.";
    Exceptions::amanzi_throw(msg);    
  }
  if (sr > 1.0) {
    msg << "Flow PK: residual saturation is greater than 1.";
    Exceptions::amanzi_throw(msg);    
  }
}


/* ****************************************************************
* Verify string for the relative permeability model.
**************************************************************** */
void RelativePermeability::VerifyStringMualemBurdine(const std::string name)
{
  Errors::Message msg;
  if (name != "Mualem" && name != "Burdine") {
    msg << "Flow PK: supported relative permeability models are Mualem and Burdine.";
    Exceptions::amanzi_throw(msg);
  }
}

}  // namespace AmanziFlow
}  // namespace Amanzi

