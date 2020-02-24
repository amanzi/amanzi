/*
  Flow PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_VWCONTENT_EVALUATOR_FACTORY_HH_
#define AMANZI_VWCONTENT_EVALUATOR_FACTORY_HH_

#include <string>

#include "Teuchos_ParameterList.hpp"

#include "VWContentEvaluator.hh"
#include "VWContentEvaluator_ConstDensity.hh"

namespace Amanzi {
namespace Flow {

class VWContentEvaluatorFactory {
 public:
  Teuchos::RCP<VWContentEvaluator> Create(
      const std::string& name, Teuchos::ParameterList& vwc_list)
  {
    if (name == "constant density") {
      Teuchos::RCP<VWContentEvaluator_ConstDensity> 
          vwc = Teuchos::rcp(new VWContentEvaluator_ConstDensity(vwc_list));
      vwc->Init_();
      return vwc;
    } else {
      if (name == "two-phase") vwc_list.set<bool>("water vapor", true);
      Teuchos::RCP<VWContentEvaluator>
          vwc = Teuchos::rcp(new VWContentEvaluator(vwc_list));
      vwc->Init_();
      return vwc;
    }
  }
};

}  // namespace Flow
}  // namespace Amanzi

#endif
