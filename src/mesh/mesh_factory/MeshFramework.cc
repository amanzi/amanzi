// -------------------------------------------------------------
/**
 * @file   MeshFramework.cc
 * @author William A. Perkins
 * @date Tue May 17 11:49:33 2011
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created March 11, 2011 by William A. Perkins
// Last Change: Tue May 17 11:49:33 2011 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

static const char* SCCS_ID = "$Id$ Battelle PNL";

#include "MeshFramework.hh"
#include "MeshException.hh"
#include "FrameworkTraits.hh"

namespace Amanzi {
namespace AmanziMesh {
  // -------------------------------------------------------------
  // framework_name
  // -------------------------------------------------------------
  std::string
  framework_name(const Framework& fw)
  {
    std::string result("unknown");
    switch (fw) {
    case (Simple):
      result = "Simple";
      break;
    case (STKMESH):
      result = "stk::mesh";
      break;
    case (MOAB):
      result = "MOAB";
      break;
    case (MSTK):
      result = "MSTK";
      break;
    }
    return result;
  }

  // -------------------------------------------------------------
  // default_preference
  // -------------------------------------------------------------
  FrameworkPreference
  default_preference(void)
  {
    FrameworkPreference result;

    // order is important here, it is the order in which the framework
    // is chosen, if there is a choice

    if (framework_available(STKMESH)) result.push_back(STKMESH);
    if (framework_available(Simple)) result.push_back(Simple);
    if (framework_available(MOAB)) result.push_back(MOAB);
    if (framework_available(MSTK)) result.push_back(MSTK);
  
    return result;
  }

  // -------------------------------------------------------------
  // available_preference

  // -------------------------------------------------------------
  /** 
   * This routine removes entries from the request preferences if they
   * are not available.
   * 
   * @param request 
   * 
   * @return new preference list with only available frameworks
   */  
  FrameworkPreference
  available_preference(const FrameworkPreference& request)
  {
    FrameworkPreference result;
    FrameworkPreference defpref(default_preference());

    for (FrameworkPreference::const_iterator i = request.begin();
         i != request.end(); i++) {
      if (std::find(defpref.begin(), defpref.end(), *i) != defpref.end()) {
        result.push_back(*i);
      }
    }
    return result;
  }
      

} // namespace AmanziMesh
} // namespace Amanzi
