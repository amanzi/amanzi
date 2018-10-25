/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: William Perkins
*/

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

    if (framework_available(MSTK)) result.push_back(MSTK);
    if (framework_available(STKMESH)) result.push_back(STKMESH);
    if (framework_available(MOAB)) result.push_back(MOAB);
    if (framework_available(Simple)) result.push_back(Simple);
  
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
