// -------------------------------------------------------------
/**
 * @file   MeshFramework.cc
 * @author William A. Perkins
 * @date Sun Mar 13 20:02:17 2011
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created March 11, 2011 by William A. Perkins
// Last Change: Sun Mar 13 20:02:17 2011 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

static const char* SCCS_ID = "$Id$ Battelle PNL";

#include "MeshFramework.hh"
#include "FrameworkTraits.hh"

// -------------------------------------------------------------
// framework_name
// -------------------------------------------------------------
std::string
framework_name(const Framework& fw)
{
  std::string result("unknown")
  switch (fw) {
  case (Simple):
    result = "Simple";
    break;
  case (STK):
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
  FrameworkPreference pref = {
    Simple, STK, MOAB, MSTK
  };


  
  
  pref.push_back(Simple);

#ifdef HAVE_STK_MESH
  pref.push_back(STK);
#endif

#ifdef HAVE_MOAB_MESH
  pref.push_back(MOAB);
#endif

#ifdef HAVE_MSTK_MESH
  pref.push_back(MSTK);
#endif

  return pref;
}
