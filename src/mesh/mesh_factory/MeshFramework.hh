/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: William Perkins
*/

#ifndef _MeshFramework_hh_
#define _MeshFramework_hh_

#include <string>
#include <vector>

namespace Amanzi {
namespace AmanziMesh {

  /// A type to identify available mesh frameworks
  enum Framework {
    Simple = 1,
    MOAB,
    STKMESH,
    MSTK
  };

  /// A type with which an ordered preference list of Framework 
  typedef std::vector<Framework> FrameworkPreference;

  /// Get a name for a given framework
  extern std::string framework_name(const Framework& fw);

  /// Generate the default framework preferences
  extern FrameworkPreference default_preference(void);

  /// Modify framework preferences to get those available
  extern FrameworkPreference available_preference(const FrameworkPreference& p);

} // close namespace AmanziMesh
} // close namespace Amanzi

#endif
