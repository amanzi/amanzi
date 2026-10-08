/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
*/

#ifndef AMANZI_PK_MPC_LSCHEME_HELPERS_HH_
#define AMANZI_PK_MPC_LSCHEME_HELPERS_HH_

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "Teuchos_ParameterList.hpp"

#include "Data_Helpers.hh"
#include "Checkpoint.hh"
#include "Key.hh"
#include "Visualization.hh"

namespace Amanzi {

struct LSchemeDataPK {
  double last_step_increment;
  double last_step_residual;

  double safety_factor;
  double safety_factor_min;
  double safety_factor_max;
  int patience = 0;

  double seq_error[3];
  int ns_itrs[3];
  int num_itrs = 0;

  // reset, typically at the begging of time step
  void reset();

  void shift();

  // modify safety factor
  double update();

  // i/o
  void print(std::ostream& os) const;
};
 
typedef std::map<Key, LSchemeDataPK> LSchemeData;


namespace Helpers {

template<>
inline void WriteVis<LSchemeData>(const Visualization& vis,
                                  const Key& fieldname,
                                  const std::vector<std::string>* subfieldnames,
                                  const LSchemeData& t) {};

template<>
inline void WriteCheckpoint<LSchemeData>(const Checkpoint& chkp,
                                         const Key& fieldname,
                                         const std::vector<std::string>* subfieldnames,
                                         const LSchemeData& t) {};

template<>
inline bool ReadCheckpoint<LSchemeData>(const Checkpoint& chkp,
                                        const Key& fieldname,
                                        const std::vector<std::string>* subfieldnames,
                                        LSchemeData& t) { return true; }

template<>
inline bool Initialize<LSchemeData>(Teuchos::ParameterList& plist,
                                    LSchemeData& t,
                                    const Key& fieldname,
                                    const std::vector<std::string>* subfieldnames) { return true; }
} // namespace Helpers
} // namespace Amanzi

#endif
