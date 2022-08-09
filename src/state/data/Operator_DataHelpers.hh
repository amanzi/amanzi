/*
  Copyright 2010-202x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon
*/

// Helpers for putting Operator data structures into State

#pragma once

#include "Teuchos_ParameterList.hpp"

#include "Key.hh"
#include "Visualization.hh"
#include "Checkpoint.hh"

namespace Amanzi {

namespace Operators {
class BCs;
class Op;
class Operator;
}


namespace Helpers {

// ======================================================================
// Specializations for Operators::BCs
// ======================================================================
inline void
UserWriteVis(const Visualization& vis, const Key& fieldname,
             const std::vector<std::string>* subfieldnames,
             const Operators::BCs& bc) {}

inline void
UserWriteCheckpoint(const Checkpoint& chkp, const Key& fieldname,
                    const std::vector<std::string>* subfieldnames,
                    const Operators::BCs& bc) {}

inline
void UserReadCheckpoint(const Checkpoint& chkp, const Key& fieldname,
                   const std::vector<std::string>* subfieldnames,
                   Operators::BCs& bc) {}

inline
bool UserInitialize(Teuchos::ParameterList& plist,
                           Operators::BCs& bc, const Key& fieldname,
                           const std::vector<std::string>* subfieldnames) {
  return true;
}


// ======================================================================
// Specializations for Op
// ======================================================================
inline
void UserWriteVis(const Visualization& vis, const Key& fieldname,
                             const std::vector<std::string>* subfieldnames,
                             const Operators::Op& vec) {}

inline
void UserWriteCheckpoint(const Checkpoint& chkp, const Key& fieldname,
                                    const std::vector<std::string>* subfieldnames,
                                    const Operators::Op& vec) {}

inline
void UserReadCheckpoint(const Checkpoint& chkp, const Key& fieldname,
                                   const std::vector<std::string>* subfieldnames,
                                   Operators::Op& vec) {}

inline
bool UserInitialize(Teuchos::ParameterList& plist,
                               Operators::Op& t, const Key& fieldname,
                               const std::vector<std::string>* subfieldnames) {
  return true;
}


// ======================================================================
// Specializations for Operator
// ======================================================================
inline
void UserWriteVis(const Visualization& vis, const Key& fieldname,
                                   const std::vector<std::string>* subfieldnames,
                                   const Operators::Operator& global_operator) {}

inline
void UserWriteCheckpoint(const Checkpoint& chkp, const Key& fieldname,
                                          const std::vector<std::string>* subfieldnames,
                                          const Operators::Operator& global_operator) {}

inline
void UserReadCheckpoint(const Checkpoint& chkp, const Key& fieldname,
                                         const std::vector<std::string>* subfieldnames,
                                         Operators::Operator& global_operator) {}

inline
bool UserInitialize(Teuchos::ParameterList& plist,
                                     Operators::Operator& t, const Key& fieldname,
                                     const std::vector<std::string>* subfieldnames) {
  return true;
}

inline
void UserAssign(Operators::Operator& dest,
        const Operators::Operator& source) {
  Errors::Message msg("Operators::Operator: assignment operator not supported.");
  Exceptions::amanzi_throw(msg);
}



} // namespace Helpers
} // namespace Amanzi
