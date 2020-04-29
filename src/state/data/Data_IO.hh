/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! Data_Initializers are functions that know how to initialize data.

/*
  Initialize from parameter list.
*/

#pragma once

#include "Teuchos_ParameterList.hpp"

#include "errors.hh"

namespace Amanzi {

class Visualization;
class Checkpoint;

namespace Data_Initializers {


//
// Null functions for things that don't vis or checkpoint
// ======================================================================
template<class T>
inline void
UserWriteVis(const Visualization& vis,
             const Teuchos::ParameterList& attrs, const T& t) {}

template<class T>
inline void
UserWriteCheckpoint(const Checkpoint& vis,
                    const Teuchos::ParameterList& attrs, const T& t) {}

template<class T>
inline void
UserReadCheckpoint(const Checkpoint& vis,
                   const Teuchos::ParameterList& attrs, T& t) {}


} // namespace Data_Initializers
} // namespace Amanzi

