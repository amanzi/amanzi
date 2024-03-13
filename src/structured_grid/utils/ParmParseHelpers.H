/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#ifndef BL_PARMPARSE_HELPERS_H
#define BL_PARMPARSE_HELPERS_H

#include <ParmParse.H>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

namespace BoxLib {

  void Initialize_ParmParse(Teuchos::ParameterList& params);

}; // namespace BoxLib

#endif /*BL_PARMPARSE_HELPERS_H*/
