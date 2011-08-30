#ifndef BL_PARMPARSE_HELPERS_H
#define BL_PARMPARSE_HELPERS_H

#include <ParmParse.H>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_Array.hpp"	
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

namespace BoxLib {

  void Initialize_ParmParse(const Teuchos::ParameterList& params);

}; // namespace BoxLib

#endif /*BL_PARMPARSE_HELPERS_H*/
