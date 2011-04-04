#ifndef _BDF2_PLISTVALIDATOR_HPP_
#define _BDF2_PLISTVALIDATOR_HPP_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterEntryValidator.hpp"
#include "Teuchos_ParameterEntry.hpp"
#include "Teuchos_Array.hpp"

namespace BDF2 {

  // this is a concrete implementation of the Teuchos::ParamterEntryValidator
  // for the purpose of validating the parameter list that is used in  
  // BDF2::Dae
  
  class PListValidator : public Teuchos::ParameterEntryValidator 
  {
  public: 
    ~PListValidator() {};

    PListValidator() {};

    void printDoc(std::string const& docString, std::ostream& out) const;

    Teuchos::RCP<const Teuchos::Array<std::string> > validStringValues() const;

    void validate(Teuchos::ParameterEntry  const& entry, std::string const& paramName, 
		  std::string const& sublistName) const;
  };
}

#endif // _BDF2_PLISTVALIDATOR_HPP_


