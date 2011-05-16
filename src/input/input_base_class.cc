#include <string>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterEntry.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "input_base_class.hh"

InputBaseClass::InputBaseClass(const std::string list_name)
{
  
  typedef Teuchos::AnyNumberParameterEntryValidator::AcceptedTypes AcceptedTypes;

  /* Initialize the list name */
  paramListName_ = list_name;
 
  /* Initialize the validators */
  this->IntOnlyValidator_ =
    Teuchos::rcp( new Teuchos::AnyNumberParameterEntryValidator(
                  Teuchos::AnyNumberParameterEntryValidator::PREFER_INT,
                  AcceptedTypes(false).allowInt(true).allowDouble(false).allowString(false)
                ));
  this->DoubleOnlyValidator_ =
    Teuchos::rcp( new Teuchos::AnyNumberParameterEntryValidator(
                  Teuchos::AnyNumberParameterEntryValidator::PREFER_DOUBLE,
                  AcceptedTypes(false).allowInt(false).allowDouble(true).allowString(false)
                ));
  
  

}

std::string InputBaseClass::getXmlString() const
{
  std::ostringstream ostream(std::ostringstream::in);
  Teuchos::RCP<const Teuchos::ParameterList> plist = getParameterList();
  Teuchos::writeParameterListToXmlOStream(*plist,ostream);
  return ostream.str();
} 

void InputBaseClass::setValidIntParameter(std::string const &keyword,
                                          int value,
                                          std::string const &doc_string)
{
   Teuchos::RCP<Teuchos::ParameterList> plist = getMyNonconstParamList();

   Teuchos::setIntParameter(keyword,
                            value,
                            doc_string,
                            &*plist);
}

void InputBaseClass::setValidDoubleParameter(std::string const &keyword,
                                             double value,
                                             std::string const &doc_string)
{
   Teuchos::RCP<Teuchos::ParameterList> plist = getMyNonconstParamList();

   Teuchos::setDoubleParameter(keyword,
                               value,
                               doc_string,
                               &*plist);
}
   






