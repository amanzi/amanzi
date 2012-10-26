
#include <string>
#include <sstream>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "test_input.hh"

typedef Teuchos::ParameterList::PrintOptions PLPrintOptions;


const std::string TestInput::num_sublist_keyword_      = "Number of sublists";
const std::string TestInput::int_parameter_keyword_    = "Integer Parameter";
const std::string TestInput::double_parameter_keyword_ = "Double Parameter";
const std::string TestInput::bool_parameter_keyword_   = "Boolean Parameter";
const std::string TestInput::sublist_root_name_        = "Sublists"; 

/* Default constructor */
TestInput::TestInput()
  : InputBaseClass("Test Parameters")
{

  /* Default parameter values */
  Int_  = 0;
  Dbl_  = 0.0;
  Bool_ = false; 

  num_sublist_ = 0;   /* Number of block regions */

  /* Set the private parameter list */
  Teuchos::RCP<Teuchos::ParameterList> plist =
      defineValidTestParameterList();

  setMyParamList(plist);

}

/* Constructor with all parameters set */
TestInput::TestInput(int x0, double x1, bool x2) 
  : InputBaseClass("Test Parameters")
{
  Int_  = x0;
  Dbl_  = x1;
  Bool_ = x2;

  num_sublist_ = 0;

  /* Set the private parameter list */
  Teuchos::RCP<Teuchos::ParameterList> plist =
      defineValidTestParameterList();

  setMyParamList(plist);

}

/* Override setParameterList from ParameterListAcceptor */
void TestInput::setParameterList(Teuchos::RCP<Teuchos::ParameterList> const
                                    &plist)
{
  TEST_FOR_EXCEPTION(is_null(plist),std::exception,"Parameter list is null");

  /* Validate the incoming list */
  Teuchos::RCP<const Teuchos::ParameterList> plist_valid = getValidParameters();
  plist->validateParametersAndSetDefaults(*plist_valid,0);

  /* Set the private parameter list from AcceptorDefaultBase */
  setMyParamList(plist);

  /* Set private variables */
  Int_  = plist->get<int>(int_parameter_keyword_);
  Dbl_  = plist->get<double>(double_parameter_keyword_);
  Bool_ = plist->get<bool>(bool_parameter_keyword_);

  num_sublist_ = plist->get<int>(num_sublist_keyword_);


}

Teuchos::RCP<const Teuchos::ParameterList>
TestInput::getValidParameters() const
{
  static Teuchos::RCP<const Teuchos::ParameterList> validParams;

  if ( is_null(validParams) ) {

    validParams = defineValidTestParameterList();

  }

  return validParams;
}

std::string TestInput::generate_sublist_name(int label) const
{
  std::stringstream ss;
  ss << "Sublist " << label;
  return ss.str();
}

void TestInput::add_sublist(int x0, double x1, bool x2)
{
  /* Grab the internal plist ... will change */
  Teuchos::RCP<Teuchos::ParameterList> plist = getNonconstParameterList();

  /* Update the number of sublists */
  std::string  new_sublist_name = generate_sublist_name(num_sublist_);
  num_sublist_++;
  plist->set<int>(num_sublist_keyword_,num_sublist_);

  /* Now add the new sublist */
  Teuchos::ParameterList& new_list =
      plist->sublist(sublist_root_name_,true).sublist(new_sublist_name);

  new_list.set<int>(int_parameter_keyword_,x0);
  new_list.set<double>(double_parameter_keyword_,x1);
  new_list.set<bool>(bool_parameter_keyword_,x2);

  setParameterList(plist);

}

/* Private */

Teuchos::RCP<Teuchos::ParameterList> 
TestInput::defineValidTestParameterList() const
{
  Teuchos::RCP<Teuchos::ParameterList> 
      plist = Teuchos::rcp(new Teuchos::ParameterList(paramListName_));
  Teuchos::ParameterEntry entry;

  /* Int parameter */
  entry.setValue(Int_,
                 true,
                 "",
                 this->IntOnlyValidator_);
  plist->setEntry(int_parameter_keyword_,entry);

  entry.setValue(Dbl_,
                 true,
                 "",
                 this->DoubleOnlyValidator_);
  plist->setEntry(double_parameter_keyword_,entry);

  /* bool Parameter */
  plist->set<bool>(bool_parameter_keyword_,Bool_);

  /* Number of sublists */
  entry.setValue(num_sublist_,
                 true,
                 "Number of sublists",
                 this->IntOnlyValidator_);
  plist->setEntry(num_sublist_keyword_,entry);
  plist->sublist(sublist_root_name_);


  return plist;
}
