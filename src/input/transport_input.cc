
#include <string>
#include <sstream>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "transport_input.hh"

typedef Teuchos::ParameterList::PrintOptions PLPrintOptions;

const std::string TransportInput::verbosity_level_keyword_       = "verbosity level";
const std::string TransportInput::cfl_keyword_                   = "CFL";
const std::string TransportInput::enable_internal_tests_keyword_ = "enable internal tests";
const std::string TransportInput::internal_tol_keyword_          = "internal tests tolerance";
const std::string TransportInput::max_dt_keyword_                = "maximal time step";
const std::string TransportInput::num_bc_keyword_                = "number of BCs";
const std::string TransportInput::bc_sublist_name_               = "Transport BCs";



/* Default constructor */
TransportInput::TransportInput()
  : InputBaseClass("Transport")
{

  /* Default parameter values */
  verbose_         = 0;
  internal_tests_  = false;
  tests_tolerance_ = 1.0e-6;
  cfl_             = 1.0;
  dT_debug_        = 1.0e+99; 
  num_bc_          = 0;  

  /* Set the private parameter list */
  Teuchos::RCP<Teuchos::ParameterList> plist =
      defineValidTransportParameterList();

  setMyParamList(plist);

}

/* Override setParameterList from ParameterListAcceptor */
void TransportInput::setParameterList(Teuchos::RCP<Teuchos::ParameterList> const
                                    &plist)
{
  //TEST_FOR_EXPECT(is_null(plist));

  /* Validate the incoming list */
  Teuchos::RCP<const Teuchos::ParameterList> plist_valid = getValidParameters();
  plist->validateParametersAndSetDefaults(*plist_valid,0);

  /* Set the private parameter list from AcceptorDefaultBase */
  setMyParamList(plist);

  /* Set private variables */
  verbose_         = plist->get<int>(verbosity_level_keyword_);
  internal_tests_  = plist->get<bool>(enable_internal_tests_keyword_);
  tests_tolerance_ = plist->get<double>(internal_tol_keyword_);
  cfl_             = plist->get<double>(cfl_keyword_);
  dT_debug_        = plist->get<double>(max_dt_keyword_); 
  num_bc_          = plist->get<int>(num_bc_keyword_);  
  

}

Teuchos::RCP<const Teuchos::ParameterList>
TransportInput::getValidParameters() const
{
  static Teuchos::RCP<const Teuchos::ParameterList> validParams;

  if ( is_null(validParams) ) {

    validParams = defineValidTransportParameterList();

  }

  return validParams;
}

void TransportInput::set_cfl(double value)
{
           
  Teuchos::RCP<Teuchos::ParameterList> plist = getNonconstParameterList();

  plist->set<double>(cfl_keyword_,value);
  setParameterList(plist);
}

void TransportInput::set_max_dt(double value)
{
  Teuchos::RCP<Teuchos::ParameterList> plist = getNonconstParameterList();

  plist->set<double>(max_dt_keyword_,value);
  setParameterList(plist);
}

void TransportInput::set_verbosity(int value)
{
  Teuchos::RCP<Teuchos::ParameterList> plist = getNonconstParameterList();

  plist->set<int>(verbosity_level_keyword_,value);
  setParameterList(plist);
}




/* Private */

Teuchos::RCP<Teuchos::ParameterList> 
TransportInput::defineValidTransportParameterList() const
{
  Teuchos::RCP<Teuchos::ParameterList> 
      plist = Teuchos::rcp(new Teuchos::ParameterList(paramListName_));
  Teuchos::ParameterEntry entry;

  /* CFL */
  entry.setValue(cfl_,
                 true,
                 "The Courant number",
                 this->DoubleOnlyValidator_);
  plist->setEntry(cfl_keyword_,entry);

  /* Internal tests */
  plist->set<bool>(enable_internal_tests_keyword_,internal_tests_);

  entry.setValue(tests_tolerance_,
                 true,
                 "Internal test tolerance??? Need more information",
                 this->DoubleOnlyValidator_);
  plist->setEntry(internal_tol_keyword_, entry);

  /* Verbosity */
  entry.setValue(verbose_,
                 true,
                 "Print out useful diagnostic information",
                 this->IntOnlyValidator_);
  plist->setEntry(verbosity_level_keyword_, entry);

  /* Maximum time step */
  entry.setValue(dT_debug_,
                 true,
                 "Maximum time-step for the transport",
                 this->DoubleOnlyValidator_);
  plist->setEntry(max_dt_keyword_, entry);

  /* Number of boundary conditions */
  entry.setValue(num_bc_,
                 true,
                 "Number of transport boundary conditions",
                 this->IntOnlyValidator_);
  plist->setEntry(num_bc_keyword_, entry);

  plist->sublist(bc_sublist_name_);

  return plist;
}











