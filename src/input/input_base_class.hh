#ifndef INPUT_IFACE_HH
#define INPUT_IFACE_HH

#include <string>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"

/* Google standard macro to disallow copy and = assignments */
#define DISALLOW_COPY_AND_ASSIGN(TypeName) \
  TypeName(const TypeName&);               \
  void operator=(const TypeName&)
/*
 *
 * Base Class to handle input data
 * This is a sub-class of ParameterListAcceptorDefaultBase. Classes
 * built from this base class must implement virtual functions
 *  setParametersList and getValidParameters 
 *
 */
class InputBaseClass : public Teuchos::ParameterListAcceptorDefaultBase
{

 

  public:

    /* Constructor */
    InputBaseClass(const std::string );
   

    /* Return input information as an XML string */
    std::string getXmlString() const;

  protected:
    /*
    *  Methods to set parameters with validators
    *  Made these methods protected because they alter the
    *  private Parameter List
    */
    void setValidIntParameter(std::string const & keyword, 
                              int value,
                              std::string const & doc_string);

    void setValidDoubleParameter(std::string const & keyword, 
                                 double value,
                                 std::string const & doc_string);

     

  
  protected:

    std::string paramListName_;

    /*
    *  Validators
    *  When Trilinos 10.7 is available, a validator factory will be 
    *  available  and we will add useful validators here. For now
    *  only, validators that restrict type int and double
    *  are provided 
    */
    Teuchos::RCP<Teuchos::AnyNumberParameterEntryValidator>
                     IntOnlyValidator_;
    Teuchos::RCP<Teuchos::AnyNumberParameterEntryValidator>
                     DoubleOnlyValidator_;

  private:
    DISALLOW_COPY_AND_ASSIGN(InputBaseClass);


};

#endif
