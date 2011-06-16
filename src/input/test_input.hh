#ifndef TEST_INPUT_HPP
#define TEST_INPUT_HPP

#include <vector>
#include <string>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"

#include "input_base_class.hh"

/*
   This is a test input class. Designed for language
   interface tests. Do NOT use.
*/
class TestInput : public InputBaseClass
{
  public:

    /* Constructors */
    TestInput();
    TestInput(int myint, double mydbl, bool mybool);

    /* Methods */
    void add_sublist(int x0,double x1, bool x2);


    /* Overridden methods from ParameterListAcceptorDefaultBase */
    void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const & );
    Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;


  private:
  
    /* Keywords */
    static const std::string num_sublist_keyword_;
    static const std::string int_parameter_keyword_;
    static const std::string double_parameter_keyword_;
    static const std::string bool_parameter_keyword_;
    static const std::string sublist_root_name_;

    int    num_sublist_;
    int    Int_;
    double Dbl_ ;
    bool   Bool_;
    
    Teuchos::RCP<Teuchos::ParameterList>
        defineValidTestParameterList() const; 

    std::string generate_sublist_name(int) const;


};

        
#endif
