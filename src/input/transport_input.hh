#ifndef TRANSPORT_INPUT_HPP
#define TRANSPORT_INPUT_HPP

#include <vector>
#include <string>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"

#include "input_base_class.hh"

class TransportInput : public InputBaseClass
{
  public:

    /* Constructors */
    TransportInput();

    /* Methods */


    /* Overridden methods from ParameterListAcceptorDefaultBase */
    void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const & );
    Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

    /* Methods that alter internal parameters */
    void set_cfl(double);
    void set_max_dt(double);
    void set_verbosity(int);



  private:
  
    /* Keywords */
    static const std::string verbosity_level_keyword_;
    static const std::string cfl_keyword_;
    static const std::string enable_internal_tests_keyword_;
    static const std::string internal_tol_keyword_;
    static const std::string max_dt_keyword_;
    static const std::string num_bc_keyword_;
    static const std::string bc_sublist_name_;

    int    verbose_;
    bool   internal_tests_;
    double tests_tolerance_;
    double cfl_;
    double dT_debug_;
    int    num_bc_;
    
    Teuchos::RCP<Teuchos::ParameterList> defineValidTransportParameterList()
        const;






};

        
#endif
