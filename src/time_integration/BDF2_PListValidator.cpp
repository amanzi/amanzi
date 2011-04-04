#include "BDF2_PListValidator.hpp"

#include "Teuchos_Tuple.hpp"

#include "errors.hh"
#include "exceptions.hh"

#include <ostream>

namespace BDF2 {

  void 
  PListValidator::printDoc(std::string const& docString, std::ostream& out) const 
  {
    out << docString << std::endl;
  }

  
  Teuchos::RCP<const Teuchos::Array<std::string> >
  PListValidator::validStringValues() const
  {
    Teuchos::RCP<Teuchos::Array<std::string> > valid_strings; 
    valid_strings = Teuchos::rcp(new Teuchos::Array<std::string> (5) );
    *valid_strings = Teuchos::tuple<std::string>(
      "Nonlinear solver max iterations",
      "Nonlinear solver tolerance",
      "NKA max vectors",
      "NKA drop tolerance",
      "Verbose");

    return valid_strings;
  }

  void 
  PListValidator::validate(Teuchos::ParameterEntry  const& entry, std::string const& paramName, 
			   std::string const& sublistName) const
  {

    // boilerplate error intro
    std::ostringstream oss;
    
    oss << "\nThe BDF2 parameter\n";
    oss << "  ";
    oss << sublistName;
    oss << "\n    ";
    
    // check the parameters of type int 
    // since both of them must be positive 
    // we can combine the check
    
    if ( paramName == "Nonlinear solver max iterations" || 
	 paramName == "NKA max vectors" ) 
      {
	// validate the type of the entry
	if ( ! entry.isType<int>() ) 
	  {
	    oss << paramName;
	    oss << "\n";
	    oss << "has the wrong type, it must be of type int.\n";
	    
	    Errors::Message msg(oss.str().c_str());
	    Exceptions::amanzi_throw(msg);
	  }
	else
	  {
	    int value = entry.getValue(&value);
	    if ( value <= 0) 
	      {
		oss << paramName << " = " << value;
		oss << "\n";
		oss << "must be positive\n";		
		
		Errors::Message msg(oss.str().c_str());
		Exceptions::amanzi_throw(msg);
	      }
	  }
	
      }
  
    // check the parameters of type double 
    // since both of them must be positive 
    // we can combine the check
    
    if ( paramName == "Nonlinear solver tolerance" ||
	 paramName == "NKA drop tolerance" ) 
      {
	// validate the type of the entry
	if ( ! entry.isType<double>() ) 
	  {
	    oss << paramName;
	    oss << "\n";
	    oss << "has the wrong type, it must be of type double.\n";
	    
	    Errors::Message msg(oss.str().c_str());
	    Exceptions::amanzi_throw(msg);
	  }
	else
	  {
	    double value = entry.getValue(&value);
	    if ( value <= 0.0) 
	      {
		oss << paramName << " = " << value;
		oss << "\n";
		oss << "must be positive\n";		
		
		Errors::Message msg(oss.str().c_str());
		Exceptions::amanzi_throw(msg);
	      }
	    
	    if ( paramName == "Nonlinear solver tolerance" )
	      {
		double value = entry.getValue(&value);
		if ( value > 1.0)
		  {
		    oss << paramName << " = " << value;
		    oss << "\n";
		    oss << "must smaller than 1.0.\n";		
		    
		    Errors::Message msg(oss.str().c_str());
		    Exceptions::amanzi_throw(msg);		    
		    
		  }
	      }
	  }
      }
    
    
    // check the type of the bool entry
    
    if ( paramName == "Verbose" )
      {
	if ( ! entry.isType<bool>() ) 
	  {
	    oss << paramName;
	    oss << "\n";
	    oss << "has the wrong type, it must be of type bool.\n";
	    
	    Errors::Message msg(oss.str().c_str());
	    Exceptions::amanzi_throw(msg);
	  }	
      }
  }
}
