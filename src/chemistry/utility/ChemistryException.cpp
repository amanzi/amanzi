/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/*******************************************************************************
 **
**  Description:
**      Simple Exception class for chemistry errors
**
**  Notes:
**    - see
**       - Eckel, Bruce and Chuck Allison. 2004. Thinking In C++ 
**           Volume 2: Practical Programming. Pearson Prentice Hall. Ch. 1
**       - Stroustrup, Bjarne. 2000. The C++ programing language. Special Edition
**           Addison Wesley. Ch.8.3 and 14
**
**  Example Use:
**    std::ostringstream error_stream;
**    error_stream << "ERROR: MyClass::MyFunction(): \n";
**    error_stream << "ERROR: something has gone wrong....\n";
**    throw ChemistryException(error_stream.str());  
**
**    try {
**      ...
**    }
**    catch (ChemistryException& geochem_error) {
**      geochem_error.PrintErrorMessage();
**      std::cerr << geochem_error.what() << std::endl;
**    }
**
*******************************************************************************/
#include <string>
#include <iostream>

#include "ChemistryException.hpp" 

ChemistryException::ChemistryException() 
    : runtime_error("ERROR"),
      error_status_(ChemistryException::kOkay)
{
    /* end ChemistryException() */
}

ChemistryException::ChemistryException(std::string error_message) 
    : runtime_error(error_message),
      error_status_(ChemistryException::kOkay)
{
    /* end ChemistryException(std::string) */
}

ChemistryException::ChemistryException(std::string error_message, Status error_status) 
    : runtime_error(error_message),
      error_status_(error_status)
{
    /* end ChemistryException(std::string) */
}

ChemistryException::~ChemistryException() throw ()
{
    /* end ~ChemistryException() */
}


void ChemistryException::PrintErrorMessage(void) const
{
    std::cout << what() << std::endl;
    /* end print_message() */
}
