#include "dbc.hh"

#include <sstream>

namespace {

std::string form_message (std::string const & assertion,
                          std::string const & file,
                          int line)
{
    std::ostringstream message;
    message << "Assertion: " << assertion << " failed in file: " << file << ", line: " << line << std::endl;

    return message.str();

}


}

void abort_program (std::string const & cond, 
                    std::string const & file,  
                    int line)
{
    DBC_assertion assrt (cond, file, line);

#ifdef _THROW
    throw assrt;
#else
    std::cout << assrt.what () << std::endl;
    abort ();
#endif

}

DBC_assertion::DBC_assertion (std::string const & assrt,
                              std::string const & file,
                              int line_number) : assertion (form_message (assrt, file, line_number)){ };


