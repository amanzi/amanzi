#ifndef _EXODUS_ERROR_HH_
#define _EXODUS_ERROR_HH_

#include <sstream>



namespace ExodusII
{


class ExodusError : public std::exception
{
    inline char const * what () const throw ();

public:
    
    int retval_;
    ExodusError (int retval) : retval_ (retval) { }
};

char const * ExodusError::what () const throw ()
{
    std::stringstream ss;
    ss << "Unknown Exodus Error (";
    ss << retval_ << ")" <<  std::endl;
    return ss.str ().c_str ();
}

}

#endif
