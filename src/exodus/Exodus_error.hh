#ifndef _EXODUS_ERROR_HH_
#define _EXODUS_ERROR_HH_

namespace ExodusII
{


class ExodusError : public std::exception
{
    char const * what () const throw ()  { return "Unknown Exodus Error."; } 

public:
    
    int retval_;
    ExodusError (int retval) : retval_ (retval) { }
};


}

#endif
