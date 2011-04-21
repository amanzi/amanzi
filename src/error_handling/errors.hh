#ifndef _ERRORS_H_
#define _ERRORS_H_

#include "exceptions.hh"

#include <sstream>

namespace Errors
{

class Message : public Exceptions::Amanzi_exception
{

    std::string message_;

public:

    explicit Message () : message_ () {  }
    explicit Message (const char* message) : message_ (message) { }
    explicit Message (const std::string& message) : message_ (message) { }
    ~Message () throw ();

    const char* what () const throw () { return message_.c_str (); }

    void add_data (const char* data) { message_ += data;  }

};

Message& operator<<(Message &message, const char* data);

}
#endif /* _ERRORS_H_ */
