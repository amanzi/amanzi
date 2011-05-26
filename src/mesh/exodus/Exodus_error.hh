#ifndef _EXODUS_ERROR_HH_
#define _EXODUS_ERROR_HH_

#include "errors.hh"
#include <sstream>


namespace Amanzi {
namespace Exodus {

class ExodusError : public Errors::Message {
 public:
  explicit ExodusError(void) : Errors::Message() {};
  explicit ExodusError(const char* message) : Errors::Message(message) {};
  ~ExodusError(void) throw() {};
};
  
} // close namespace Exodus
}

#endif
