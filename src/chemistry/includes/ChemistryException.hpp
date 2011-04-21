/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef __CHEMISTRY_EXCEPTION__

#define __CHEMISTRY_EXCEPTION__

#include <string>

#include "errors.hh"

class ChemistryException : public Errors::Message
{
 public:
  enum Status {
    kOkay, kInvalidSolution, kMaxIterationsExceeded, kUnrecoverableError, kRecoverableError
  };

  ChemistryException();
  ChemistryException(std::string error_message);
  ChemistryException(std::string error_message, Status error_status);
  virtual ~ChemistryException() throw();

  void PrintErrorMessage(void) const;
  Status error_status(void) const { return error_status_; };

 protected:
  Status error_status_;
 private:

};
#endif   /* __CHEMISTRY_EXCEPTION__ */
