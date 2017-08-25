/* -*-  mode: c++; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_EXCEPTION_HH_
#define AMANZI_CHEMISTRY_EXCEPTION_HH_

#include <string>

#include "errors.hh"

namespace Amanzi {
namespace AmanziChemistry {

class ChemistryException : public Errors::Message {
 public:
  ChemistryException();
  explicit ChemistryException(std::string error_message);
  virtual ~ChemistryException() throw();

 public:
  static const std::string kChemistryError;

 protected:
  static const std::string kDefaultMessage;
};


class ChemistryInvalidInput : public ChemistryException {
 public:
  explicit ChemistryInvalidInput(std::string error_message);
  virtual ~ChemistryInvalidInput() throw();

 protected:
  ChemistryInvalidInput();
};


class ChemistryUnrecoverableError : public ChemistryException {
 public:
  explicit ChemistryUnrecoverableError(std::string error_message);
  virtual ~ChemistryUnrecoverableError() throw();

 protected:
  ChemistryUnrecoverableError();
};


class ChemistryMaxIterationsReached : public ChemistryException {
 public:
  explicit ChemistryMaxIterationsReached(std::string error_message);
  virtual ~ChemistryMaxIterationsReached() throw();

 protected:
  ChemistryMaxIterationsReached();
};


class ChemistryInvalidSolution : public ChemistryException {
 public:
  explicit ChemistryInvalidSolution(std::string error_message);
  virtual ~ChemistryInvalidSolution() throw();

 protected:
  ChemistryInvalidSolution();
};


class ChemistryMemorySizeError : public ChemistryException {
 public:
  explicit ChemistryMemorySizeError(std::string error_message);
  virtual ~ChemistryMemorySizeError() throw();

 protected:
  ChemistryMemorySizeError();
};

}  // namespace AmanziChemistry
}  // namespace Amanzi
#endif  // AMANZI_CHEMISTRY_EXCEPTION_HH_
