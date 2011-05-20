/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_EXCEPTION_HH_
#define AMANZI_CHEMISTRY_EXCEPTION_HH_

#include <string>

#include "errors.hh"

class ChemistryException : public Errors::Message
{
 public:
  ChemistryException();
  ChemistryException(std::string error_message);
  virtual ~ChemistryException() throw();

  static const std::string kChemistryError;

 protected:
  static const std::string kDefaultMessage;
 private:
};


class ChemistryInvalidInput : public ChemistryException
{
 public:
  ChemistryInvalidInput(std::string error_message);
  virtual ~ChemistryInvalidInput() throw();

 protected:
  ChemistryInvalidInput();

 private:
};


class ChemistryUnrecoverableError : public ChemistryException
{
 public:
  ChemistryUnrecoverableError(std::string error_message);
  virtual ~ChemistryUnrecoverableError() throw();

 protected:
  ChemistryUnrecoverableError();

 private:
};


class ChemistryMaxIterationsReached : public ChemistryException
{
 public:
  ChemistryMaxIterationsReached(std::string error_message);
  virtual ~ChemistryMaxIterationsReached() throw();

 protected:
  ChemistryMaxIterationsReached();

 private:
};


class ChemistryInvalidSolution : public ChemistryException
{
 public:
  ChemistryInvalidSolution(std::string error_message);
  virtual ~ChemistryInvalidSolution() throw();

 protected:
  ChemistryInvalidSolution();

 private:
};

#endif  // AMANZI_CHEMISTRY_EXCEPTION_HH_
