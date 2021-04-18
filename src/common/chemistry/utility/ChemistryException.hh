/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Jeffrey Johnson <jnjohnson@lbl.gov>

  Simple exception class for chemistry errors.
  Notes:
    - see Eckel, Bruce and Chuck Allison. 2004. Thinking In C++
      Volume 2: Practical Programming. Pearson Prentice Hall. Ch. 1
    - Stroustrup, Bjarne. 2000. The C++ programing language. Special Edition
      Addison Wesley. Ch.8.3 and 14

  Example Use:
    std::ostringstream error_stream;
    error_stream << "ERROR: MyClass::MyFunction(): \n";
    error_stream << "ERROR: something has gone wrong....\n";
    throw ChemistryException(error_stream.str());

    try {
      ...
    }
    catch (ChemistryException& geochem_error) {
      geochem_error.PrintErrorMessage();
      std::cerr << geochem_error.what() << std::endl;
    }
*/

#ifndef AMANZI_CHEMISTRY_EXCEPTION_HH_
#define AMANZI_CHEMISTRY_EXCEPTION_HH_

#include <string>
#include <iostream>

#include "errors.hh"

namespace Amanzi {
namespace AmanziChemistry {

const std::string kChemistryError = "CHEMISTRY_ERROR: ";
const std::string kDefaultMessage = "An unknown error has occured.";

class ChemistryException : public Errors::Message {
 public:
  ChemistryException() : Message(kChemistryError) { this->add_data(kDefaultMessage.c_str()); }
  explicit ChemistryException(std::string msg) : Message(kChemistryError) { this->add_data(msg.c_str()); }
  virtual ~ChemistryException() noexcept {};
};


class ChemistryInvalidInput : public ChemistryException {
 public:
  explicit ChemistryInvalidInput(std::string msg) : ChemistryException(msg) {};
  virtual ~ChemistryInvalidInput() noexcept {};

 protected:
  ChemistryInvalidInput() : ChemistryException() {};
};


class ChemistryUnrecoverableError : public ChemistryException {
 public:
  explicit ChemistryUnrecoverableError(std::string msg) : ChemistryException(msg) {};
  virtual ~ChemistryUnrecoverableError() noexcept {};

 protected:
  ChemistryUnrecoverableError() : ChemistryException() {};
};


class ChemistryMaxIterationsReached : public ChemistryException {
 public:
  explicit ChemistryMaxIterationsReached(std::string msg) : ChemistryException(msg) {};
  virtual ~ChemistryMaxIterationsReached() noexcept {};

 protected:
  ChemistryMaxIterationsReached() : ChemistryException() {};
};


class ChemistryInvalidSolution : public ChemistryException {
 public:
  explicit ChemistryInvalidSolution(std::string msg) : ChemistryException(msg) {};
  virtual ~ChemistryInvalidSolution() noexcept {};

 protected:
  ChemistryInvalidSolution() : ChemistryException() {};
};


class ChemistryMemorySizeError : public ChemistryException {
 public:
  explicit ChemistryMemorySizeError(std::string msg) : ChemistryException(msg) {};
  virtual ~ChemistryMemorySizeError() noexcept {};

 protected:
  ChemistryMemorySizeError() : ChemistryException() {};
};

}  // namespace AmanziChemistry
}  // namespace Amanzi
#endif
