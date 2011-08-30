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
#include "chemistry_exception.hh"

#include <string>
#include <iostream>

namespace amanzi {
namespace chemistry {

/*******************************************************************************
 **
 ** ChemistryException
 **
 *******************************************************************************/

const std::string ChemistryException::kChemistryError = "CHEMISTRY_ERROR: ";
const std::string ChemistryException::kDefaultMessage = "An unknown error has occured.";

ChemistryException::ChemistryException()
    : Message(ChemistryException::kChemistryError) {
  this->add_data(kDefaultMessage.c_str());
  /* end ChemistryException() */
}

ChemistryException::ChemistryException(std::string error_message)
    : Message(ChemistryException::kChemistryError) {
  this->add_data(error_message.c_str());
  /* end ChemistryException(std::string) */
}

ChemistryException::~ChemistryException() throw() {
  /* end ~ChemistryException() */
}

/*******************************************************************************
 **
 ** ChemistryInvalidInput
 **
 *******************************************************************************/

ChemistryInvalidInput::ChemistryInvalidInput()
    : ChemistryException() {
  /* end ChemistryInvalidInput() */
}

ChemistryInvalidInput::ChemistryInvalidInput(std::string error_message)
    : ChemistryException(error_message) {
  /* end ChemistryInvalidInput(std::string) */
}

ChemistryInvalidInput::~ChemistryInvalidInput() throw() {
  /* end ~ChemistryInvalidInput() */
}

/*******************************************************************************
 **
 ** ChemistryUnrecoverableError
 **
 *******************************************************************************/

ChemistryUnrecoverableError::ChemistryUnrecoverableError()
    : ChemistryException() {
  /* end ChemistryUnrecoverableError() */
}

ChemistryUnrecoverableError::ChemistryUnrecoverableError(std::string error_message)
    : ChemistryException(error_message) {
  /* end ChemistryUnrecoverableError(std::string) */
}

ChemistryUnrecoverableError::~ChemistryUnrecoverableError() throw() {
  /* end ~ChemistryUnrecoverableError() */
}

/*******************************************************************************
 **
 ** ChemistryMaxIterationsReached
 **
 *******************************************************************************/

ChemistryMaxIterationsReached::ChemistryMaxIterationsReached()
    : ChemistryException() {
  /* end ChemistryMaxIterationsReached() */
}

ChemistryMaxIterationsReached::ChemistryMaxIterationsReached(std::string error_message)
    : ChemistryException(error_message) {
  /* end ChemistryMaxIterationsReached(std::string) */
}

ChemistryMaxIterationsReached::~ChemistryMaxIterationsReached() throw() {
  /* end ~ChemistryMaxIterationsReached() */
}

/*******************************************************************************
 **
 ** ChemistryInvalidSolution
 **
 *******************************************************************************/

ChemistryInvalidSolution::ChemistryInvalidSolution()
    : ChemistryException() {
  /* end ChemistryInvalidSolution() */
}

ChemistryInvalidSolution::ChemistryInvalidSolution(std::string error_message)
    : ChemistryException(error_message) {
  /* end ChemistryInvalidSolution(std::string) */
}

ChemistryInvalidSolution::~ChemistryInvalidSolution() throw() {
  /* end ~ChemistryInvalidSolution() */
}

}  // namespace chemistry
}  // namespace amanzi
