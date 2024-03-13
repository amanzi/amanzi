/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <stdio.h>

#include "errors.hh"
#include "exceptions.hh"
#include "dbc.hh"

#include <xercesc/util/XMLString.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/framework/LocalFileInputSource.hpp>
#include <xercesc/sax/ErrorHandler.hpp>
#include <xercesc/sax/SAXParseException.hpp>
#include <xercesc/validators/common/Grammar.hpp>

class AmanziErrorHandler : public xercesc::ErrorHandler
{
 private:
  void reportParseException(const xercesc::SAXParseException& ex) {
    char* msg = xercesc::XMLString::transcode(ex.getMessage());
    fprintf(stderr, " at line %lu column %lu, %s\n",
            (unsigned long)ex.getLineNumber(), (unsigned long)ex.getColumnNumber(), msg);
    xercesc::XMLString::release(&msg);
  }

 public:
  void warning(const xercesc::SAXParseException& ex) {
    //reportParseException(ex);
    char* msg = xercesc::XMLString::transcode(ex.getMessage());
    fprintf(stderr, "WARNING at line %lu column %lu, %s\n",
            (unsigned long)ex.getLineNumber(), (unsigned long)ex.getColumnNumber(), msg);
    xercesc::XMLString::release(&msg);
  }

  void error(const xercesc::SAXParseException& ex) {
    //reportParseException(ex);
    char* msg = xercesc::XMLString::transcode(ex.getMessage());
    fprintf(stderr, "ERROR at line %lu column %lu, %s\n",
            (unsigned long)ex.getLineNumber(), (unsigned long)ex.getColumnNumber(), msg);
    xercesc::XMLString::release(&msg);
    Exceptions::amanzi_throw(Errors::Message("Errors occured while parsing the input file. Aborting."));
  }

  void fatalError(const xercesc::SAXParseException& ex) {
    //reportParseException(ex);
    char* msg = xercesc::XMLString::transcode(ex.getMessage());
    fprintf(stderr, "FATAL ERROR at line %lu column %lu, %s\n",
            (unsigned long)ex.getLineNumber(), (unsigned long)ex.getColumnNumber(), msg);
    xercesc::XMLString::release(&msg);
    Exceptions::amanzi_throw(Errors::Message("Errors occured while parsing the input file. Aborting."));
  }

  void resetErrors() {};
};
