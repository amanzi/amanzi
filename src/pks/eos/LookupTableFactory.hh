/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  EOS

  Factory of lookup tables.
*/

#ifndef AMANZI_LOOKUP_TABLE_FACTORY_HH_
#define AMANZI_LOOKUP_TABLE_FACTORY_HH_

#include <fstream>

#include "Teuchos_ParameterList.hpp"

#include "errors.hh"

#include "LookupTable_Amanzi.hh"
#include "LookupTable_FEHM.hh"

namespace Amanzi {
namespace AmanziEOS {

inline Teuchos::RCP<LookupTable>
CreateLookupTable(Teuchos::ParameterList& plist)
{
  // identify file format
  std::ifstream ifs;
  std::string filename = plist.get<std::string>("table name");

  Errors::Message msg;
  msg << "\nFailed to open/read data from file: " << filename;

  ifs.open(filename, std::ifstream::in);
  if (ifs.fail()) Exceptions::amanzi_throw(msg);

  char line[100];
  ifs.getline(line, 100);
  if (ifs.fail()) Exceptions::amanzi_throw(msg);
  std::string format = (strncmp(line, "Amanzi", 6) == 0) ? "Amanzi" : "FEHM";
  ifs.close();

  // select a file reader
  if (format == "Amanzi") {
    return Teuchos::rcp(new LookupTable_Amanzi(plist));
  } else if (format == "FEHM") {
    return Teuchos::rcp(new LookupTable_FEHM(plist));
  } else {
    Errors::Message msg;
    msg << "\nUnknown format for the lookup table: \"" << format << "\"\n";
    Exceptions::amanzi_throw(msg);
    return Teuchos::null;
  }
}


} // namespace AmanziEOS
} // namespace Amanzi

#endif
