/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

License: see $AMANZI_DIR/COPYRIGHT
Author: Ethan Coon

Basic VerboseObject for use by Amanzi code.  Trilinos's VerboseObject is
templated with the class (for no reason) and then requests that the
VerboseObject be inserted as a base class to the using class.  This is serious
code smell (composition over inheritance, especially for code reuse).

------------------------------------------------------------------------- */

#include "Teuchos_VerboseObjectParameterListHelpers.hpp"

#include "VerboseObject.hh"

namespace Amanzi {

VerboseObject::VerboseObject(std::string name, Teuchos::ParameterList& plist) 
{
  // Set up the default level.
  setDefaultVerbLevel(global_default_level);

  // Options from ParameterList

  // -- Set up the VerboseObject header.
  std::string headername = plist.sublist("VerboseObject").get<std::string>("Name",name);
  plist.sublist("VerboseObject").remove("Name");

  std::string header(headername);
  if (header.size() > line_prefix_size) {
    header.erase(line_prefix_size);
  } else if (header.size() < line_prefix_size) {
    header.append(line_prefix_size - header.size(), ' ');
  }
  setLinePrefix(header);

  // -- Show the line prefix
  bool no_pre = plist.sublist("VerboseObject").get<bool>("Hide Line Prefix", hide_line_prefix);
  plist.sublist("VerboseObject").remove("Hide Line Prefix");

  // Override from ParameterList.
  Teuchos::readVerboseObjectSublist(&plist,this);

  // out, tab
  out_ = getOStream();
  out_->setShowLinePrefix(!no_pre);
}

} // namespace Amanzi
