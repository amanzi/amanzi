/*
  Utils

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Basic VerboseObject for use by Amanzi code.  Trilinos's VerboseObject is
  templated with the class (for no reason) and then requests that the
  VerboseObject be inserted as a base class to the using class.  This is serious
  code smell (composition over inheritance, especially for code reuse).
*/

#include <ctime>

#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "AmanziComm.hh"
#include "Key.hh"
#include "VerboseObject.hh"

namespace Amanzi {

/* ******************************************************************
* Create parameter list with the given verbosity.
****************************************************************** */
VerboseObject::VerboseObject(const std::string& name, const std::string& verbosity)
{
  setDefaultVerbLevel(global_default_level);
  set_name(name);

  auto validator = Teuchos::verbosityLevelParameterEntryValidator("verbosity level");
  auto ilevel = validator->getIntegralValue(verbosity);
  setVerbLevel(ilevel);

  // out, tab
  getOStream()->setShowLinePrefix(global_hide_line_prefix);
}


/* ******************************************************************
* Read verbosity from a parameter list
****************************************************************** */
VerboseObject::VerboseObject(const std::string& name, Teuchos::ParameterList plist)
{
  // Options from ParameterList
  // Set up the default level.
  setDefaultVerbLevel(global_default_level);

  // -- Set up the VerboseObject header.
  std::string headername(name);
  if (plist.sublist("verbose object").isParameter("name")) {
    headername = plist.sublist("verbose object").get<std::string>("name");
  }
  set_name(headername);

  // -- Show the line prefix
  bool no_pre = plist.sublist("verbose object").get<bool>("hide line prefix", global_hide_line_prefix);

  // Override from ParameterList.
  Teuchos::ParameterList plist_out;
  if (plist.sublist("verbose object").isParameter("verbosity level")) {
    plist_out.sublist("VerboseObject").set("Verbosity Level",
        plist.sublist("verbose object").get<std::string>("verbosity level"));
  }

  if (plist.sublist("verbose object").isParameter("output filename")) {
    plist_out.sublist("VerboseObject").set("Output File",
        plist.sublist("verbose object").get<std::string>("output filename"));
  }

  Teuchos::readVerboseObjectSublist(&plist_out, this);

  // out, tab
  getOStream()->setShowLinePrefix(!no_pre);
}


VerboseObject::VerboseObject(const Comm_ptr_type& comm, const std::string& name,
                             Teuchos::ParameterList outer_plist) :
    comm_(comm)
{
  Teuchos::ParameterList& plist = outer_plist.sublist("verbose object");

  // Make sure we have a unique OStream to work with in this VO
  auto ostream = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  setOverridingOStream(ostream);

  // First options that act on this
  // Options from ParameterList
  // -- Set up the VerboseObject header.
  std::string headername(name);
  if (plist.isParameter("name")) {
    headername = plist.get<std::string>("name");
  }
  int width = plist.get("header width", -1);
  set_name(headername, width);

  // -- Set the verbosity level
  if (plist.isParameter("verbosity level")) {
    auto level = plist.get<std::string>("verbosity level");
    auto validator = Teuchos::verbosityLevelParameterEntryValidator("verbosity level");
    auto ilevel = validator->getIntegralValue(level);
    setVerbLevel(ilevel);
  }

  // Next options that act on OStream
  // -- Show the line prefix
  bool no_pre = plist.get<bool>("hide line prefix", global_hide_line_prefix);
  getOStream()->setShowLinePrefix(!no_pre);

  // -- Include rank in the prefix
  bool show_rank = plist.get<bool>("show rank", false);
  getOStream()->setShowProcRank(show_rank);

  // -- set the comm info
  int size = comm_->NumProc();
  int pid = comm_->MyPID();
  getOStream()->setProcRankAndSize(pid, size);

  // -- write from a different rank than 0
  int root = plist.get<int>("write on rank", (int) global_writing_rank);
  getOStream()->setOutputToRootOnly(root);
}


void VerboseObject::set_name(std::string name, int width)
{
  if (width < 0) width = global_line_prefix_size;

  if (name.size() > width) name = Keys::abbreviate(name, width);

  // hard cut/pad to size
  if (name.size() > width) {
    name.erase(width);
  } else if (name.size() < width) {
    name.append(width - name.size(), ' ');
  }
  setLinePrefix(name);
}


std::string VerboseObject::color(const std::string& name) const
{
  std::string output("");
  if (name == "red") {
    output = std::string("\033[1;31m");
  } else if (name == "green") {
    output = std::string("\033[1;32m");
  } else if (name == "yellow") {
    output = std::string("\033[1;33m");
  } else if (name == "good") {
    output = color("green");
  } else if (name == "bad") {
    output = color("red");
  }
  return output;
}


std::string VerboseObject::reset() const
{
  std::string output("");
  output = std::string("\033[0m");
  return output;
}


std::string VerboseObject::clock() const
{
  int tmp = std::clock() / CLOCKS_PER_SEC;
  int s = tmp % 60;
  tmp /= 60;

  int m = tmp % 60;
  tmp /= 60;

  int h = tmp % 100;

  std::stringstream ss;
  ss << "[" << std::setfill('0')
     << std::setw(2) << h << ":"
     << std::setw(2) << m << ":"
     << std::setw(2) << s << "]";

  return ss.str();
}

void VerboseObject::WriteWarning(Teuchos::EVerbosityLevel verbosity, const std::stringstream& data) const
{
  if (getVerbLevel() >= verbosity) {
    Teuchos::OSTab tab = getOSTab();

    std::string str(data.str());
    size_t pos = str.length() - 1;
    if (str[pos] == '\n') str.erase(pos, 1);
    *os() << color("yellow") << str.erase(pos, 1) << reset() << std::endl;
  }
}


void VerboseObject::WriteWarning(Teuchos::EVerbosityLevel verbosity, const std::string& data) const
{
  if (getVerbLevel() >= verbosity) {
    Teuchos::OSTab tab = getOSTab();
    *os() << color("yellow") << data << reset();
  }
}

}  // namespace Amanzi
