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

#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Epetra_MpiComm.h"

#include "VerboseObject_objs.hh"

namespace Amanzi {

/* ******************************************************************
* Create parameter list with the given verbosity.
****************************************************************** */
VerboseObject::VerboseObject(std::string name, const std::string& verbosity) :
    comm_(NULL)
{
  setDefaultVerbLevel(global_default_level);
  set_name(name);

  // Override from ParameterList.
  Teuchos::ParameterList plist;
  plist.sublist("VerboseObject").set("Verbosity Level", verbosity);

  Teuchos::readVerboseObjectSublist(&plist, this);

  // out, tab
  out_ = getOStream();
  out_->setShowLinePrefix(hide_line_prefix);
}


/* ******************************************************************
* Read verbosity from a parameter list
****************************************************************** */
VerboseObject::VerboseObject(std::string name, Teuchos::ParameterList& plist) :
    comm_(NULL)
{
  // Set up the default level.
  setDefaultVerbLevel(global_default_level);

  // Options from ParameterList

  // -- Set up the VerboseObject header.
  std::string headername(name);
  if (plist.sublist("verbose object").isParameter("name")) {
    headername = plist.sublist("verbose object").get<std::string>("name");
  }
  set_name(headername);

  // -- Show the line prefix
  bool no_pre = plist.sublist("verbose object").get<bool>("hide line prefix", hide_line_prefix);

  // Override from ParameterList.
  Teuchos::ParameterList plist_out;
  if (plist.sublist("verbose object").isParameter("verbosity level")) {
    plist_out.sublist("VerboseObject").set("Verbosity Level",
        plist.sublist("verbose object").get<std::string>("verbosity level"));
  }

  Teuchos::readVerboseObjectSublist(&plist_out, this);

  // out, tab
  out_ = getOStream();
  out_->setShowLinePrefix(!no_pre);
}


VerboseObject::VerboseObject(const Epetra_MpiComm* const comm, std::string name,
                             Teuchos::ParameterList& plist) :
    comm_(comm)
{
  int root = -1;
  // Check if we are in the mode of writing only a specific rank.
  if (plist.sublist("verbose object").isParameter("write on rank")) {
    root = plist.sublist("verbose object").get<int>("write on rank");
  }

  // Init the basics
  // Set up the default level.
  setDefaultVerbLevel(global_default_level);

  // Options from ParameterList

  // -- Set up the VerboseObject header.
  std::string headername(name);
  if (plist.sublist("verbose object").isParameter("name")) {
    std::string headername = plist.sublist("verbose object").get<std::string>("name");
  }
  set_name(headername);

  // -- Show the line prefix
  bool no_pre = plist.sublist("verbose object").get<bool>("hide line prefix", hide_line_prefix);

  // Override from ParameterList.
  Teuchos::ParameterList plist_out;
  if (plist.sublist("verbose object").isParameter("verbosity level")) {
    plist_out.sublist("VerboseObject").set("Verbosity Level",
        plist.sublist("verbose object").get<std::string>("verbosity level"));
  }

  Teuchos::readVerboseObjectSublist(&plist_out, this);

  // out, tab
  out_ = getOStream();
  out_->setShowLinePrefix(!no_pre);

  // Set up a local FancyOStream
  if (root >= 0) {
    int size = comm_->NumProc();
    int pid = comm_->MyPID();
    Teuchos::RCP<Teuchos::FancyOStream> newout = Teuchos::rcp(new Teuchos::FancyOStream(out_->getOStream()));
    newout->setProcRankAndSize(pid,size);
    newout->setOutputToRootOnly(root);
    setOStream(newout);

    std::stringstream headerstream;
    headerstream << pid << ": " << getLinePrefix();
    std::string header = headerstream.str();
    if (header.size() > line_prefix_size) {
      header.erase(line_prefix_size);
    } else if (header.size() < line_prefix_size) {
      header.append(line_prefix_size - header.size(), ' ');
    }

    setLinePrefix(header);
    out_ = getOStream();
    out_->setShowLinePrefix(!no_pre);
  }
}


void VerboseObject::set_name(std::string name)
{
  std::string header(name);
  if (header.size() > line_prefix_size) {
    header.erase(line_prefix_size);
  } else if (header.size() < line_prefix_size) {
    header.append(line_prefix_size - header.size(), ' ');
  }
  setLinePrefix(header);
}


std::string VerboseObject::color(std::string name) const
{ 
  std::string output("");
#ifdef __linux
  if (name == "red") {
    output = std::string("\033[1;31m");
  } else if (name == "green") {
    output = std::string("\033[1;32m");
  } else if (name == "yellow") {
    output = std::string("\033[1;33m");
  }
#endif
  return output;
}


std::string VerboseObject::reset() const
{ 
  std::string output("");
#ifdef __linux
  output = std::string("\033[0m");
#endif
  return output;
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
