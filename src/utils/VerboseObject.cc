/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! <MISSING_ONELINE_DOCSTRING>

#include <ctime>

#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "VerboseObject.hh"

namespace Amanzi {

/* ******************************************************************
 * Create parameter list with the given verbosity.
 ****************************************************************** */
VerboseObject::VerboseObject(const std::string& name,
                             const std::string& verbosity,
                             const Comm_ptr_type& comm)
    : comm_(comm)
{
  Teuchos::ParameterList plist;
  plist.sublist("verbose object").set("verbosity level", verbosity);
  Init_(name, plist);
}


/* ******************************************************************
 * Read verbosity from a parameter list
 ****************************************************************** */
VerboseObject::VerboseObject(const std::string& name,
                             Teuchos::ParameterList& plist,
                             const Comm_ptr_type& comm)
    : comm_(comm)
{
  Init_(name, plist);
}

/* ******************************************************************
 * Read verbosity from a parameter list
 ****************************************************************** */
VerboseObject::VerboseObject(const std::string& name,
                             const Teuchos::ParameterList& plist,
                             const Comm_ptr_type& comm)
    : comm_(comm)
{
  Teuchos::ParameterList new_plist(plist);
  Init_(name, new_plist);
}


void
VerboseObject::Init_(const std::string& name, Teuchos::ParameterList& plist)
{
  // Options from ParameterList
  // Set up the default level.
  setDefaultVerbLevel(global_default_level);

  // Set the fancy ostream
  out_ = getOStream();

  // set up with comm
  if (comm_.get()) {
    int size = comm_->getSize();
    int pid = comm_->getRank();
    if (size > 1) { // this creates a new OStream which can result in shuffling
                    // of lines.  If we don't have to do this, don't!
      out_->setProcRankAndSize(pid, size);

      // Options from ParameterList
      int root = -2;
      // Check if we are in the mode of writing only a specific rank.
      if (plist.sublist("verbose object").isParameter("write on rank")) {
        root = plist.sublist("verbose object").get<int>("write on rank");
      } else if (plist.sublist("verbose object").isParameter("write on all")) {
        root = -1;
      }

      // Set up a local FancyOStream
      if (root >= -1) {
        Teuchos::RCP<Teuchos::FancyOStream> newout =
            Teuchos::rcp(new Teuchos::FancyOStream(out_->getOStream()));
        newout->setProcRankAndSize(pid, size);
        newout->setOutputToRootOnly(root);
        setOStream(newout);
        out_ = getOStream();
      }
    }
  }
  out_->setMaxLenLinePrefix(global_line_prefix_size);

  // -- Set up the VerboseObject header.
  std::string headername(name);
  if (plist.sublist("verbose object").isParameter("name")) {
    headername = plist.sublist("verbose object").get<std::string>("name");
  }
  setLinePrefix(headername);

  // -- Show the line prefix
  bool no_pre = plist.sublist("verbose object")
                  .get<bool>("hide line prefix", global_hide_line_prefix);
  out_->setShowLinePrefix(!no_pre);

  // Override from ParameterList.
  Teuchos::ParameterList plist_out;
  if (plist.sublist("verbose object").isParameter("verbosity level")) {
    plist_out.sublist("VerboseObject")
      .set("Verbosity Level",
           plist.sublist("verbose object").get<std::string>("verbosity level"));
  }
  Teuchos::readVerboseObjectSublist(&plist_out, this);
}    


std::string VerboseObject::color(const std::string& name) const
{ 
  std::string output("");
#ifdef __linux
  if (name == "red") {
    output = std::string("\033[1;31m");
  } else if (name == "green") {
    output = std::string("\033[1;32m");
  } else if (name == "yellow") {
    output = std::string("\033[1;33m");
  } else {
    output = std::string("\033[0m");
  }
#endif
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

} // namespace Amanzi
