/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*!

This allows control of log-file verbosity for a wide variety of objects
and physics.

.. _verbose-object-spec:
.. admonition:: verbose-object-spec

   * `"verbosity level`" ``[string]`` **GLOBAL_VERBOSITY**, `"low`",
     `"medium`", `"high`", `"extreme`" The default is set by the global
     verbosity spec, (fix me!)  Typically, `"low`" prints out minimal
     information, `"medium`" prints out errors and overall high level
     information, `"high`" prints out basic debugging, and `"extreme`" prints
     out local debugging information.

   * `"write on rank`" ``[int]`` **0** VerboseObjects only write on a single
     rank -- by deafult the 0th rank.  However, sometimes it is useful for
     debugging to write from another rank due to a need for cleaner output or
     writing a specific cell/entity information.

   * `"output filename`" ``[string]`` **optional** Redirect this output to a
     specific file rather than writing to screen.  Note this will be done
     by-the-instance, so this may not catch as much as one might think.

Example:

.. code-block:: xml

  <ParameterList name="verbose object">
    <Parameter name="verbosity level" type="string" value="medium"/>
    <Parameter name="name" type="string" value="my header"/>
    <Parameter name="hide line prefix" type="bool" value="false"/>
    <Parameter name="write on rank" type="int" value="0"/>
  </ParameterList>

*/

/*

Developer notes:

Basic VerboseObject for use by Amanzi code.  Trilinos's VerboseObject is
templated with the class (for no reason) and then requests that the
VerboseObject be inserted as a base class to the using class.  This is serious
code smell (composition over inheritance, especially for code reuse).

I would prefer to get rid of the call to getOSTab(), but I can't figure out how
to do it.

Usage:

  class MyClass {
   public:
    MyClass(Teuchos::ParameterList& plist) {
      vo = Teuchos::rcp(new VerboseObject("my_class", plist);
      Teuchos::OSTab tab = vo.getOSTab();

      if (vo.os_OK(Teuchos::VERB_MEDIUM)) {
        *vo.os() << "my string to print" << std::endl;
      }
    }

   protected:
    Teuchos::RCP<VerboseObject> vo;
  }

*/


#ifndef AMANZI_VERBOSE_OBJECT_HH_
#define AMANZI_VERBOSE_OBJECT_HH_

#include <string>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "AmanziTypes.hh"

namespace Amanzi {

class VerboseObject : public Teuchos::VerboseObject<VerboseObject> {
 public:
  // Constructors
  VerboseObject(const std::string& name, const std::string& verbosity);
  VerboseObject(const std::string& name, Teuchos::ParameterList plist);
  VerboseObject(const Comm_ptr_type& comm, const std::string& name, Teuchos::ParameterList plist);

  // this constructor is fragile and should be frowned upon, but is useful for
  // maintaining backward compatibility with the Epetra stack
  VerboseObject(const Comm_type& comm, const std::string& name, Teuchos::ParameterList plist)
    : VerboseObject(Teuchos::rcpFromRef(comm), name, plist)
  {}

  // NOTE: Default destructor, copy construct should be ok.

  // Is process 0 and verbosity is included?
  inline bool os_OK(Teuchos::EVerbosityLevel verbosity) const;

  // is the verbosity included?
  inline bool includesVerbLevel(Teuchos::EVerbosityLevel verbosity) const;

  // Get the stream (errors if !os_OK()).
  inline Teuchos::RCP<Teuchos::FancyOStream> os() const;

  // Simple one-line wrapper
  inline void Write(Teuchos::EVerbosityLevel verbosity, const std::string& data) const;
  void WriteWarning(Teuchos::EVerbosityLevel verbosity, const std::stringstream& data) const;
  void WriteWarning(Teuchos::EVerbosityLevel verbosity, const std::string& data) const;

 public:
  // The default global verbosity level.
  static Teuchos::EVerbosityLevel global_default_level;

  // Show or hide line prefixes
  static bool global_hide_line_prefix;

  // Size of the left column of names.
  static unsigned int global_line_prefix_size;

  // static rank to write on
  static unsigned int global_writing_rank;

  // Color output for developers
  std::string color(const std::string& name) const;
  std::string reset() const;
  std::string clock() const;

  // width parameter provides width of the header
  void set_name(std::string name, int width = -1);

 protected:
  Comm_ptr_type comm_;
};


bool
VerboseObject::os_OK(Teuchos::EVerbosityLevel verbosity) const
{
  return getOStream().get() && includesVerbLevel(verbosity);
};

bool
VerboseObject::includesVerbLevel(Teuchos::EVerbosityLevel verbosity) const
{
  return Teuchos::includesVerbLevel(getVerbLevel(), verbosity, true);
}


Teuchos::RCP<Teuchos::FancyOStream>
VerboseObject::os() const
{
  return getOStream();
};


void
VerboseObject::Write(Teuchos::EVerbosityLevel verbosity, const std::string& data) const
{
  if (getVerbLevel() >= verbosity) {
    Teuchos::OSTab tab = getOSTab();
    *os() << data;
  }
}

} // namespace Amanzi

#endif
