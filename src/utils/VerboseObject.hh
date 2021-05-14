/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
//! VerboseObject: a controller for writing log files on multiple cores with varying verbosity.

/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Basic VerboseObject for use by Amanzi code.  Trilinos's VerboseObject is
  templated with the class (for no reason) and then requests that the
  VerboseObject be inserted as a base class to the using class.  This is serious
  code smell (composition over inheritance, especially for code reuse).

  I would prefer to get rid of the call to getOSTab(), but I can't figure out
  how to do it.

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

  Parameters:

  <ParameterList name="my class">
    <ParameterList name="verbose object">
      <Parameter name="verbosity level" type="string" value="medium"/>
      <Parameter name="name" type="string" value="my header"/>
      <Parameter name="hide line prefix" type="bool" value="false"/>
      <Parameter name="write on rank" type="int" value="0"/>
    </ParameterList>
  </ParameterList>
*/

/*!

This allows control of log-file verbosity for a wide variety of objects
and physics.

* `"verbosity level`" ``[string]`` **GLOBAL_VERBOSITY**, `"low`", `"medium`", `"high`", `"extreme`"

   The default is set by the global verbosity spec, (fix me!)  Typically,
   `"low`" prints out minimal information, `"medium`" prints out errors and
   overall high level information, `"high`" prints out basic debugging, and
   `"extreme`" prints out local debugging information.

Note: while there are other options, users should typically not need them.
Instead, developers can use them to control output.
   
Example:

.. code-block:: xml

  <ParameterList name="verbose object">
    <Parameter name="verbosity level" type="string" value="medium"/>
    <Parameter name="name" type="string" value="my header"/>
    <Parameter name="hide line prefix" type="bool" value="false"/>
    <Parameter name="write on rank" type="int" value="0"/>
  </ParameterList>

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
  VerboseObject(const Comm_ptr_type& comm, const std::string& name,
                Teuchos::ParameterList plist);

  // this constructor is fragile and should be frowned upon, but is useful for
  // maintaining backward compatibility with the Epetra stack
  VerboseObject(const Comm_type& comm, const std::string& name,
                Teuchos::ParameterList plist) :
      VerboseObject(Teuchos::rcpFromRef(comm), name, plist) {}

  // NOTE: Default destructor, copy construct should be ok.

  // Is process 0 and verbosity is included?
  inline bool os_OK(Teuchos::EVerbosityLevel verbosity) const;

  // Get the stream (errors if !os_OK()).
  inline Teuchos::RCP<Teuchos::FancyOStream> os() const;

  // Simple one-line wrapper
  inline
  void Write(Teuchos::EVerbosityLevel verbosity, const std::string& data) const;
  void WriteWarning(Teuchos::EVerbosityLevel verbosity, const std::stringstream& data) const;
  void WriteWarning(Teuchos::EVerbosityLevel verbosity, const std::string& data) const;

 public:
  // The default global verbosity level.
  static Teuchos::EVerbosityLevel global_default_level;

  // Show or hide line prefixes
  static bool global_hide_line_prefix;

  // Size of the left column of names.
  static unsigned int global_line_prefix_size;

  // Color output for developers
  std::string color(const std::string& name) const;
  std::string reset() const;
  std::string clock() const;

  void set_name(const std::string& name);

 protected:
  Teuchos::RCP<Teuchos::FancyOStream> out_;
  Comm_ptr_type comm_;
};


bool VerboseObject::os_OK(Teuchos::EVerbosityLevel verbosity) const {
  return out_.get() &&
      includesVerbLevel(getVerbLevel(), verbosity, true);
};


Teuchos::RCP<Teuchos::FancyOStream> VerboseObject::os() const {
  return out_;
};


void VerboseObject::Write(Teuchos::EVerbosityLevel verbosity, const std::string& data) const {
  if (getVerbLevel() >= verbosity) {
    Teuchos::OSTab tab = getOSTab();
    *os() << data;
  }
}

}  // namespace Amanzi

#endif
