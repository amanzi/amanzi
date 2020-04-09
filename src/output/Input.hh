/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! Interface for Input implementations.

/*
  Defines an interface for vis writing.
*/

#ifndef AMANZI_INPUT_HH_
#define AMANZI_INPUT_HH_

#include <string>
#include <vector>
#include "AmanziTypes.hh"
#include "MeshDefs.hh"
#include "CompositeVector_decl.hh"

namespace Amanzi {

// Interface for generic visualization.
class Input {
 public:
  template <typename T>
  struct reads {
    static const bool value = false;
  };

  virtual ~Input() {}

  // read data from file
  virtual void Read(const Teuchos::ParameterList& attrs, double& val) const = 0;
  virtual void Read(const Teuchos::ParameterList& attrs, int& val) const = 0;
  virtual void
  Read(const Teuchos::ParameterList& attrs, std::string& val) const = 0;
  virtual void
  Read(const Teuchos::ParameterList& attrs, Vector_type& vec) const = 0;
  virtual void
  Read(const Teuchos::ParameterList& attrs, IntVector_type& vec) const = 0;
  virtual void
  Read(const Teuchos::ParameterList& attrs, MultiVector_type& vec) const = 0;
  virtual void
  Read(const Teuchos::ParameterList& attrs, IntMultiVector_type& vec) const = 0;
  virtual void Read(const Teuchos::ParameterList& attrs,
                    CompositeVector_<double>& vec) const = 0;
  virtual void Read(const Teuchos::ParameterList& attrs,
                    CompositeVector_<int>& vec) const = 0;
};

template <>
struct Input::reads<int> {
  static const bool value = true;
};
template <>
struct Input::reads<double> {
  static const bool value = true;
};
template <>
struct Input::reads<std::string> {
  static const bool value = true;
};
template <>
struct Input::reads<Vector_type> {
  static const bool value = true;
};
template <>
struct Input::reads<IntVector_type> {
  static const bool value = true;
};
template <>
struct Input::reads<MultiVector_type> {
  static const bool value = true;
};
template <>
struct Input::reads<IntMultiVector_type> {
  static const bool value = true;
};
template <>
struct Input::reads<CompositeVector_<int>> {
  static const bool value = true;
};
template <>
struct Input::reads<CompositeVector_<double>> {
  static const bool value = true;
};

} // namespace Amanzi

#endif
