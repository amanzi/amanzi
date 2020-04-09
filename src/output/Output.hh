/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! Interface for Output implementations.

/*
  Defines an interface for writing vis.
*/

#ifndef AMANZI_OUTPUT_HH_
#define AMANZI_OUTPUT_HH_

#include <string>
#include <vector>
#include "AmanziTypes.hh"
#include "MeshDefs.hh"
#include "CompositeVector_decl.hh"

namespace Amanzi {

// Interface for generic visualization.
class Output {
 public:
  // Trait for whether a type can be written.
  //
  // This allows easier extension.
  template <typename T>
  struct writes {
    static const bool value = false;
  };

  virtual ~Output() {}

  // open and close files
  virtual void CreateFile(double time, int cycle) = 0;
  virtual void FinalizeFile() = 0;
  virtual std::string Filename() const = 0;

  // how nice it would be to template virtual methods...
  virtual void
  Write(const Teuchos::ParameterList& attrs, const int& val) const = 0;
  virtual void
  Write(const Teuchos::ParameterList& attrs, const double& val) const = 0;
  virtual void
  Write(const Teuchos::ParameterList& attrs, const std::string& val) const = 0;

  virtual void
  Write(const Teuchos::ParameterList& attrs, const Vector_type& vec) const = 0;
  virtual void Write(const Teuchos::ParameterList& attrs,
                     const IntVector_type& vec) const = 0;

  virtual void Write(const Teuchos::ParameterList& attrs,
                     const MultiVector_type& vec) const = 0;
  virtual void Write(const Teuchos::ParameterList& attrs,
                     const IntMultiVector_type& vec) const = 0;

  virtual void Write(const Teuchos::ParameterList& attrs,
                     const CompositeVector_<int>& vec) const = 0;
  virtual void Write(const Teuchos::ParameterList& attrs,
                     const CompositeVector_<double>& vec) const = 0;
};


template <>
struct Output::writes<int> {
  static const bool value = true;
};
template <>
struct Output::writes<double> {
  static const bool value = true;
};
template <>
struct Output::writes<std::string> {
  static const bool value = true;
};
template <>
struct Output::writes<Vector_type> {
  static const bool value = true;
};
template <>
struct Output::writes<IntVector_type> {
  static const bool value = true;
};
template <>
struct Output::writes<MultiVector_type> {
  static const bool value = true;
};
template <>
struct Output::writes<IntMultiVector_type> {
  static const bool value = true;
};
template <>
struct Output::writes<CompositeVector_<int>> {
  static const bool value = true;
};
template <>
struct Output::writes<CompositeVector_<double>> {
  static const bool value = true;
};


} // namespace Amanzi

#endif
