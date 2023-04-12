/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! Interface for Input implementations.
/*

Defines an interface for reading checkpoint files.

Developer notes:

The design of this library within Amanzi is based on two orthogonal concepts:

1. how to read data within a file (File objects)
2. where to read data within a file (Input objects).

The former is lower level, the latter includes file structure, metadata, etc.
Input classes use File objects to do the actual writing.  Often a given
Input type (e.g. OutputXDMF) uses more than one File objects (e.g. FileHDF5
and FileXDMF), specifying a layout.


*/

#ifndef AMANZI_INPUT_HH_
#define AMANZI_INPUT_HH_

#include <string>
#include <vector>
#include "AmanziTypes.hh"
#include "MeshDefs.hh"
#include "CompositeVector_decl.hh"
#include "OutputUtils.hh"

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
  virtual void read(const Teuchos::ParameterList& attrs, double& val) const = 0;
  virtual void read(const Teuchos::ParameterList& attrs, int& val) const = 0;
  virtual void read(const Teuchos::ParameterList& attrs, std::string& val) const = 0;
  virtual void read(const Teuchos::ParameterList& attrs, Vector_type& vec) const = 0;
  virtual void read(const Teuchos::ParameterList& attrs, IntVector_type& vec) const = 0;

  virtual void read(const Teuchos::ParameterList& attrs, AmanziGeometry::Point& val) const {
    Teuchos::Array<double> arr(val.dim());
    readArray_(attrs, arr);
    for (int i=0; i!=val.dim(); ++i) val[i] = arr[i];
  }

  virtual void read(const Teuchos::ParameterList& attrs, MultiVector_type& vec) const {
    readMultiVector_(attrs, vec);
  }
  virtual void read(const Teuchos::ParameterList& attrs, IntMultiVector_type& vec) const {
    readMultiVector_(attrs, vec);
  }

  virtual void
  read(const Teuchos::ParameterList& attrs, CompositeVector_<int_type>& vec) const {
    readCompositeVector_(attrs, vec);
  }
  virtual void
  read(const Teuchos::ParameterList& attrs, CompositeVector_<double_type>& vec) const {
    readCompositeVector_(attrs, vec);
  }

 protected:

  template<typename SerialArray>
  void readArray_(const Teuchos::ParameterList& attrs, SerialArray& arr) const {
    for (int i=0; i!=arr.size(); ++i) {
      Teuchos::ParameterList attrs_i(attrs);
      attrs_i.setName(attrs.name() + "_" + std::to_string(i));
      read(attrs_i, arr[i]);
    }
  }

  template<typename Scalar>
  void readMultiVector_(const Teuchos::ParameterList& attrs, MultiVector_type_<Scalar>& vec) const {
    std::vector<std::string> names = OutputUtils::names(attrs, vec.getNumVectors());
    for (int i = 0; i != vec.getNumVectors(); ++i) {
      Teuchos::ParameterList attrs_i(attrs);
      attrs_i.setName(names[i]);
      read(attrs_i, *vec.getVectorNonConst(i));
    }
  }

  template<typename Scalar>
  void readCompositeVector_(const Teuchos::ParameterList& attrs, CompositeVector_<Scalar>& vec) const {
    for (const auto& compname : vec) {
      std::string vname = attrs.name() + "." + compname;
      Teuchos::ParameterList attrs_comp(attrs);
      attrs_comp.setName(vname);
      read(attrs_comp, *vec.getComponent(compname, false));
    }
  }
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
