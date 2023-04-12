/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! Interface for Output implementations.
/*

Defines an interface for writing vis or checkpoint files.

Developer notes:

The design of this library within Amanzi is based on two orthogonal concepts:

1. how to read and write data within a file (File objects)
2. where to read and write data within a file (Input/Output objects).

The former is lower level, the latter includes file structure, metadata, etc.
Input/Output classes use File objects to do the actual writing.

Output types can be tied to their file type, working with one and only one File
object (e.g. OutputSilo or OutputHDF5), or use more than one (e.g. OutputXDMF
uses a FileHDF5 and FileXDMF), specifying a layout and metadata within those
files.

*/

#ifndef AMANZI_OUTPUT_HH_
#define AMANZI_OUTPUT_HH_

#include <string>
#include <vector>
#include "AmanziTypes.hh"
#include "MeshDefs.hh"
#include "CompositeVector_decl.hh"
#include "OutputUtils.hh"

namespace Amanzi {


// Interface for generic visualization.
class Output {
 public:

  // lambda/function for generating a filename
  using FilenameFormatter = std::function<std::string(int)>;

  // Trait for whether a type can be written.
  //
  // This allows easier extension.
  template <typename T>
  struct writes {
    static const bool value = false;
  };

  virtual ~Output() {}

  // open and close files
  virtual void createTimestep(double time, int cycle) = 0;
  virtual void finalizeTimestep() = 0;

  // must return the last written filename, called AFTER finalizeTimestep()
  virtual std::string getFilename(int cycle) const = 0;

  // how nice it would be to template virtual methods...
  virtual void write(const Teuchos::ParameterList& attrs, const bool& val) const {
    write(attrs, (int)val);
  }
  virtual void write(const Teuchos::ParameterList& attrs, const int& val) const = 0;
  virtual void write(const Teuchos::ParameterList& attrs, const double& val) const = 0;
  virtual void write(const Teuchos::ParameterList& attrs, const std::string& val) const = 0;
  virtual void write(const Teuchos::ParameterList& attrs, const Teuchos::Array<int>& val) const {
    writeArray_(attrs, val);
  }
  virtual void write(const Teuchos::ParameterList& attrs, const Teuchos::Array<double>& val) const {
    writeArray_(attrs, val);
  }
  virtual void write(const Teuchos::ParameterList& attrs, const Teuchos::Array<bool>& val) const {
    writeArray_(attrs, val);
  }
  virtual void write(const Teuchos::ParameterList& attrs, const AmanziGeometry::Point& val) const {
    Teuchos::Array<double> arr(val.dim());
    for (int i=0; i!=val.dim(); ++i) arr[i] = val[i];
    writeArray_(attrs, arr);
  }

  virtual void write(const Teuchos::ParameterList& attrs, const Map_type& map) const {
    write(attrs, OutputUtils::asVector(Teuchos::rcpFromRef(map)));
  }
  virtual void write(const Teuchos::ParameterList& attrs, const Vector_type& vec) const = 0;
  virtual void write(const Teuchos::ParameterList& attrs, const IntVector_type& vec) const = 0;

  virtual void write(const Teuchos::ParameterList& attrs, const MultiVector_type& vec) const {
    writeMultiVector_(attrs, vec);
  }
  virtual void write(const Teuchos::ParameterList& attrs, const IntMultiVector_type& vec) const {
    writeMultiVector_(attrs, vec);
  }

  virtual void
  write(const Teuchos::ParameterList& attrs, const CompositeVector_<int_type>& vec) const {
    writeCompositeVector_(attrs, vec);
  }
  virtual void
  write(const Teuchos::ParameterList& attrs, const CompositeVector_<double_type>& vec) const {
    writeCompositeVector_(attrs, vec);
  }

 protected:
  template<typename Scalar>
  void writeMultiVector_(const Teuchos::ParameterList& attrs, const MultiVector_type_<Scalar>& vec) const {
    std::vector<std::string> names = OutputUtils::names(attrs, vec.getNumVectors());
    for (int i = 0; i != vec.getNumVectors(); ++i) {
      Teuchos::ParameterList attrs_i(attrs);
      attrs_i.setName(names[i]);
      write(attrs_i, *vec.getVector(i));
    }
  }

  template<typename Scalar>
  void writeCompositeVector_(const Teuchos::ParameterList& attrs, const CompositeVector_<Scalar>& vec) const {
    for (const auto& compname : vec) {
      std::string vname = attrs.name() + "." + compname;
      Teuchos::ParameterList attrs_comp(attrs);
      attrs_comp.setName(vname);
      write(attrs_comp, *vec.getComponent(compname, false));
    }
  }

  template<typename SerialArray>
  void writeArray_(const Teuchos::ParameterList& attrs, const SerialArray& arr) const {
    for (int i=0; i!=arr.size(); ++i) {
      Teuchos::ParameterList attrs_i(attrs);
      attrs_i.setName(attrs.name() + "_" + std::to_string(i));
      write(attrs_i, arr[i]);
    }
  }

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

template <>
struct Output::writes<AmanziMesh::MeshCache<MemSpace_kind::HOST>> {
  static const bool value = true;
};


} // namespace Amanzi

#endif
