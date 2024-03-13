/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
  Output

  Interface for Output wrappers.
  Defines an interface for vis and checkpoint read and writing.
*/

#ifndef AMANZI_OUTPUT_HH_
#define AMANZI_OUTPUT_HH_

#include <string>
#include <vector>

#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"

#include "Mesh.hh"

namespace Amanzi {

class Output {
 public:
  virtual ~Output() {}

  // open and close files
  virtual void InitializeCycle(double time, int cycle, const std::string& tag) = 0;
  virtual void FinalizeCycle() = 0;

  // write data to file
  virtual void WriteVector(const Epetra_Vector& vec,
                           const std::string& name,
                           const AmanziMesh::Entity_kind& kind) const = 0;
  virtual void WriteMultiVector(const Epetra_MultiVector& vec,
                                const std::vector<std::string>& names,
                                const AmanziMesh::Entity_kind& kind) const = 0;

  // can we template this?
  virtual void WriteAttribute(const double& val, const std::string& name) const = 0;
  virtual void WriteAttribute(const int& val, const std::string& name) const = 0;
  virtual void WriteAttribute(const std::string& val, const std::string& name) const = 0;

  // read data from file
  virtual void ReadVector(Epetra_Vector& vec, const std::string& name) const = 0;
  virtual void
  ReadMultiVector(Epetra_MultiVector& vec, const std::vector<std::string>& names) const = 0;

  virtual void ReadAttribute(double& val, const std::string& name) const = 0;
  virtual void ReadAttribute(int& val, const std::string& name) const = 0;
  virtual void ReadAttribute(std::string& val, const std::string& name) const = 0;
};

} // namespace Amanzi

#endif
