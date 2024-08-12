/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*

Developer Note:

This virtual class and its deriving implementations are meant to be used to
read serial, dense data from file.  They are NOT intended to read distributed
data.

*/

#pragma once

#include <string>

namespace Teuchos {
template<typename LO, typename Scalar> class SerialDenseVector;
template<typename LO, typename Scalar> class SerialDenseMatrix;
}

namespace Amanzi {

class Reader {
 public:
  virtual ~Reader() {}

  virtual bool hasVariableOrGroup(const std::string& name) const = 0;
  virtual void read(const std::string& varname, Teuchos::SerialDenseVector<int, double>& vec, int index = -1) const = 0;
  virtual void read(const std::string& varname, Teuchos::SerialDenseVector<int, int>& vec, int index = -1) const = 0;
  virtual void read(const std::string& varname, Teuchos::SerialDenseMatrix<int, double>& vec, int index = -1) const = 0;
  virtual void read(const std::string& varname, Teuchos::SerialDenseMatrix<int, int>& vec, int index = -1) const = 0;
};

} // namespace Amanzi
