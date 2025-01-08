/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*

Developer Note:

This virtual class and its deriving implementations are meant to be used to
read serial, dense data from file.  They are NOT intended to read distributed
data.

*/

#pragma once

#include <string>
#include <memory>
#include <vector>

#include "Teuchos_Array.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

namespace Amanzi {

class Reader {
 public:
  virtual ~Reader() {}

  virtual bool hasVariableOrGroup(const std::string& name) const = 0;

  void read(const std::string& varname, std::vector<double>& vec, int index = -1) const
  {
    Teuchos::Array<double> vec2;
    read(varname, vec2, index);
    vec = vec2.toVector();
  }
  void read(const std::string& varname, std::vector<int>& vec, int index = -1) const
  {
    Teuchos::Array<int> vec2;
    read(varname, vec2, index);
    vec = vec2.toVector();
  }

  virtual void
  read(const std::string& varname, Teuchos::Array<double>& vec, int index = -1) const = 0;
  virtual void read(const std::string& varname, Teuchos::Array<int>& vec, int index = -1) const = 0;
  virtual void read(const std::string& varname,
                    Teuchos::SerialDenseMatrix<int, double>& vec,
                    int index = -1) const = 0;
  virtual void read(const std::string& varname,
                    Teuchos::SerialDenseMatrix<int, int>& vec,
                    int index = -1) const = 0;
};


// factory, uses filename extension to figure out the type
std::unique_ptr<Reader>
createReader(const std::string& filename);

} // namespace Amanzi
