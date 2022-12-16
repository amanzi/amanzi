/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#ifndef TABULAR_FUNCTION_H_
#define TABULAR_FUNCTION_H_

#include <iostream>

#include <REAL.H>
#include <Array.H>

class TabularFunction {
public:

    enum Form {LINEAR, CONSTANT};

    TabularFunction() {}
    TabularFunction(const Array<double>&      x,
                    const Array<double>&      y,
                    const Array<std::string>& form);

    Real operator() (Real x) const;

    Array<double> x_, y_;
    Array<Form> form_;
};

#endif  // TABULAR_FUNCTION_H_
