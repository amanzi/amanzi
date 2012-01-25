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
