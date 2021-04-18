/*
  Functions

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Factory for vector functions which are composed of multiple scalar functions.
  The expected plist is of the form:

  <ParameterList name="constuctor plist">
    <Parameter name="number of dofs" type="int" value="N">
    <ParameterList name="dof 1 function">
      <ParameterList name="function-constant">
        ...
      </ParameterList>
    </ParameterList>

    <ParameterList name="dof 2 function">
      <ParameterList name="function-linear">
        ...
      </ParameterList>
    </ParameterList>

    ...
  </ParameterList>

  Where each of the "function-xxx" lists are valid input to the
  function-factory Create() method (see ./function-factory.hh).

  To simplify input for the case where all scalar functions are 
  constant, parameter "values" is supported

  <ParameterList name="constuctor plist">
    <Parameter name="constant values" type="Array(double)" value="{1.0, 2.1, 3.2}">
  </ParameterList>
*/

#include "dbc.hh"
#include "AmanziTypes.hh"
#include "MultiFunction.hh"

namespace Amanzi {

MultiFunction::MultiFunction(
        const std::vector<Teuchos::RCP<const Function> >& functions) :
    functions_(functions) {
  values_ = new double[functions_.size()];
};

MultiFunction::MultiFunction(const Teuchos::RCP<const Function>& function) :
    functions_(1, function)  {
  values_ = new double[1];
};


MultiFunction::MultiFunction(Teuchos::ParameterList& plist) {
  FunctionFactory factory;

  if (plist.isParameter("number of dofs")) {
    if (plist.isType<int>("number of dofs")) {
      int ndofs = plist.get<int>("number of dofs");

      if (ndofs < 1) {
        // ERROR -- invalid number of dofs
        AMANZI_ASSERT(0);
      }

      for (int lcv = 1; lcv != (ndofs+1); ++lcv) {
        std::stringstream sublist_name;
        sublist_name << "dof " << lcv << " function";
        functions_.push_back(Teuchos::rcp(factory.Create(plist.sublist(sublist_name.str()))));
      }
    } else {
      // ERROR -- invalid number of dofs
      AMANZI_ASSERT(0);
    }

  } else if (plist.isParameter("constant values")) {
    auto values = plist.get<Teuchos::Array<double> >("constant values").toVector();
    int ndofs = values.size();

    for (int i = 0; i < ndofs; ++i) {
      Teuchos::ParameterList flist;
      flist.sublist("function-constant").set<double>("value", values[i]);
      functions_.push_back(Teuchos::rcp(factory.Create(flist)));
    }

  } else {
    // assume it is a single dof function
    functions_.push_back(Teuchos::rcp(factory.Create(plist)));
  };

  values_ = new double[functions_.size()];
}


MultiFunction::~MultiFunction() {
  delete [] values_;
};


int MultiFunction::size() const {
  return functions_.size();
};


double* MultiFunction::operator()(const std::vector<double>& xt) const {
  for (int i=0; i!=size(); ++i) {
    values_[i] = (*functions_[i])(xt);
  }
  return values_;
};



} // namespace

