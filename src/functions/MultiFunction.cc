/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! Function from R^d to R^n.

/*!

  A MultiFunction is simply an array of functions, which allow Functions to
  be used for MultiVectors.

  Factory for vector functions which are composed of multiple scalar functions.
  The expected plist is of the form:

  <ParameterList name="constuctor plist">
    <Parameter name="number of dofs">
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

  Where each of the "Function X" lists are valid input to the
  function-factory Create() method (see ./function-factory.hh).

*/

#include "dbc.hh"
#include "MultiFunction.hh"

namespace Amanzi {

MultiFunction::MultiFunction(
  const std::vector<Teuchos::RCP<const Function>>& functions)
  : functions_(functions)
{
  Kokkos::resize(values_, functions_.size());
  // values_ = new double[functions_.size()];
};

MultiFunction::MultiFunction(const Teuchos::RCP<const Function>& function)
  : functions_(1, function)
{
  Kokkos::resize(values_, 1);
  // values_ = new double[1];
};


MultiFunction::MultiFunction(Teuchos::ParameterList& plist)
{
  FunctionFactory factory;

  if (plist.isParameter("number of dofs")) {
    if (plist.isType<int>("number of dofs")) {
      int ndofs = plist.get<int>("number of dofs");

      if (ndofs < 1) {
        // ERROR -- invalid number of dofs
        AMANZI_ASSERT(0);
      }

      for (int lcv = 1; lcv != (ndofs + 1); ++lcv) {
        std::stringstream sublist_name;
        sublist_name << "dof " << lcv << " function";
        functions_.push_back(
          Teuchos::rcp(factory.Create(plist.sublist(sublist_name.str()))));
      }
    } else {
      // ERROR -- invalid number of dofs
      AMANZI_ASSERT(0);
    }
  } else {
    // assume it is a single dof function
    functions_.push_back(Teuchos::rcp(factory.Create(plist)));
  };

  // values_ = new double[functions_.size()];
  Kokkos::resize(values_, functions_.size());
}


MultiFunction::~MultiFunction(){
  // delete [] values_;
};


int
MultiFunction::size() const
{
  return functions_.size();
};


Kokkos::View<double*,Kokkos::HostSpace>
MultiFunction::operator()(const Kokkos::View<double*,Kokkos::HostSpace>& xt) const
{
  for (int i = 0; i != size(); ++i) { values_.view_host()[i] = (*functions_[i])(xt); }
  return values_.view_host();
};


} // namespace Amanzi
