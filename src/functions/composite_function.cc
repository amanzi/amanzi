/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $AMANZI_DIR/COPYRIGHT
Author Ethan Coon

Factory for vector functions which are composed of multiple scalar functions.

The expected plist is of the form:

<ParameterList name="constuctor plist">
  <Parameter name="Number of DoFs">

  <ParameterList name="Function 1">
    <ParameterList name="function-constant">
      ...
    </ParameterList>
  </ParameterList>

  <ParameterList name="Function 2">
    <ParameterList name="function-linear">
      ...
    </ParameterList>
  </ParameterList>

  ...
</ParameterList>


Where each of the "Function X" lists are valid input to the
function-factory Create() method (see ./function-factory.hh).

------------------------------------------------------------------------- */

#include "composite_function.hh"

namespace Amanzi {

// Register the factory
RegisteredVectorFunctionFactory CompositeFunction::reg_("composite function", &CreateCompositeFunction);

CompositeFunction::CompositeFunction(
        const std::vector<Teuchos::RCP<const Function> >& functions) :
    functions_(functions) {
  values_ = new double[functions.size()];
};

CompositeFunction::CompositeFunction(const Teuchos::RCP<const Function>& function) :
    functions_(1, function)  {
values_ = new double[1];
};

CompositeFunction::~CompositeFunction() {
  delete [] values_;
};

VectorFunction* CompositeFunction::Clone() const {
  return new CompositeFunction(functions_);
};


int CompositeFunction::size() const {
  return functions_.size();
};


double* CompositeFunction::operator() (const double* xt) const {
  for (int i=0; i!=size(); ++i) {
    values_[i] = (*functions_[i])(xt);
  }
  return values_;
};


Teuchos::RCP<VectorFunction>
CreateCompositeFunction(Teuchos::ParameterList& plist) {
  std::vector<Teuchos::RCP<const Function> > functions;
  FunctionFactory factory;

  if (plist.isParameter("Number of DoFs")) {
    if (plist.isType<int>("Number of DoFs")) {
      int ndofs = plist.get<int>("Number of DoFs");

      if (ndofs < 1) {
        // ERROR -- invalid number of dofs
      }

      for (int lcv=1; lcv!=(ndofs+1); ++lcv) {
        std::stringstream sublist_name;
        sublist_name << "DoF " << lcv << " Function";
        functions.push_back(Teuchos::rcp(factory.Create(plist.sublist(sublist_name.str()))));
      }
    } else {
      // ERROR -- invalid number of dofs
    }
  } else {
    // assume it is a single dof function
    functions.push_back(Teuchos::rcp(factory.Create(plist)));
  };
  return Teuchos::rcp(new CompositeFunction(functions));
}


} // namespace

