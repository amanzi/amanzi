/* -------------------------------------------------------------------------
ATS & Amanzi

License: see $AMANZI_DIR/COPYRIGHT
Author Ethan Coon

Function from R^d to R^n.
------------------------------------------------------------------------- */

#ifndef AMANZI_MULTIVECTOR_FUNCTION_HH_
#define AMANZI_MULTIVECTOR_FUNCTION_HH_

#include <vector>
#include "Teuchos_RCP.hpp"
#include "Function.hh"
#include "FunctionFactory.hh"

namespace Amanzi {

class MultiFunction {

public:
  MultiFunction(const std::vector<Teuchos::RCP<const Function> >& functions);
  MultiFunction(const Teuchos::RCP<const Function>& function);
  MultiFunction(Teuchos::ParameterList& plist);

  ~MultiFunction();

  int size() const;
  Kokkos::View<double*> operator()(const Kokkos::View<double*>& xt) const;

  void apply(const Kokkos::View<double**>& in, Kokkos::View<double**>& out) const {
    for(int i = 0 ; i < size(); ++i){
      Kokkos::View<double*> out_i = Kokkos::subview(out,Kokkos::ALL,i); 
      functions_[i]->apply(in,out_i); 
    }
  }

 private:
  std::vector<Teuchos::RCP<const Function> > functions_;
  Kokkos::View<double*> values_;
};

} // namespace

#endif
