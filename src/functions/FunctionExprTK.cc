#include "FunctionExprTK.hh"
#include "errors.hh"

namespace Amanzi {

FunctionExprTK::FunctionExprTK(int n, const std::string& formula)
{
  if (n < 1 || n > 4) {
    Errors::Message m;
    m << "The number of pameters is incorrect";
    Exceptions::amanzi_throw(m);
  }

  exprtk_ = std::make_shared<Utils::ExprTK>();
  if (!exprtk_->Initialize(n, formula)) {
    Errors::Message m;
    m << "Initialization of expression has failed.";
    Exceptions::amanzi_throw(m);
  }
}


double FunctionExprTK::operator()(const std::vector<double>& args) const
{
  return (*exprtk_)(args);
}

}  // namespace Amanzi
