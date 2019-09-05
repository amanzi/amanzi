#include <math.h>

#include "FunctionStandardMath.hh"
#include "errors.hh"

namespace Amanzi {

FunctionStandardMath::FunctionStandardMath(char op[10],
        double amplitude, double parameter, double shift) :
    amplitude_(amplitude), parameter_(parameter), shift_(shift) {

  if (strcmp(op,"cos")==0){
    op_ = COS; 
  }else if (strcmp(op,"sin")==0){
    op_ = SIN; 
  }else if (strcmp(op,"tan")==0){
    op_ = TAN; 
  }else if (strcmp(op,"acos")==0){
    op_ = ACOS; 
  }else if (strcmp(op,"asin")==0){
    op_ = ASIN; 
  }else if (strcmp(op,"atan")==0){
    op_ = ATAN; 
  }else if (strcmp(op,"cosh")==0){
    op_ = COSH; 
  }else if (strcmp(op,"sinh")==0){
    op_ = SINH; 
  }else if (strcmp(op,"tanh")==0){
    op_ = TANH; 
  }else if (strcmp(op,"exp")==0){
    op_ = EXP; 
  }else if (strcmp(op,"log")==0){
    op_ = LOG; 
  }else if (strcmp(op,"log10")==0){
    op_ = LOG10; 
  }else if (strcmp(op,"sqrt")==0){
    op_ = SQRT; 
  }else if (strcmp(op,"ceil")==0){
    op_ = CEIL; 
  }else if (strcmp(op,"fabs")==0){
    op_ = FABS; 
  }else if (strcmp(op,"floor")==0){
    op_ = FLOOR; 
  }else if (strcmp(op,"mod")==0){
    op_ = MOD; 
  }else if (strcmp(op,"pow")==0){
    op_ = POW; 
  }else{
    std::stringstream m;
    m << "Invalid or unknown standard math function " << op;
    Errors::Message message(m.str());
    Exceptions::amanzi_throw(message);
  }
}


double FunctionStandardMath::operator()(const Kokkos::View<double*>& x) const
{
  double x0 = x[0] - shift_;
  switch(op_){
    case COS:
      return amplitude_ * cos(parameter_ * x0);
      break; 
    case SIN: 
      return amplitude_ * sin(parameter_ * x0);
      break; 
    case TAN:  
      return amplitude_ * tan(parameter_ * x0);
      break; 
    case ACOS: 
      return amplitude_ * acos(parameter_ * x0);
      break; 
    case ASIN: 
      return amplitude_ * asin(parameter_ * x0);
      break; 
    case ATAN:
      return amplitude_ * atan(parameter_ * x0);
      break; 
    case COSH: 
      return amplitude_ * cosh(parameter_ * x0);
      break; 
    case SINH:
      return amplitude_ * sinh(parameter_ * x0);
      break; 
    case TANH: 
      return amplitude_ * tanh(parameter_ * x0);
      break; 
    case EXP: 
      return amplitude_ * exp(parameter_ * x0);
      break; 
    case LOG: 
      if (x0 <= 0) InvalidDomainError_(x[0]);
      return amplitude_ * log(parameter_ * x0);
      break; 
    case LOG10: 
      if (x0 <= 0) InvalidDomainError_(x[0]);
      return amplitude_ * log10(parameter_ * x0);
      break; 
    case SQRT: 
      if (x0 < 0) InvalidDomainError_(x[0]);
      return amplitude_ * sqrt(parameter_ * x0);
      break; 
    case CEIL: 
      return amplitude_ * ceil(x0);
      break; 
    case FABS: 
      return amplitude_ * fabs(x0);
      break; 
    case FLOOR: 
      return amplitude_ * floor(x0);
      break; 
    case POW:
      return amplitude_ * pow(x0, parameter_);
      break; 
    case MOD: 
      return fmod(x0, parameter_);
      break; 
    default: 
      std::stringstream m;
      m << "Invalid or unknown standard math function " << op_;
      Errors::Message message(m.str());
      Exceptions::amanzi_throw(message);
  } 
  return 0.0;
}


void FunctionStandardMath::InvalidDomainError_(double x) const {
  std::stringstream m;
  m << "Value " << x << " is not in the domain of operator " << op_ << ".";
  Errors::Message message(m.str());
  Exceptions::amanzi_throw(message);
}

} // namespace Amanzi
