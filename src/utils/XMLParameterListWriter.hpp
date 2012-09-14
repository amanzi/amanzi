#ifndef XMLPARAMETERLISTWRITER_H
#define XMLPARAMETERLISTWRITER_H

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLObject.hpp"
#include "Teuchos_Utils.hpp"
#include "Teuchos_ArrayView.hpp"

namespace Teuchos
{

#if 0
  template<>
  class ToStringTraits<double> {
  public:
    static std::string toString (const double& t) {
      std::ostringstream os;
      os.setf (std::ios::scientific);
      os.precision (17);
      os << t;
      return os.str();
    }
  };
#endif

  inline std::string Amanzi_toString(const Teuchos::Array<double>& a) {
    std::ostringstream ss;
    
    ss.setf (std::ios::scientific);
    ss.precision (17);

    ss << "{";
    
    for (int i=0; i < a.size(); ++i) {
      ss << a[i];
      if (i < a.size()-1) ss << ", ";
    }
    ss << "}";
    
    return ss.str();
  }
    

  
  class  Amanzi_XMLParameterListWriter
  {
  public:
    Amanzi_XMLParameterListWriter();
    XMLObject toXML(const ParameterList& p) const ;
    
  private:
    
    XMLObject toXML(const ParameterEntry& p) const ;
  };
}
#endif

