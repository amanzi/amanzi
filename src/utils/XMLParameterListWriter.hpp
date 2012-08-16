#ifndef XMLPARAMETERLISTWRITER_H
#define XMLPARAMETERLISTWRITER_H

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLObject.hpp"
#include "Teuchos_Utils.hpp"

namespace Teuchos
{
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

