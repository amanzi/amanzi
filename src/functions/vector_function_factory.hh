/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $AMANZI_DIR/COPYRIGHT
Author Ethan Coon

Factory for functions from R^d to R^n.

See ATS_DIR/src/factory/factory.hh for an in-depth discussion of this
approach.

USAGE:

------------------------------------------------------------------------- */

#ifndef ATS_FACTORY_BY_FUNCTION_HH_
#define ATS_FACTORY_BY_FUNCTION_HH_

#include <map>
#include <string>
#include "Teuchos_ParameterList.hpp"

#include "errors.hh"

#include "vector_function.hh"

namespace Amanzi {

class VectorFunctionFactory {

public:
  typedef std::map<std::string, Teuchos::RCP<VectorFunction> (*)(Teuchos::ParameterList&)> map_type;

  static Teuchos::RCP<VectorFunction> Create(Teuchos::ParameterList& plist) {
    // Pull the type.  Note that this defaults to composite function,
    // which handles the (old-form) case of a single function on a
    // single dof with no meta-data allowed (or it breaks the
    // function-factory).
    std::string s;
    if (plist.isParameter("Function type")) {
      s = plist.get<std::string>("Function type");
    } else {
      s = "composite function";
    }

    typename map_type::iterator iter = GetMap()->find(s);
    if (iter == GetMap()->end()) {
      Errors::Message m;
      m << "Factory: unknown type: " << s.c_str();
      Exceptions::amanzi_throw(m);
    }
    return iter->second(plist);
  }

protected:
  static map_type* GetMap() {
    if (!map_) {
      map_ = new map_type;
    }
    return map_;
  }

private:
  static map_type* map_;
};


class RegisteredVectorFunctionFactory : public VectorFunctionFactory {
public:
  // Constructor for the registered factory.  Needs some error checking in
  // case a name s is already in the map? (i.e. two implementations trying to
  // call themselves the same thing) --etc
  typedef Teuchos::RCP<VectorFunction> (*FactoryFunction)(Teuchos::ParameterList&);

  RegisteredVectorFunctionFactory(const std::string& s, FactoryFunction functor) {
    VectorFunctionFactory::GetMap()->insert(std::pair<std::string, FactoryFunction>(s,functor));
  }
};

} // namespace


#endif

