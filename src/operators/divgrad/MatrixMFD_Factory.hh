/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

   ATS
   Author: Ethan Coon

   MatrixMFD factory.

   See a more thorough factory discussion in $ATS_DIR/src/factory/factory.hh.

   Simplest usage:

   // pk_implementation.hh
   #include "pk.hh"
   class DerivedMatrixMFD : public MatrixMFD {
     DerivedMatrixMFD(Teuchos::ParameterList& plist,
               const Teuchos::RCP<TreeVector>& solution);
     ...
   private:
     static RegisteredMatrixMFDFactory<MatrixMFD,DerivedMatrixMFD> factory_; // my factory entry
     ...
   };

   ------------------------------------------------------------------------- */

#ifndef ATS_MatrixMFD_FACTORY_HH_
#define ATS_MatrixMFD_FACTORY_HH_

#include <iostream>
#include <map>
#include <string>
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Mesh.hh"
#include "MatrixMFD.hh"

namespace Amanzi {
namespace Operators {

class MatrixMFD_Factory {

public:
  typedef std::map<std::string, MatrixMFD* (*)(Teuchos::ParameterList&,
                     const Teuchos::RCP<const AmanziMesh::Mesh>&)> map_type;

  static Teuchos::RCP<MatrixMFD> 
  CreateMatrixMFD(Teuchos::ParameterList& plist,
                  const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) {
    std::string s = plist.get<string>("MatrixMFD type");
    typename map_type::iterator iter = GetMap()->find(s);
    if (iter == GetMap()->end()) {
      std::cout << "cannot get item of type: " << s << std::endl;

      for (typename map_type::iterator iter=GetMap()->begin();
           iter!=GetMap()->end(); ++iter) {
        std::cout << "  option: " << iter->first << std::endl;
      }
      return Teuchos::null;
    }
    return Teuchos::rcp(iter->second(plist, mesh));
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


template<typename T> MatrixMFD* 
CreateT(Teuchos::ParameterList& plist,
        const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) {
  return new T(plist, mesh);
}


template<typename T>
class RegisteredMatrixMFD_Factory : public MatrixMFD_Factory {
public:
  // Constructor for the registered factory.  Needs some error checking in
  // case a name s is already in the map? (i.e. two implementations trying to
  // call themselves the same thing) --etc
  RegisteredMatrixMFD_Factory(const std::string& s) {
    GetMap()->insert(std::pair<std::string,MatrixMFD* (*)(Teuchos::ParameterList&, const Teuchos::RCP<const AmanziMesh::Mesh>&)>(s, &CreateT<T>));
  }
};

} // namespace
} // namespace

#endif
