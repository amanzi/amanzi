/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
      Markus Berndt
*/

//! Visualization writes data and meshes to files.

/*
  Writes to file for visualization using generic Output object.
*/

#ifndef AMANZI_STATE_VISUALIZATION_HH_
#define AMANZI_STATE_VISUALIZATION_HH_

#include <string>
#include <vector>

#include "Teuchos_ParameterList.hpp"

#include "MeshDefs.hh"
#include "CompositeVector.hh"

#include "IOEvent.hh"
#include "Output.hh"

namespace Amanzi {

namespace AmanziMesh {
class Mesh;
}

class Visualization : public IOEvent {
 public:
  Visualization(Teuchos::ParameterList& plist,
                const Teuchos::RCP<AmanziMesh::Mesh>& mesh);

  void CreateFile(const double& time, const int& cycle);
  void FinalizeFile();

  // user-provided writing
  template <typename T>
  typename std::enable_if<!Output::writes<T>::value>::type
  Write(const Teuchos::ParameterList& attrs, const T& t) const
  {
    UserWriteVis(*this, attrs, t);
  }

  // output-provided writing
  template <typename T>
  typename std::enable_if<Output::writes<T>::value>::type
  Write(const Teuchos::ParameterList& attrs, const T& t) const
  {
    output_->Write(attrs, t);
  }

 protected:
  std::unique_ptr<Output> output_;
};

} // namespace Amanzi

#endif
