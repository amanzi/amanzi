#ifndef DATADEBUG_HH__
#define DATADEBUG_HH__

#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_VerboseObject.hpp"

#include "Mesh.hh"
#include "VerboseObject_objs.hh"

namespace Amanzi {

class DataDebug {
 public:
  explicit DataDebug(Teuchos::RCP<AmanziMesh::Mesh> mesh);
  ~DataDebug() {}

  void
  write_region_data(std::string& region_name, const Epetra_Vector& data, std::string& description);
  void write_region_statistics(std::string& region_name,
                               const Epetra_Vector& data,
                               std::string& description);

 private:
  Teuchos::RCP<AmanziMesh::Mesh> mesh_;
};

} // namespace Amanzi

#endif
