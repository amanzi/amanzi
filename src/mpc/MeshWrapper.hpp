#ifndef __MeshWrapper_hpp__
#define __MeshWrapper_hpp__

#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"
#include "STKMesh1D.hpp"
#include "STK1DDisc.hpp"
#include "DataLayout.hpp"


class MeshWrapper {

public:
  MeshWrapper( Teuchos::RCP<STKMesh1D>, Teuchos::RCP<STK1DDisc>, Teuchos::RCP<DataLayout>);
  ~MeshWrapper() {};

  Teuchos::RCP<const Epetra_Vector> get_element_volumes() const { return element_volumes; }

private:
  
  Teuchos::RCP<STKMesh1D> mesh1D;
  Teuchos::RCP<STK1DDisc> disc1D;

  Teuchos::RCP<Epetra_Vector> element_volumes;
  Teuchos::RCP<const DataLayout> data_layout;

};

#endif
