#ifndef __Darcy_problem__
#define __Darcy_problem__

#include "Epetra_Vector.h"
#include "Mesh_maps.hh"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_Map.h"
#include "Flow_State.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_RowMatrix.h"
#include "mimetic_hex.hpp"

class Darcy_problem {

public:

  Darcy_problem(Teuchos::RCP<Flow_State> FS_ ):
    FS(FS_) {};
  ~Darcy_problem() {};

  void ComputeF(const Epetra_Vector & x, Epetra_Vector & f);

  const Teuchos::RCP<Epetra_Map> get_NL_map () const { return NL_map; };
  const Teuchos::RCP<Epetra_CrsMatrix> get_PrecMat () const { return PrecMat; };

private:
  const Teuchos::RCP<Flow_State> FS;
  const Teuchos::RCP<Epetra_Map> NL_map;
  const Teuchos::RCP<Epetra_CrsMatrix> PrecMat;
  mimetic_hex MD[];
  double K[]; // array of (scalar) diffusion coefficients on cells

};

#endif
