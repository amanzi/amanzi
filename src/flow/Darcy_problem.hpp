#ifndef __Darcy_problem__
#define __Darcy_problem__

#include "Epetra_Vector.h"
#include "Mesh_maps_simple.hh"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_Map.h"
#include "Flow_State.hpp"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_RowMatrix.h"
#include "mimetic_hex.hpp"
#include "Flow_BCs.hpp"

class Darcy_problem {

public:

  Darcy_problem(Teuchos::RCP<Flow_State> FS_ ):
    FS(FS_) {};
  ~Darcy_problem() {};

  void ComputeF(const Epetra_Vector & x, Epetra_Vector & f);

  const Teuchos::RCP<Epetra_Map> get_NL_map () const { return NL_map; };
  const Teuchos::RCP<Epetra_CrsMatrix> get_PrecMat () const { return PrecMat; };

  void initialize();

private:

  Teuchos::RCP<Flow_State> FS;
  Teuchos::RCP<Epetra_Map> NL_map;
  Teuchos::RCP<Epetra_FECrsMatrix> PrecMat;
  Teuchos::RCP<Flow_BCs> FBC;
  std::vector<mimetic_hex> MD;
  
  double rho_; // constant fluid density
  std::vector<double> K_; // (scalar) diffusion coefficients on cells
  std::vector<double> area_; // face areas
  
  // Private stuff related to boundary conditions.
  
  enum bc_types {
    PRESSURE_CONSTANT = 1,
    NO_FLOW,
    DARCY_CONSTANT
  };
    
  struct bc_spec {
    bc_types type;
    int num_faces;
    std::vector<int> faces;
    std::vector<double> aux;
    double value;
  };
  std::vector<struct bc_spec> bc_;
  
  void BC_setup (std::vector<flow_bc> & bc);
  void FBC_initial_pass(double p_face[]);
  void FBC_final_pass(double f_face[]);
  
};

#endif
