#ifndef RICHARDS_MODEL_EVALUATOR_HPP
#define RICHARDS_MODEL_EVALUATOR_HPP

#include "Teuchos_ParameterList.hpp"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Comm.h"
#include "Epetra_CrsGraph.h"

#include "BDF2_fnBase.hpp"
#include "RichardsProblem.hpp"

class RichardsModelEvaluator : public BDF2::fnBase {
public:

  // Constructor
  RichardsModelEvaluator(RichardsProblem *problem, 
			 Teuchos::RCP<DiffusionMatrix> &matrix,
			 Teuchos::ParameterList &plist, 
			 const Epetra_Map &map); 

  // Initialization
  void initialize(Teuchos::RCP<Epetra_Comm> &epetra_comm_ptr, Teuchos::ParameterList &params);


  void fun(const double t, const Epetra_Vector& u, const Epetra_Vector& udot, Epetra_Vector& f);
  void precon(const Epetra_Vector& u, Epetra_Vector& Pu);
  double enorm(const Epetra_Vector& u, const Epetra_Vector& du);
  void update_precon(const double t, const Epetra_Vector& up, const double h, int& errc);

  bool is_admissible(const Epetra_Vector& up);

private:
  
  RichardsProblem* problem_;

  // The diffusion matrix
  Teuchos::RCP<DiffusionMatrix> D;

  Epetra_Map map_;

  // ML preconditioner for the Schur complement system.
  ML_Epetra::MultiLevelPreconditioner *MLprec;

  Teuchos::ParameterList ML_plist;

};

#endif // RICHARDS_MODEL_EVALUATOR_HPP 
