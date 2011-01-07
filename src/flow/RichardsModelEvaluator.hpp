#ifndef RICHARDS_MODEL_EVALUATOR_HPP
#define RICHARDS_MODEL_EVALUATOR_HPP

#include "EpetraExt_ModelEvaluator.h"
#include "Teuchos_ParameterList.hpp"
#include "Rythmos_ConfigDefs.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Comm.h"
#include "Epetra_CrsGraph.h"


class RichardsModelEvaluator : public EpetraExt::ModelEvaluator {
public:

  // Constructor
  RichardsModelEvaluator(Teuchos::RCP<Epetra_Comm> &epetra_comm_ptr, Teuchos::ParameterList &params);

  // Initialization
  void initialize(Teuchos::RCP<Epetra_Comm> &epetra_comm_ptr, Teuchos::ParameterList &params);

  Teuchos::RCP<const Epetra_Map> get_x_map() const;

  Teuchos::RCP<const Epetra_Map> get_f_map() const;

  Teuchos::RCP<const Epetra_Vector> get_x_init() const;

  Teuchos::RCP<const Epetra_Vector> get_x_dot_init() const;

  Teuchos::RCP<Epetra_Operator> create_W() const;

  InArgs createInArgs() const;

  OutArgs createOutArgs() const;

  void evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const;

  Teuchos::RCP<Epetra_Vector> get_exact_solution( double t ) const;

private:

    // Epetra Comm:
    Teuchos::RCP<Epetra_Comm> epetra_comm_ptr_;
    // Epetra Map:
    Teuchos::RCP<const Epetra_Map> epetra_map_ptr_;
    
    // Global number of unknowns:
    int numElements_;

    Teuchos::RCP<Epetra_CrsGraph> W_graph_;

    // // This ModelEvaluator object will own a TransientInterface object
    // Teuchos::RCP<TransientInterface> problemInterfacePtr_;
};

#endif // RICHARDS_MODEL_EVALUATOR_HPP 
