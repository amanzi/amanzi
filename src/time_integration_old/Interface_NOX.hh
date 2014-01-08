/*
This is the flow component of the Amanzi code. 
License: BSD
Authors: Neil Carlson (version 1) 
         Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
         Daniil Svyatskiy (version 3) (dasvyat@lanl.gov)
*/

#ifndef __INTERFACE_NOX_HPP__
#define __INTERFACE_NOX_HPP__


#include "exceptions.hh"
#include "NOX_Epetra_Interface_Required.H"
#include "NOX_Epetra_Interface_Jacobian.H"
#include "NOX_Epetra_Interface_Preconditioner.H"
#include "NOX_StatusTest_Generic.H"
#include "BDF1_TI.hh"
#include "Epetra_Map.h"
#include "Epetra_Operator.h"
#include "Mesh.hh"
#include "BDF1_State.hh"


class BDF1State;

namespace Amanzi {
namespace AmanziFlow {

class Interface_NOX : public NOX::Epetra::Interface::Required,
                      public NOX::Epetra::Interface::Jacobian,
                      public NOX::Epetra::Interface::Preconditioner {
 public:
  Interface_NOX(BDF2::fnBase* FPK_, BDF1State* state_, const Epetra_Vector& uprev, double time_, double dt) : 
		FPK(FPK_), u0(uprev), lag_prec_(3), lag_count_(0) {
                        state = state_;
			time = time_;
			deltaT = dt;
			fun_eval = 0;
			fun_eval_time = 0;
		}
  ~Interface_NOX() {}

  // required interface members
  bool computeF(const Epetra_Vector& x, Epetra_Vector& f, FillType flag);
  bool computeJacobian(const Epetra_Vector& x, Epetra_Operator& J) { assert(false); }
  bool computePreconditioner(const Epetra_Vector& x, Epetra_Operator& M, Teuchos::ParameterList* params=NULL);
  void printTime();

  inline void setPrecLag(int lag_prec) { lag_prec_ = lag_prec;}
  inline void resetPrecLagCounter() { lag_count_ = 0; }
  inline int getPrecLag() const { return lag_prec_; }
  inline int getPrecLagCounter() const { return lag_count_; }

 private:
  BDF2::fnBase* FPK;
  BDF1State *state;
  const Epetra_Vector& u0;	// value at the previous time step

  double deltaT, time;		// time step
  int lag_prec_;  // the preconditioner is lagged this many times before it is recomputed
  int lag_count_; // this counts how many times the preconditioner has been lagged
  int fun_eval;
  double  fun_eval_time;
};

/// class for precondtioning JFNK solver
class JFNK_Preconditioner : public Epetra_Operator {
	public:			    
                /// constructor
		JFNK_Preconditioner(BDF2::fnBase* FPK_
                            ):
                           FPK(FPK_) {mesh = FPK->mesh();}
                 /// destructor
		~JFNK_Preconditioner() {}
		
		 /// required methods
	int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const { Y = X; return 0;}
	int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {return FPK->ApllyPrecInverse(X,Y); }
	bool UseTranspose() const { return false; }
	int SetUseTranspose(bool trans) { return 1; }
	
	const Epetra_Comm& Comm() const { return *(mesh->get_comm()); }
	const Epetra_Map& OperatorDomainMap() const { return FPK->super_map(); }
	const Epetra_Map& OperatorRangeMap() const { return FPK->super_map(); }

	const char* Label() const { return strdup("Preconditioner Test"); }
	double NormInf() const { return 0.0; }
	bool HasNormInf() const { return false; }
	
	private:
// 		Teuchos::RCP<Flow_State> FS;
                BDF2::fnBase* FPK;
		Teuchos::RCP<const AmanziMesh::Mesh> mesh;
//		Epetra_Map map;
};


class PK_enorm : public NOX::StatusTest::Generic {
        
        public:
         //! Constructor
                PK_enorm(BDF2::fnBase* FPK_, double tolerance): FPK(FPK_) { specifiedTolerance = tolerance;}
         //! Destructor.
                ~PK_enorm() {}

        // derived
                NOX::StatusTest::StatusType 
                        checkStatus(const NOX::Solver::Generic& problem,
                                NOX::StatusTest::CheckType checkType);

        // derived
                NOX::StatusTest::StatusType getStatus() const;

                ostream& print(ostream& stream, int indent = 0) const;
        
        private:
                BDF2::fnBase* FPK;
                
                 //! Vector containing the update for the current outer iteration 
                Teuchos::RCP<NOX::Abstract::Vector> updateVectorPtr;
                
                 //! %Status
                NOX::StatusTest::StatusType status;

                //! Tolerance required for convergence.
                double specifiedTolerance;

                //! Initial tolerance
                double initialTolerance;
                
                //! Norm of Error
                
                double normF;

                //! Ostream used to print errors
                NOX::Utils utils;
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
