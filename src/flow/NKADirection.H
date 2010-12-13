#ifndef __NKADIRECTION_H__
#define __NKADIRECTION_H__


#include "NOX_Common.H"

#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "Teuchos_ParameterList.hpp"
#include "NOX_Direction_Generic.H"
#include "NOX_GlobalData.H"
#include "NOX_Utils.H"

#include "NKA.H"


class NKADirection : public NOX::Direction::Generic {


public:

  NKADirection(const Teuchos::RCP<NOX::GlobalData>&, 
	       Teuchos::ParameterList&, const NOX::Abstract::Vector&);

  ~NKADirection();

  bool reset (const Teuchos::RCP<NOX::GlobalData>&, 
	      Teuchos::ParameterList&);
  
  bool compute (NOX::Abstract::Vector&, NOX::Abstract::Group&,
		const NOX::Solver::Generic&);

  bool compute (NOX::Abstract::Vector&, NOX::Abstract::Group&,
		const NOX::Solver::LineSearchBased&);

private:
  
  void throwError(const string&, const string&);

  // Printing Utils
  Teuchos::RCP<NOX::Utils> utils;

  // Global data pointer.  
  // Keep this so any stored parameter list remains valid.
  Teuchos::RCP<NOX::GlobalData> globalDataPtr;

  Teuchos::ParameterList *paramPtr; 

  nka *state;

  Teuchos::RCP<NOX::Abstract::Vector> tmpVecPtr;

};

#endif
