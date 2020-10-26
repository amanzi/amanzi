/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (amklinv@sandia.gov)
           Ethan Coon (coonet@ornl.gov)
*/
//!  Direct solvers via Trilinos.
/*!

.. warning:: undocumented

*/

#ifndef  AMANZI_AMESOS2_OPERATOR_HH_
#define  AMANZI_AMESOS2_OPERATOR_HH_

#include <cmath>

#include "Amesos2.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

//#include "CompositeMatrix.hh"
#include "errors.hh"
#include "VerboseObject.hh"

#include "Inverse.hh"
#include "InverseDefs.hh"

namespace Amanzi {

class VerboseObject;

namespace AmanziSolvers {

/* ******************************************************************
* Auxiliary base class.
****************************************************************** */
class DirectMethodAmesos2 :
      public Inverse<Matrix_type,Matrix_type,Vector_type,Map_type>
{
 private:
  using Inv = Inverse<Matrix_type,Matrix_type,Vector_type,Map_type>;
  
 public:
  DirectMethodAmesos2() :
      Inv(),
      inited_(false),
      updated_(false),
      computed_(false)
  {};


  virtual void set_inverse_parameters(Teuchos::ParameterList& plist) override;
  virtual void initializeInverse() override;
  virtual void computeInverse() override;
  virtual int applyInverse(const Vector_type&, Vector_type&) const override;
  virtual void update(const Teuchos::RCP<Matrix_type>&) override final {}; 

  virtual int returned_code() const override { return returned_code_; }
  virtual std::string returned_code_string() const override { return "success"; }

 protected:
  using Inv::m_;
  using Inv::h_;  

  Teuchos::ParameterList plist_;
  Teuchos::RCP<VerboseObject> vo_;

  Teuchos::RCP<Amesos2::Solver<Matrix_type,MultiVector_type> > solver_;
  
  std::string solver_name_;
  mutable int returned_code_;

  bool inited_;
  bool updated_;
  bool computed_;
};


}  // namespace AmanziSolvers
}  // namespace Amanzi

#endif
