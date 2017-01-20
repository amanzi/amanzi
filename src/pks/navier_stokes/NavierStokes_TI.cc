/*
  Navier Stokes PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "NavierStokes_PK.hh"

namespace Amanzi {
namespace NavierStokes {

/* ******************************************************************
* Calculate f(u, du/dt) = a d(s(u))/dt + A*u - rhs.
****************************************************************** */
void NavierStokes_PK::Functional(double t_old, double t_new, 
                                 Teuchos::RCP<TreeVector> u_old, Teuchos::RCP<TreeVector> u_new, 
                                 Teuchos::RCP<TreeVector> f)
{ 
}

}  // namespace NavierStokes
}  // namespace Amanzi
 
