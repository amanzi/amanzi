#include <iostream>
#include "UnitTest++.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Epetra_MpiComm.h"
#include "Epetra_SerialComm.h"

#include "Epetra_BlockMap.h"
#include "Epetra_Vector.h"

#include "BDF2_SolutionHistory.hpp"


TEST(SolutionHistory_1) {

  Epetra_Comm *comm = new Epetra_SerialComm();

  // create an Epetra_Vector
  Epetra_BlockMap map(10,1,0,*comm);
  
  Epetra_Vector x(map);

  
  x.PutScalar(0.0);
  // create a solution history of size three
  BDF2::SolutionHistory SH(3, 0.0, x);

  CHECK_EQUAL(SH.history_size(),1);
  CHECK_EQUAL(SH.most_recent_time(),0.0);
  
  x.PutScalar(1.0);
  SH.record_solution(1.0,x);

  CHECK_EQUAL(SH.history_size(),2);
  CHECK_EQUAL(SH.most_recent_time(),1.0);

  
  x.PutScalar(4.0);
  SH.record_solution(2.0,x);  

  CHECK_EQUAL(SH.history_size(),3);
  CHECK_EQUAL(SH.most_recent_time(),2.0);
  
  x.PutScalar(9.0);
  SH.record_solution(3.0,x);  

  CHECK_EQUAL(SH.history_size(),3);
  CHECK_EQUAL(SH.most_recent_time(),3.0);

  // check that the most recent vector in fact is the
  // one that contains  9.0
  Epetra_Vector y(map);
  SH.most_recent_solution(y);
  x.PutScalar(9.0);
  y.Update(-1.0,x,1.0); // compute the difference between y and x 
  double *norminf = new double[1];
  y.NormInf(norminf); // compute the max norm of the computed difference
  CHECK_EQUAL(norminf[0],0.0);

  
  std::vector<double> h;
  SH.time_deltas(h);
  
  CHECK_EQUAL(h.size(),SH.history_size()-1);
  CHECK_EQUAL(h[0], 1.0);
  CHECK_EQUAL(h[1], 2.0);

  
  // interpolate 1st order
  SH.interpolate_solution(2.5, y, 1);
  
  // we should get 6.5
  x.PutScalar(6.5);
  x.Update(-1.0,y,1.0);
  x.NormInf(norminf);
  CHECK_EQUAL(norminf[0],0.0);
  
  // interpolate 2nd order
  SH.interpolate_solution(2.5, y, 2);
  
  // from the quadratic through (1,1),(2,4),(3,9)
  // we should get 6.25
  x.PutScalar(6.25);
  x.Update(-1.0,y,1.0);
  x.NormInf(norminf);
  CHECK_EQUAL(norminf[0],0.0);  

  // interpolate maximum order (should be 2nd)
  SH.interpolate_solution(2.5, y);
  
  // from the quadratic through (1,1),(2,4),(3,9)
  // we should get 6.25
  x.PutScalar(6.25);
  x.Update(-1.0,y,1.0);
  x.NormInf(norminf);
  CHECK_EQUAL(norminf[0],0.0);  
}


TEST(SolutionHistory_2) {

  Epetra_Comm *comm = new Epetra_SerialComm();

  // create an Epetra_Vector
  Epetra_BlockMap map(10,1,0,*comm);
  
  Epetra_Vector x(map);
  Epetra_Vector xdot(map);
  
  x.PutScalar(0.0);
  xdot.PutScalar(0.0);
  // create a solution history of size three
  BDF2::SolutionHistory SH(4, 0.0, x, xdot);

  x.PutScalar(1.0);
  xdot.PutScalar(2.0);
  SH.record_solution(1.0,x,xdot);
  
  Epetra_Vector y(map);
  
  // interpolate 1st order
  SH.interpolate_solution(0.5, y, 2);

  double *norminf = new double[1];
  x.PutScalar(0.25);
  x.Update(-1.0,y,1.0);
  x.NormInf(norminf);
  CHECK_EQUAL(norminf[0],0.0);    

  
  // interpolate maximum order (3rd in this case)
  SH.interpolate_solution(0.5, y);

  x.PutScalar(0.25);
  x.Update(-1.0,y,1.0);
  x.NormInf(norminf);
  CHECK_EQUAL(norminf[0],0.0);    
  

}
