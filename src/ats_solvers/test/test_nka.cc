#include <iostream>
#include "UnitTest++.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Epetra_MpiComm.h"
#include "Epetra_SerialComm.h"

#include "Epetra_BlockMap.h"
#include "Epetra_Vector.h"

#include "TreeVector.hh"
#include "NonlinarKrylovAccelerator.hh"

TEST(NKA_NONLINEAR) {

  std::cout << "NKA_NONLINEAR..." << std::endl;

  Epetra_Comm *comm = new Epetra_SerialComm(); 
  // create an Epetra_Vector
  Epetra_BlockMap map(3,1,0,*comm);
  
  Teuchos::RCP<Epetra_Vector> x = Teuchos::rcp(new Epetra_Vector(map));
  Teuchos::RCP<Amanzi::TreeVector>  x_tree = Teuchos::rcp(new Amanzi::TreeVector(std::string("x vec"),x));
  
  Amanzi::NonlinearKrylovAccelerator FPA(10,0.0,*x_tree);
  FPA.nka_restart();
  FPA.nka_relax();

  // set up a simple fixed point iteration...
  
  Teuchos::RCP<Epetra_Vector> f = Teuchos::rcp(new Epetra_Vector(map));
  Teuchos::RCP<Amanzi::TreeVector>  f_tree = Teuchos::rcp(new Amanzi::TreeVector(std::string("f vec"),f));
  
  Teuchos::RCP<Epetra_Vector> dx = Teuchos::rcp(new Epetra_Vector(map));
  Teuchos::RCP<Amanzi::TreeVector>  dx_tree = Teuchos::rcp(new Amanzi::TreeVector(std::string("dx vec"),dx));
  

  // initial value
  Epetra_Vector* fe = (* (*f_tree)[0])(0);
  Epetra_Vector* xe = (* (*x_tree)[0])(0);

  (*xe)[0] = 0.0; 
  (*xe)[1] = 0.0; 
  (*xe)[2] = 0.0;

  double norm(0.0);
  int nka_iterations(0);
  do {
    nka_iterations++;      
    
    // function evaluation f <-- f(x^n)
    
    (*fe)[0] = (*xe)[0] - (cos((*xe)[0]) - sin((*xe)[1])) / 3.0;
    (*fe)[1] = (*xe)[1] - (cos((*xe)[0]) - 2.0*sin((*xe)[1])) / 3.0;
    (*fe)[2] = (*xe)[2] - (cos((*xe)[0]) - 3.0*sin((*xe)[2])) / 3.0;

    // compute the NKA correction  dx <-- NKA(f)
    FPA.nka_correction(*f_tree, *dx_tree);
    
    // do the fixed point step   x^{n+1} <-- x^n - dx
    x_tree->Update(-1.0, *dx_tree, 1.0);
    
    dx_tree->NormInf(&norm);

    std::cout << std::setprecision(15) << std::scientific;

    std::cout << "iterate = " << nka_iterations;
    std::cout << ", error = " << norm << std::endl;
  } while (norm > 1e-14);
  std::cout << "final error = " << norm << std::endl;

  CHECK_EQUAL(nka_iterations,9);  

  delete comm;
  std::cout << "NKA_NONLINEAR... DONE." << std::endl;
}

TEST(NKA_LINEAR) {

  // This is a linear problem, NKA solves this in N steps where
  // N is the number of unknowns.

  std::cout << "NKA_LINEAR..." << std::endl;

  Epetra_Comm *comm = new Epetra_SerialComm(); 
  // create an Epetra_Vector
  Epetra_BlockMap map(3,1,0,*comm);
  
  Teuchos::RCP<Epetra_Vector> x = Teuchos::rcp(new Epetra_Vector(map));
  Teuchos::RCP<Amanzi::TreeVector>  x_tree = Teuchos::rcp(new Amanzi::TreeVector(std::string("x vec"),x));
  
  Amanzi::NonlinearKrylovAccelerator FPA(3,0.01,*x_tree);
  FPA.nka_restart();

  // set up a simple fixed point iteration...
  
  Teuchos::RCP<Epetra_Vector> f = Teuchos::rcp(new Epetra_Vector(map));
  Teuchos::RCP<Amanzi::TreeVector>  f_tree = Teuchos::rcp(new Amanzi::TreeVector(std::string("f vec"),f));
  
  Teuchos::RCP<Epetra_Vector> dx = Teuchos::rcp(new Epetra_Vector(map));
  Teuchos::RCP<Amanzi::TreeVector>  dx_tree = Teuchos::rcp(new Amanzi::TreeVector(std::string("dx vec"),dx));
  

  // initial value
  Epetra_Vector* fe = (* (*f_tree)[0])(0);
  Epetra_Vector* xe = (* (*x_tree)[0])(0);

  (*xe)[0] = 1.0;
  (*xe)[1] = 2.0; 
  (*xe)[2] = -1.0;   

  double norm(0.0);
  int nka_iterations(0);
  do {
    nka_iterations++;      
    // function evaluation f <-- f(x^n)
    
    (*fe)[0] = 1.0 - (3.0*(*xe)[0]      -(*xe)[1]      -(*xe)[2]);
    (*fe)[1] = 2.0 - (   -(*xe)[0] + 3.0*(*xe)[1]      -(*xe)[2]);
    (*fe)[2] = -1.0 -(   -(*xe)[0]      -(*xe)[1] + 3.0*(*xe)[2]);

    // compute the NKA correction  dx <-- NKA(f)
    FPA.nka_correction(*f_tree, *dx_tree);

    // do the fixed point step   x^{n+1} <-- x^n - dx
    x_tree->Update(-1.0, *dx_tree, 1.0);
    
    dx_tree->NormInf(&norm);

    std::cout << std::setprecision(15) << std::scientific;

    std::cout << "iterate = " << nka_iterations;
    std::cout << ", error = " << norm << std::endl;
    
  } while (norm > 1e-14);

  std::cout << "final error = " << norm << std::endl;

  CHECK_EQUAL(nka_iterations,3);  
  
  delete comm;
  std::cout << "NKA_LINEAR... DONE." << std::endl;

}
