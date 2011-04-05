#include <iostream>
#include "UnitTest++.h"

#include "dbc.hh"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerbosityLevel.hpp"

#include "Epetra_MpiComm.h"
#include "Epetra_SerialComm.h"

#include "Epetra_BlockMap.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"

#include "BDF2_fnBase.hpp"
#include "BDF2_Dae.hpp"


// 2D array index algebra for 3*n arrays
#define IND(i,j) ((i)+3*(j))    

class nodal1Dfem : public BDF2::fnBase 
{
public:
  nodal1Dfem(int nnode, double x0, double x1, double u0, double u1, int pnum, double d_, 
	     double atol_, double rtol_) :
    u_left(u0), u_right(u1), prob(pnum), d(d_),rtol(rtol_), atol(atol_)
  {
    // create the Epetra map for the nodal values
    Epetra_Comm* comm = new Epetra_SerialComm();
    nodal_map = new Epetra_BlockMap(nnode, 1, 0,*comm);
    cell_map = new Epetra_BlockMap(nnode-1, 1, 0, *comm);
    cell_map_3 = new Epetra_BlockMap(nnode-1, 3, 0, *comm); // 3 elements per entry
    nodal_map_3 = new Epetra_BlockMap(nnode, 3, 0, *comm);

    // create equi-spaced mesh
    mesh = new Epetra_Vector(*nodal_map);

    n = nnode;

    double h = (x1-x0)/(n-1);
    dx = new Epetra_Vector(*cell_map);
    dx->PutScalar(h);
    
    (*mesh)[0] = x0;
    for (int i=1; i<=n-1; i++)
      {
	(*mesh)[i] = (*mesh)[i-1] + h;
      }
    (*mesh)[n-1] = x1;

    jac = new Epetra_Vector(*nodal_map_3);

  }

  void fun(double t, Epetra_Vector& u, Epetra_Vector& udot, Epetra_Vector& f) 
  {
    ASSERT(udot.MyLength() == n);
    ASSERT(u.MyLength() == n);
    ASSERT(f.MyLength() == n);

    Epetra_Vector a(*cell_map);

    
    eval_diff_coef (u, a);
    
    f[1] = ( (*dx)[0]/3.0)*udot[1] + (a[0]/(*dx)[0])*(u[1] - u_left);
    for (int j = 1; j<n-2; j++)
      {
	f[j] = f[j] + ((*dx)[j]/6.0)*(2.0*udot[j] + udot[j+1]) - (a[j]/(*dx)[j])*(u[j+1] - u[j]);
	f[j+1] =      ((*dx)[j]/6.0)*(udot[j] + 2.0*udot[j+1]) + (a[j]/(*dx)[j])*(u[j+1] - u[j]);
      }
    f[n-2] = f[n-2] + ((*dx)[n-2]/3.0)*udot[n-2] - (a[n-2]/(*dx)[n-2])*(u_right - u[n-2]);
    
    f[0] = u[0] - u_left;
    f[n-1] = u[n-1] - u_right;
    
  }
    
  
  void precon(Epetra_Vector& u, Epetra_Vector& Pu) 
  {
    ASSERT(u.MyLength() == Pu.MyLength());

    tdsolve (*jac, u, 0);
    
    Pu = u;
  }
  
  double enorm(Epetra_Vector& u, Epetra_Vector& du) 
  {

    double en = 0.0;
    for (int j=0; j<u.MyLength(); j++)
      {
	double tmp = abs(du[j])/(atol+rtol*abs(u[j]));
	en = std::max<double>(en, tmp);
      }

    return  en;
    
    
  }


  void update_precon(double t, Epetra_Vector& up, double h, int& errc) 
  {
    
    // Jacobian of the linear term in udot.
    eval_mass_matrix (*dx, *jac);

    
    jac->Scale(1.0 / h);

    // Jacobian of the linear diffusion term in u.
    Epetra_Vector a(*cell_map);
    eval_diff_coef (up, a);
    
    double tmp;
    
    for (int j = 0; j<=n-2; j++)
      { 
	tmp = a[j]/(*dx)[j];
	(*jac)[IND(1,j)]   = (*jac)[IND(1,j)] + tmp;
	(*jac)[IND(2,j)]   = (*jac)[IND(2,j)] - tmp;
	(*jac)[IND(0,j+1)] = (*jac)[IND(0,j+1)] - tmp;
	(*jac)[IND(1,j+1)] = (*jac)[IND(1,j+1)] + tmp;
      }

    // Dirichlet BC at the left-most node.
    (*jac)[IND(1,0)] = 1.0;
    (*jac)[IND(2,0)] = 0.0;
    (*jac)[IND(0,1)] = 0.0;

    // Dirichlet BC at the right-most node.
    (*jac)[IND(1,n-1)] = 1.0;
    (*jac)[IND(0,n-1)] = 0.0;
    (*jac)[IND(2,n-2)] = 0.0;

    
    tdfactor (*jac, 0);

    errc = 0;
    

  }

  bool is_admissible(Epetra_Vector& u)
  {
    return true;
  }



  void eval_diff_coef (Epetra_Vector& u, Epetra_Vector& a)
  {
    ASSERT(a.MyLength() == n-1);
    
    switch (prob)
      {
      case 1: 
	a.PutScalar(d);
	break;
      case 2:
	for (int j=0; j<n-1; j++)
	  {
	    a[j] = d + std::max<double>(0.0, 0.5*(u[j]+u[j+1]));	  }
	break;
      }	
    
  }
  

  void eval_mass_matrix(Epetra_Vector& dx, Epetra_Vector& m)
  {
    ASSERT(dx.MyLength() == n-1);
    ASSERT(3*(dx.MyLength()+1) == m.MyLength());
    
    m[IND(1,0)] = 0.0;
    for (int j=0; j<dx.MyLength(); j++)
      {
	m[IND(1,j)] = m[IND(1,j)] + dx[j]/3.0;
	m[IND(2,j)] = dx[j]/6.0;
	m[IND(0,j+1)] = dx[j]/6.0;
	m[IND(1,j+1)] = dx[j]/3.0;
      } 
  }


  void tdfactor(Epetra_Vector& a, int is)
  {

    for (int j=1+is; j < n-is; j++)
      {
	a[IND(2,j-1)] = a[IND(2,j-1)]/a[IND(1,j-1)];
	a[IND(1,j)]   = a[IND(1,j)] - a[IND(0,j)]*a[IND(2,j-1)];
      }

  }

  
  void tdsolve(Epetra_Vector& a, Epetra_Vector& x, int is)
  {
    ASSERT(a.MyLength() == 3*(x.MyLength()));

    // Forward substitution.
    x[is] = x[is]/a[IND(1,is)];
    for (int j = 1+is; j<=n-1-is; j++)
      {
	x[j] = (x[j] - a[IND(0,j)]*x[j-1])/a[IND(1,j)];
      }

    // Backward substitution.
    for (int j = n-2-is; j>=is; j--)
      {
	x[j] = x[j] - a[IND(2,j)]*x[j+1];
      }
  }

  
  void compute_udot(double t, Epetra_Vector& u, Epetra_Vector& udot)
  {
    Epetra_Vector a(*cell_map);
    Epetra_Vector m(*nodal_map_3);

    ASSERT(u.MyLength() == n);
    ASSERT(udot.MyLength() == n);
    
    eval_diff_coef(u, a);
    udot[1] = - (a[0]/(*dx)[0])*(u[1] - u_left);
    for (int j = 1; j< n-2; j++)
      {
	udot[j]   = udot[j] + (a[j]/(*dx)[j])*(u[j+1] - u[j]);
	udot[j+1] =         - (a[j]/(*dx)[j])*(u[j+1] - u[j]);
      }
    udot[n-2] = udot[n-2] + (a[n-2]/(*dx)[n-2])*(u_right - u[n-2]);

    // Evaluate the mass matrix.
    eval_mass_matrix(*dx, m);



    // Solve for UDOT on the interior nodes only.

    tdfactor (m, 1);

    tdsolve (m, udot, 1);

      // Time independent Dirichlet BV.
    udot[0] = 0.0;
    udot[n-1] = 0.0;   

  } 

  Epetra_BlockMap* nodal_map;
  Epetra_BlockMap* cell_map;
  Epetra_BlockMap* nodal_map_3;
  Epetra_BlockMap* cell_map_3;


  Epetra_Vector* mesh;
  Epetra_Vector* dx;
  
  double u_left, u_right;
  
  Epetra_Vector* jac;
    
  int prob; 
  int n;
  double d;
  double atol, rtol;

};




TEST(Nodal_1D_FEM) {

  // create the parameter list for BDF2
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::rcp(new Teuchos::ParameterList());
  plist->set("Nonlinear solver max iterations", 10);
  plist->set("Nonlinear solver tolerance", 0.01);
  plist->set("NKA max vectors",5);
  plist->set("NKA drop tolerance",0.05);
  Teuchos::ParameterList& verblist = plist->sublist("VerboseObject");
  
  verblist.set("Verbosity Level","none");

  // set the parameters for this problem
  int nnodes = 201;
  double x0  = 0.0;
  double x1  = 1.0;
  double u0  = 0.0;
  double u1  = 0.0;
  int problem_number = 2;
  double diff_coef = 0.0002;

  // set parameters for the error function
  double atol = 1.0e-5;
  double rtol = 0.0;
  

  // create the PDE problem
  nodal1Dfem NF (nnodes, x0, x1, u0, u1, problem_number, diff_coef, atol, rtol);
  
  // create the time stepper
  BDF2::Dae TS( NF, *NF.nodal_map);
  TS.setParameterList(plist);
  
  // create the initial condition
  Epetra_Vector u(*NF.nodal_map);
  for (int j=0; j< u.MyLength(); j++)
    u[j] = sin(4.0*atan(1.0)* (*NF.mesh)[j]);


  // initial time
  double t=0.0;

  // final time
  double tout = 0.2;

  // create udot and compute its initial value
  Epetra_Vector *udot  = new Epetra_Vector(*NF.nodal_map);
  NF.compute_udot(t, u, *udot);

  // initial time step
  double h = 1.0e-5;
  double hnext;

  // initialize the state of the time stepper
  TS.set_initial_state(t, u, *udot);

  // iterate until the final time
  int i=0;
  double tlast;
  do {
    
    TS.bdf2_step(h,0.0,20,u,hnext);

    TS.commit_solution(h,u);
    
    TS.write_bdf2_stepping_statistics();

    h = hnext;
    i++;

    tlast=TS.most_recent_time();
  } while (tout >= tlast);

  
  CHECK_EQUAL(i,147);
  CHECK_CLOSE(tlast,0.2044366451977221 ,1.0e-15);
}
