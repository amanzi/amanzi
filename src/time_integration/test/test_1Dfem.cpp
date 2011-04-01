#include <iostream>
#include "UnitTest++.h"

#include "dbc.hh"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Epetra_MpiComm.h"
#include "Epetra_SerialComm.h"

#include "Epetra_BlockMap.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"

#include "BDF2_fnBase.hpp"
#include "BDF2_Dae.hpp"


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
    nodal_map = new Epetra_BlockMap(nnode, 1, 1,*comm);
    cell_map = new Epetra_BlockMap(nnode-1, 1, 1, *comm);
    cell_map_3 = new Epetra_BlockMap(nnode-1, 3, 1, *comm); // 3 elements per entry

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

    jac = new Epetra_Vector(*cell_map_3,n-1);
    m = new Epetra_Vector(*cell_map_3,n-1);

  }

  void fun(double t, Epetra_Vector& u, Epetra_Vector& udot, Epetra_Vector& f) 
  {
    Epetra_Vector a(*cell_map);

    
    eval_diff_coef (u, a);
    
    f[1] = ( (*dx)[0]/3.0)*udot[1] + (a[0]/(*dx)[0])*(u[1] - u_left);
    for (int j = 1; j<n-2; j++)
      {
	f[j] = f[j] + ((*dx)[j]/6.0)*(2.0*udot[j] + udot[j+1]) - (a[j]/(*dx)[j])*(u[j+1] - u[j]);
	f[j+1] =      ((*dx)[j]/6.0)*(udot[j] + 2.0*udot[j+1]) + (a[j]/(*dx)[j])*(u[j+1] - u[j]);
      }
    f[n-2] = f[n-2] + ((*dx)[n-2]/3.0)*udot[n-2] - (a[n-2]/(*dx)[n-2])*(u_right - u[n-2]);
    
    f[1] = u[1] - u_left;
    f[n-1] = u[n-1] - u_right;
    
  }
    
  
  void precon(Epetra_Vector& u, Epetra_Vector& Pu) 
  {
    // Precondition the function.
    double **jacv;
    jac->ExtractView(jacv);

    double **uv;
    u.ExtractView(uv);

    tdsolve (jacv[0], uv[0], n);
    Pu = u;
      
  }
  
  double enorm(Epetra_Vector& u, Epetra_Vector& du) 
  {
    Epetra_Vector tmp1(u);
    Epetra_Vector tmp2(u);
    Epetra_Vector tmp3(u);



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
	(*jac)[IND(1,j)] = (*jac)[IND(1,j)] + tmp;
	(*jac)[IND(2,j)] = (*jac)[IND(2,j)] - tmp;
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

    double **jacv;
    jac->ExtractView(jacv);
    
    tdfactor (jacv[0], n);

    errc = 0;
    

  }












  //private:
  
  void eval_diff_coef (Epetra_Vector& u, Epetra_Vector& a)
  {
    switch (prob)
      {
      case 1: 
	a.PutScalar(d);
	break;
      case 2:
	for (int j=0; j<n-1; j++)
	  {
	    a[j] = d + std::max<double>(0.0, 0.5*(u[j]+u[j+1]));
	  }
	break;
	
      }	
    
  }
  

  void eval_mass_matrix(Epetra_Vector& dx, Epetra_Vector& m)
  {
    ASSERT(dx.MyLength() == n-1);
    
    m[IND(1,0)] = 0.0;
    for (int j=0; j<n-1; j++)
      {
	m[IND(1,j)] = m[IND(1,j)] + dx[j]/3.0;
	m[IND(2,j)] = dx[j]/6.0;
	m[IND(0,j+1)] = dx[j]/6.0;
	m[IND(1,j+1)] = dx[j]/3.0;
      } 
  }


  void tdfactor(double* a, int len)
  {

    for (int j=1; j < len; j++)
      {
	a[IND(2,j-1)] = a[IND(2,j-1)]/a[IND(1,j-1)];
	a[IND(1,j)]   = a[IND(1,j)] - a[IND(0,j)]*a[IND(2,j-1)];
      }

  }

  
  void tdsolve(double* a, double* x, int len)
  {
    
    // Forward substitution.
    x[0] = x[0]/a[IND(1,0)];
    for (int j = 1; j<=len-1; j++)
      {
	x[j] = (x[j] - a[IND(0,j)]*x[j-1])/a[IND(1,j)];
      }
    
    // Backward substitution.
    for (int j = len-2; j>=0; j--)
      {
	x[j] = x[j] - a[IND(2,j)]*x[j+1];
      }

  }

  
  void compute_udot(double t, Epetra_Vector& u, Epetra_Vector& udot)
  {
    Epetra_Vector a(*cell_map);
    
    eval_diff_coef(u, a);
    udot[1] = - (a[0]/(*dx)[0])*(u[1] - u_left);
    for (int j = 1; j< n-2; j++)
      {
	udot[j]   = udot[j] + (a[j]/(*dx)[j])*(u[j+1] - u[j]);
	udot[j+1] =         - (a[j]/(*dx)[j])*(u[j+1] - u[j]);
      }
    udot[n-2] = udot[n-2] + (a[n-2]/(*dx)[n-2])*(u_right - u[n-2]);

    // Evaluate the mass matrix.
    eval_mass_matrix(*dx, *m);

    // Solve for UDOT on the interior nodes only.
    double **mv;
    m->ExtractView(mv);

    double **udotv;
    udot.ExtractView(udotv);
    
    tdfactor (mv[0]+1, n-1);
    tdsolve (mv[0]+1, udotv[0]+1, n-1);

      // Time independent Dirichlet BV.
    udot[0] = 0.0;
    udot[n-1] = 0.0;
   
  } 






  
  Epetra_BlockMap* nodal_map;
  Epetra_BlockMap* cell_map;
  Epetra_BlockMap* cell_map_3;

  Epetra_Vector* mesh;
  Epetra_Vector* m;
  Epetra_Vector* dx;
  
  double u_left, u_right;
  
  Epetra_Vector* jac;
    
  int prob; 
  int n;
  double d;
  double atol, rtol;

};




TEST(Nodal_1D_FEM) {

  nodal1Dfem NF (11, 0.0, 1.0, 10.0, 5.0, 2, 1.0, 0.0, 1.0e-5);
  
  Epetra_Vector u(*NF.nodal_map);
  
  for (int j=0; j< u.MyLength(); j++)
    u[j] = sin(4.0*atan(1.0)* (*NF.mesh)[j]);

  double t=0.0;
  
  BDF2::Dae TS( NF, *NF.nodal_map, 10, 0.01, 2, 0.01);
  
  Epetra_Vector udot(*NF.nodal_map);
  NF.compute_udot(t, u, udot);

  double h = 1.0e-5;
  double hnext;

  TS.set_initial_state(t, u, udot);

  //TS.bdf2_step(h,0.0000001,10,u,hnext);
  
	       

}
