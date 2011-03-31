#ifndef _BDF2_DAE_HPP_
#define _BDF2_DAE_HPP_

#include "Epetra_Vector.h"
#include "Epetra_BlockMap.h"

#include "NKA.H"

#include "BDF2_State.hpp"
#include "BDF2_fnBase.hpp"


namespace BDF2 {

  class Dae {

  public:

    Dae(fnBase& fn_, Epetra_BlockMap& map, int mitr, double ntol, int mvec, double vtol); 
    Dae(fnBase& fn_, Epetra_BlockMap& map);

    void set_initial_state(const double t, const Epetra_Vector& x, const Epetra_Vector& xdot);

    void commit_solution(const double h, const Epetra_Vector& u);
    void select_step_size(const std::vector<double>& dt, const double perr, double& h); 

    void trap_step_one(double h, Epetra_Vector& u, double& hnext, int& errc);
    void bdf2_step_gen(double h, Epetra_Vector& u, double& hnext, int& errc, bool ectrl);
    void bdf2_step_simple(double h, Epetra_Vector& u, double& hnext, int& errc, bool ectrl);

    void bdf2_step(double& h, double hmin, int mtries, Epetra_Vector& u, double& hnext); 

    void solve_bce(double t, double h, Epetra_Vector& u0, Epetra_Vector& u, int& errc);


  private:

    double rmin; 
    double rmax; 
    double margin;

    State state;
    nka* fpa;
    fnBase& fn;


    // constants
    const static double RMIN = 0.25;
    const static double RMAX = 4.0;
    const static double MARGIN = 3.0;

  };

}

#endif // _BDF2_DAE_HPP_
