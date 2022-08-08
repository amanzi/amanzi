#ifndef AMANZI_SOLUTION_HISTORY_H_
#define AMANZI_SOLUTION_HISTORY_H_

// This class is based on Neil Carlson's SOLUTION_HISTORY
// module that is part of LANL's Truchas code.
// Modified for Amanzi.

#include <type_traits>

#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"
#include "Epetra_Map.h"

#include "dbc.hh"
#include "errors.hh"
#include "State.hh"

namespace Amanzi {

  // The Amanzi::SolutionHistory class stores the solution history of the time
  // stepper and provides interpolation methods that are used in the DAE
  // object.

template<class Vector>
class SolutionHistory {
 public:
  SolutionHistory(const std::string& name, int mvec, double t, const Vector& x, const Teuchos::RCP<State>& S = Teuchos::null);
  SolutionHistory(const std::string& name, int mvec, double t, const Vector& x, const Vector& xdot,
                  const Teuchos::RCP<State>& S = Teuchos::null);

  // Flushes the accumulated solution vectors from an existing history
  // structure, and records the solution vector X with time index T as
  // the initial solution vector of a new history.  If XDOT is specified
  // it is also recorded as the solution vector time derivative at time
  // index T.
  void FlushHistory(double t, const Vector& x);
  void FlushHistory(double t, const Vector& x, const Vector& xdot);

  // Records the vector X with time index T as the most recent solution
  // vector in the history structure.  If the vector XDOT is present,
  // it is recorded as the solution vector time derivative at the same time
  // index.  The oldest solution vector (or the two oldest in the case XDOT
  // is present) is discarded once the history is fully populated with MVEC
  // vectors.  Note that when only one of a X/XDOT pair of vectors is
  // discarded, it is the derivative vector that gets discarded.
  void RecordSolution(double t, const Vector& x);
  void RecordSolution(double t, const Vector& x, const Vector& xdot);

  // Computes the interpolated (or extrapolated) vector X at time T using
  // polynomial interpolation from the set of solution vectors maintained
  // by the history.  ORDER, if present, specifies the interpolation
  // order using the ORDER + 1 most recent solution vectors; 1 for linear
  // interpolation, 2 for quadratic, etc.  It is an error to request an
  // order for which there is insufficient data.  If not specified, the
  // maximal interpolation order is used given the available data; once
  // the history is fully populated, the order of interpolation is MVEC-1.
  void InterpolateSolution(double t, Vector& x);
  void InterpolateSolution(double t, Vector& x, unsigned int order);

  // Function returns the most recent solution vector
  // maintained by the history.
  void MostRecentSolution(Vector& x);

  // Function returns the the time index T associated with the most
  // recent solution vector maintained by the history THIS.
  double MostRecentTime();

  // Function returns an array H of time index differences.  The first
  // element of H is the difference between the most recent time and the
  // penultimate time.  The second element is the difference between the
  // most recent time and the antepenultimate time, and so forth.  The
  // length of the result equals one less than the number of solution
  // vectors being maintained by the history.
  void TimeDeltas(std::vector<double>& h);

  // Returns the number of solution vectors currently
  // maintained in the history structure THIS.  The number will be
  // between 0 and the value of MVEC used to create the structure.
  int history_size() { return nvec_; }

  // copies pointers from our d_ array to the State
  void MoveToState();

 protected:
  void Initialize_(int mvec, const Vector& initvec);

 protected:
  unsigned int nvec_;
  std::vector<double> times_;
  std::vector<Teuchos::RCP<Vector> > d_; // divided differences
  Teuchos::RCP<State> S_;
  std::string name_;
};


/* ******************************************************************
* Constructors
****************************************************************** */
template<class Vector>
SolutionHistory<Vector>::SolutionHistory(const std::string& name, int mvec, double t, const Vector& x,
        const Teuchos::RCP<State>& S)
  : S_(S),
    nvec_(0),
    name_(name)
{
  Initialize_(mvec, x);
  RecordSolution(t, x);
}


template<class Vector>
SolutionHistory<Vector>::SolutionHistory(const std::string& name, int mvec, double t, const Vector& x, const Vector& xdot,
        const Teuchos::RCP<State>& S)
  : S_(S),
    nvec_(0),
    name_(name)
{
  Initialize_(mvec, x);
  RecordSolution(t, x, xdot);
}


/* ******************************************************************
* Initialize and allocate memory
****************************************************************** */
template<class Vector>
void SolutionHistory<Vector>::Initialize_(int mvec, const Vector& initvec) {
  nvec_ = 0;
  d_.resize(mvec);
  times_.resize(mvec);

  for (int j=0; j<mvec; j++)
    d_[j] = Teuchos::rcp(new Vector(initvec));

  // require data in state
  if (S_ != Teuchos::null) {
    for (int j=0; j<mvec; j++) {
      S_->Require<Vector>(name_+"_"+std::to_string(j), Tags::DEFAULT, name_);
    }
    MoveToState();
  }
}


/* ******************************************************************
* Modify history
****************************************************************** */
template<class Vector>
void SolutionHistory<Vector>::FlushHistory(double t, const Vector& x) {
  nvec_ = 0;
  RecordSolution(t, x);
}


template<class Vector>
void SolutionHistory<Vector>::FlushHistory(double t, const Vector& x, const Vector& xdot) {
  nvec_ = 0;
  RecordSolution(t, x, xdot);
}


/* *****************************************************************
* Update history.
****************************************************************** */
template<class Vector>
void SolutionHistory<Vector>::RecordSolution(double t, const Vector& x) {
  // update the number of vectors
  nvec_++;
  if (nvec_ > d_.size()) nvec_ = d_.size();

  // shift the times and history vectors,
  // while storing the pointer to the last one
  Teuchos::RCP<Vector> tmp = d_[nvec_-1];
  for (int j = nvec_ - 1; j >= 1; j--) {
    times_[j] = times_[j - 1];
    d_[j] = d_[j - 1];
  }

  // insert the new vector
  times_[0] = t;
  d_[0] = tmp;
  *d_[0] = x;

  // update the divided differences
  for (unsigned int j = 1; j <= nvec_ - 1; j++) {
    if (times_[0] - times_[j] == 0.0) {
      Errors::Message message("SolutionHistory: Time step is too small.");
      Exceptions::amanzi_throw(message);
    }
    double div = 1.0 / (times_[0] - times_[j]);
    d_[j]->Update(div, *d_[j - 1], -div);
  }

  // update pointers in state
  if (S_ != Teuchos::null) MoveToState();
}


/* *****************************************************************
* Update history, including differences.
****************************************************************** */
template<class Vector>
void SolutionHistory<Vector>::RecordSolution(double t, const Vector& x, const Vector& xdot)
{
  RecordSolution(t, x);

  // update the number of vectors
  nvec_++;
  if (nvec_ > d_.size()) nvec_ = d_.size();

  if (d_.size()>1) {
    // shift the divided differences, except the first; the new vector and
    // time index are the same as the most recent.
    Teuchos::RCP<Vector> tmp = d_[nvec_-1];
    for (unsigned int j = nvec_ - 1; j >= 2; j--) {
      times_[j] = times_[j - 1];
      d_[j] = d_[j - 1];
    }

    // the first divided difference (same time index) is the specified derivative.
    times_[1] = times_[0];
    d_[1] = tmp;
    *d_[1] = xdot;

    // update the rest of the divided differences
    for (unsigned int j = 2; j <= nvec_-1; j++) {
      double div = 1.0 / (times_[0] - times_[j]);
      d_[j]->Update(div, *d_[j-1], -div);
    }
  }

  // update pointers in state
  if (S_ != Teuchos::null) MoveToState();
}


/* *****************************************************************
* Interpolation routines.
****************************************************************** */
template<class Vector>
void SolutionHistory<Vector>::InterpolateSolution(double t, Vector& x)
{
  unsigned int order = nvec_ - 1;
  InterpolateSolution(t, x, order);
}


template<class Vector>
void SolutionHistory<Vector>::InterpolateSolution(double t, Vector& x, unsigned int order)
{
  AMANZI_ASSERT(order < nvec_);
  AMANZI_ASSERT(order >= 0);

  x = *d_[order];
  for (int k = order - 1; k >= 0; k--) {
    x.Update(1.0, *d_[k], t - times_[k]);
  }
}


/* *****************************************************************
* Access members.
****************************************************************** */
template<class Vector>
void SolutionHistory<Vector>::MostRecentSolution(Vector& x)
{
  x = *d_[0];
}


template<class Vector>
double SolutionHistory<Vector>::MostRecentTime()
{
  return times_[0];
}


template<class Vector>
void SolutionHistory<Vector>::TimeDeltas(std::vector<double>& h)
{
  h.resize(nvec_ - 1);

  for (unsigned int j = 0; j <= nvec_ - 2; j++) {
    h[j] = times_[0] - times_[j + 1];
  }
}


/* ******************************************************************
* Save our history in state to allow checkpoint/restart
****************************************************************** */
template<class Vector>
void SolutionHistory<Vector>::MoveToState()
{
  if (S_ != Teuchos::null) {
    for (int j = 0; j != d_.size(); ++j) {
      std::string varname = name_ + "_" + std::to_string(j);
      S_->SetPtr<Vector>(varname, Tags::DEFAULT, name_, d_[j]);
    }
  }
}



}  // namespace Amanzi

#endif
