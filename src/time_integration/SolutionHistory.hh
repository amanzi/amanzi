#ifndef _ATS_SOLUTION_HISTORY_H_
#define _ATS_SOLUTION_HISTORY_H_

// This class is based on Neil Carlson's SOLUTION_HISTORY
// module that is part of LANL's Truchas code.

#include "TreeVector.hh"

namespace Amanzi {

// The BDF2::SolutionHistory class stores the solution history
// of the time stepper and provides interpolation methods
// that are used in the BDF2::Dae object.

class SolutionHistory {

 public:
  SolutionHistory (const int mvec, const Teuchos::RCP<const TreeVector> x); 
  SolutionHistory (const int mvec, double t, const Teuchos::RCP<const TreeVector> x);
  SolutionHistory (const int mvec, double t, const Teuchos::RCP<const TreeVector> x, 
		   Teuchos::RCP<const TreeVector> xdot);


  // Flushes the accumulated solution vectors from an existing history
  // structure, and records the solution vector X with time index T as
  // the initial solution vector of a new history.  If XDOT is specified
  // it is also recorded as the solution vector time derivative at time
  // index T.

  void flush_history(const double t, const Teuchos::RCP<const TreeVector> x);
  void flush_history(const double t, const Teuchos::RCP<const TreeVector> x, 
		     const Teuchos::RCP<const TreeVector> xdot);



  // Records the vector X with time index T as the most recent solution
  // vector in the history structure.  If the vector XDOT is present,
  // it is recorded as the solution vector time derivative at the same time
  // index.  The oldest solution vector (or the two oldest in the case XDOT
  // is present) is discarded once the history is fully populated with MVEC
  // vectors.  Note that when only one of a X/XDOT pair of vectors is
  // discarded, it is the derivative vector that gets discarded.

  void record_solution(const double t, const Teuchos::RCP<const TreeVector> x);
  void record_solution(const double t, const Teuchos::RCP<const TreeVector> x, 
		       const Teuchos::RCP<const TreeVector> xdot);



  // Computes the interpolated (or extrapolated) vector X at time T using
  // polynomial interpolation from the set of solution vectors maintained
  // by the history.  ORDER, if present, specifies the interpolation
  // order using the ORDER + 1 most recent solution vectors; 1 for linear
  // interpolation, 2 for quadratic, etc.  It is an error to request an
  // order for which there is insufficient data.  If not specified, the
  // maximal interpolation order is used given the available data; once
  // the history is fully populated, the order of interpolation is MVEC-1.

  void interpolate_solution(const double t, const Teuchos::RCP<TreeVector> x) const;
  void interpolate_solution(const double t, const Teuchos::RCP<TreeVector> x, const int order) const;



  // Function returns the most recent solution vector
  // maintained by the history.

  void most_recent_solution(Teuchos::RCP<TreeVector> x) const;



  // Function returns the the time index T associated with the most
  // recent solution vector maintained by the history THIS.

  double most_recent_time() const;



  // Function returns an array H of time index differences.  The first
  // element of H is the difference between the most recent time and the
  // penultimate time.  The second element is the difference between the
  // most recent time and the antepenultimate time, and so forth.  The
  // length of the result equals one less than the number of solution
  // vectors being maintained by the history.

  void time_deltas(std::vector<double>& h) const;



  // Returns the number of solution vectors currently
  // maintained in the history structure THIS.  The number will be
  // between 0 and the value of MVEC used to create the structure.

  int history_size() const;

  void Print(ostream& os) const;

 private:

  void initialize(const int mvec, const Teuchos::RCP<const TreeVector> v);

  int nvec;
  std::vector<double> times;                   // times
  std::vector<Teuchos::RCP<TreeVector> > d;    // divided(?) differences

};
}

#endif // _BDF2_SOLUTION_HISTORY_H_
