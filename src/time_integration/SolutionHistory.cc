/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#include "SolutionHistory.hh"
#include "dbc.hh"

namespace Amanzi {


SolutionHistory::SolutionHistory(const int mvec, const Teuchos::RCP<const TreeVector> x) {
  ASSERT(mvec>0);

  initialize(mvec, x);
}


SolutionHistory::SolutionHistory(const int mvec, const double t, const Teuchos::RCP<const TreeVector> x) {
  ASSERT(mvec>0);

  initialize(mvec, x);
  record_solution(t,x);
}


SolutionHistory::SolutionHistory(const int mvec, const double t, 
                                 const Teuchos::RCP<const TreeVector> x, 
                                 const Teuchos::RCP<const TreeVector> xdot) {
  ASSERT(mvec>0);

  initialize(mvec, x);
  record_solution(t,x,xdot);
}


void SolutionHistory::initialize(const int mvec, const Teuchos::RCP<const TreeVector> v) {
  ASSERT(mvec>0);

  nvec = 0;
  d.resize(mvec);
  times.resize(mvec);

  for (int j=0; j<mvec; j++)
    d[j] = Teuchos::rcp(new TreeVector(std::string("solhist"),*v));
}


void SolutionHistory::flush_history(const double t, const Teuchos::RCP<const TreeVector> x) {
  nvec = 0;
  record_solution(t,x);
}


void SolutionHistory::flush_history(const double t, const Teuchos::RCP<const TreeVector> x, 
                                    const Teuchos::RCP<const TreeVector> xdot) {
  nvec = 0;
  record_solution(t,x,xdot);
}


void SolutionHistory::record_solution(const double t, const Teuchos::RCP<const TreeVector> x) {

  // update the number of vectors
  nvec++;
  if (nvec > d.size()) nvec = d.size();

  // shift the times and vectors,
  // while storing the pointer to the last one
  Teuchos::RCP<TreeVector> tmp = d[nvec-1];
  for (int j=nvec-1; j>=1; j--)  {
    times[j] = times[j-1];
    d[j] = d[j-1];
  }

  // insert the new vector
  times[0] = t;
  d[0] = tmp;
  *(d[0]) = *x;
  

  // update the divided differences
  for (int j=1; j<=nvec-1; j++) {
    double div = 1.0/(times[0] - times[j]);
    d[j]->Update(div,*d[j-1],-div);
  }
}

void SolutionHistory::record_solution(const double t, const Teuchos::RCP<const TreeVector> x, 
                                      const Teuchos::RCP<const TreeVector> xdot) {

  record_solution(t,x);

  // update the number of vectors
  nvec++;
  if (nvec > d.size()) nvec = d.size();

  // shift the divided differences, except the first; the new vector and
  // time index are the same as the most recent.
  Teuchos::RCP<TreeVector> tmp = d[nvec-1];
  for (int j=nvec-1; j>=2; j--) {
    times[j] = times[j-1];
    d[j] = d[j-1];
  }

  // the first divided difference (same time index) is the specified derivative.
  times[1] = times[0];
  d[1] = tmp;
  *d[1] = *xdot;

  // update the rest of the divided differences
  for (int j=2; j<=nvec-1; j++) {
    double div = 1.0/(times[0]-times[j]);
    d[j]->Update(div,*d[j-1],-div);
  }


}


void SolutionHistory::interpolate_solution(const double t, const Teuchos::RCP<TreeVector> x) const {

  int order = nvec-1;

  interpolate_solution(t, x, order);
}


void SolutionHistory::interpolate_solution(const double t, const Teuchos::RCP<TreeVector> x, const int order) const {

  ASSERT(order<nvec);
  ASSERT(order>0);

  *x = *d[order];
  
  for (int k=order-1; k>=0; k--) {
    x->Update(1.0,*d[k],t-times[k]);
  }

}


void SolutionHistory::most_recent_solution(Teuchos::RCP<TreeVector> x) const {

  *x = *d[0];
}


double SolutionHistory::most_recent_time() const {
  return times[0];
}


void SolutionHistory::time_deltas(std::vector<double>& h)  const {

  h.resize(nvec-1);

  for (int j=0; j<=nvec-2; j++)
    h[j] = times[0] - times[j+1];
}


int SolutionHistory::history_size() const {
  return nvec;
}


void SolutionHistory::Print(ostream& os) const {
  
  std::vector<double>::const_iterator itimes = times.begin();
  for (std::vector<Teuchos::RCP<Amanzi::TreeVector> >::const_iterator it = d.begin();
       it != d.end(); ++it, ++itimes) {
    os << "time = " << *itimes << std::endl;
    (*it)->Print(os);
  }

  
}
}
