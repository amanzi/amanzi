#include "SolutionHistory.hh"
#include "dbc.hh"
#include "errors.hh"

namespace BDF2 {

SolutionHistory::SolutionHistory(int mvec, const Epetra_BlockMap& map)
{
  ASSERT(mvec>0);
  initialize(mvec, map);
}


SolutionHistory::SolutionHistory(int mvec, double t, const Epetra_Vector& x)
{
  ASSERT(mvec>0);
  initialize(mvec, x.Map());
  record_solution(t,x);
}


SolutionHistory::SolutionHistory(int mvec, double t, const Epetra_Vector& x, const Epetra_Vector& xdot)
{
  ASSERT(mvec>0);
  // ASSERT(x.Map().SameAs(xdot.Map()));  // x and xdot have the same structure

  initialize(mvec, x.Map() );
  record_solution(t,x,xdot);
}


void SolutionHistory::initialize(int mvec, const Epetra_BlockMap& map)
{
  ASSERT(mvec>0);

  nvec = 0;
  d.resize(mvec);
  times.resize(mvec);

  for (int j=0; j<mvec; j++)
    d[j] = Teuchos::rcp(new Epetra_Vector(map));
}


void SolutionHistory::flush_history(double t, const Epetra_Vector& x)
{
  nvec = 0;
  record_solution(t,x);
}


void SolutionHistory::flush_history(double t, const Epetra_Vector& x, const Epetra_Vector& xdot)
{
  nvec = 0;
  record_solution(t,x,xdot);
}


void SolutionHistory::record_solution(double t, const Epetra_Vector& x)
{
  // ASSERT(x.Map().SameAs( (*d[0]).Map() ));

  // update the number of vectors
  nvec++;
  if (nvec > d.size()) nvec = d.size();

  // shift the times and vectors,
  // while storing the pointer to the last one
  Teuchos::RCP<Epetra_Vector> tmp = d[nvec-1];
  for (int j=nvec-1; j>=1; j--) {
    times[j] = times[j-1];
    d[j] = d[j-1];
  }

  // insert the new vector
  times[0] = t;
  d[0] = tmp;
  *d[0] = x;

  // update the divided differences
  for (int j=1; j<=nvec-1; j++) {
    if (times[0] - times[j] == 0.0) {
      Errors::Message message("SolutionHistory: Time step was too small.");
      Exceptions::amanzi_throw(message);
    }
    double div = 1.0/(times[0] - times[j]);
    d[j]->Update(div,*d[j-1],-div);
  }
}


void SolutionHistory::record_solution(double t, const Epetra_Vector& x, const Epetra_Vector& xdot)
{
  //ASSERT(x.Map().SameAs( xdot.Map() ) );
  //ASSERT(x.Map().SameAs( (*d[0]).Map() ) );
  //ASSERT(xdot.Map().SameAs( (*d[0]).Map() ) );

  record_solution(t,x);

  // update the number of vectors
  nvec++;
  if (nvec > d.size()) nvec = d.size();

  if (d.size()>1) {
  // shift the divided differences, except the first; the new vector and
  // time index are the same as the most recent.
  Teuchos::RCP<Epetra_Vector> tmp = d[nvec-1];
  for (int j=nvec-1; j>=2; j--) {
    times[j] = times[j-1];
    d[j] = d[j-1];
  }

  // the first divided difference (same time index) is the specified derivative.
  times[1] = times[0];
  d[1] = tmp;
  *d[1] = xdot;

  // update the rest of the divided differences
  for (int j=2; j<=nvec-1; j++) {
    double div = 1.0/(times[0]-times[j]);
    d[j]->Update(div,*d[j-1],-div);
  }
  }
}


void SolutionHistory::interpolate_solution(double t, Epetra_Vector& x)
{
  //ASSERT(x.Map().SameAs( (*d[0]).Map() ) );

  int order = nvec-1;
  interpolate_solution(t, x, order);
}


void SolutionHistory::interpolate_solution(double t, Epetra_Vector& x, int order)
{
  //ASSERT(x.Map().SameAs( (*d[0]).Map() ) );
  ASSERT(order<nvec);
  ASSERT(order>=0);

  x = *d[order];
  for (int k=order-1; k>=0; k--) {
    x.Update(1.0,*d[k],t-times[k]);
  }
}


void SolutionHistory::most_recent_solution(Epetra_Vector& x)
{
  //ASSERT(x.Map().SameAs( (*d[0]).Map() ) );
  x = *d[0];
}


double SolutionHistory::most_recent_time()
{
  return times[0];
}


void SolutionHistory::time_deltas(std::vector<double>& h)
{
  h.resize(nvec-1);

  for (int j=0; j<=nvec-2; j++)
    h[j] = times[0] - times[j+1];
}


int SolutionHistory::history_size()
{
  return nvec;
}

}  // namespace BDF2
