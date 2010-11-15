#include "Glenn.hpp"
#include "Beaker.hpp"

Glenn::Glenn(Beaker *b) 
{
  b_ = b;
} // end Glenn() constructor

Glenn::~Glenn() 
{
  if (b_) b_ = NULL;
} // end Glenn destructor

void Glenn::solve(Beaker::BeakerComponents* components, 
                  double final_time, double ts_size,
                  const Beaker::BeakerParameters& parameters)
{

  // speciate to get initial guess (and realistic activity coefficients)
  b_->Speciate(components, parameters);
  b_->print_results();

  double time = 0.;
  // just converting seconds to years -- both obviously zero in this case
  b_->print_results(time / 365. / 24. / 3600.);
  do {
    b_->ReactionStep(components, parameters, ts_size);
    // increment time
    time += ts_size;
    b_->print_results(time / 365. / 24. / 3600.);
  } while (time < final_time);

  b_->print_results();

} // end solve()

