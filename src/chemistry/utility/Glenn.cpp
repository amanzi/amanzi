#include "Glenn.hpp"

Glenn::Glenn(Beaker *b) 
{
  b_ = b;
} // end Glenn() constructor

Glenn::~Glenn() 
{
  if (b_) b_ = NULL;
} // end Glenn destructor

void Glenn::solve(std::vector<double> &total, double final_time, double ts_size,
                  double porosity, double saturation, double water_density, 
                  double volume) 
{

  b_->initializeMolalities(1.e-9);

  double time = 0.;
  // just converting seconds to years -- both obviously zero in this case
  b_->print_results(time / 365. / 24. / 3600.);
  do {
    b_->react(total,porosity,saturation,water_density,volume,ts_size);
    // increment time
    time += ts_size;
    b_->print_results(time / 365. / 24. / 3600.);
  } while (time < final_time);

} // end solve()

