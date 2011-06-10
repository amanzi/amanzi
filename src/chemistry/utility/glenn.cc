/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include "glenn.hh"
#include "beaker.hh"

Glenn::Glenn(Beaker* b) {
  b_ = b;
}  // end Glenn() constructor

Glenn::~Glenn() {
  if (b_) {
    b_ = NULL;
  }
}  // end Glenn destructor

void Glenn::solve(Beaker::BeakerComponents* components,
                  double final_time, double ts_size,
                  const Beaker::BeakerParameters& parameters) {

  // speciate to get initial guess (and realistic activity coefficients)
  b_->Speciate(*components, parameters);
  b_->UpdateComponents(components);
  b_->print_results();

  std::vector<double> times;
  std::vector<double> A;
  std::vector<double> B;

  double time = 0.;

  times.push_back(time);
  A.push_back(components->total[0]);
  B.push_back(components->total[1]);

  // just converting seconds to years -- both obviously zero in this case
  b_->print_results(time / 365. / 24. / 3600.);
  do {

    b_->ReactionStep(components, parameters, ts_size);
    // increment time
    time += ts_size;
    b_->print_results(time / 365. / 24. / 3600.);

    times.push_back(time / 24. / 3600.);
    A.push_back(components->total[0]);
    B.push_back(components->total[1]);

  } while (time < final_time);


  b_->print_results();

  for (unsigned int i = 0; i < times.size(); i++) {
    std::cout << times[i] << ' ' << A[i] << ' ' << B[i] << std::endl;
  }

}  // end solve()
