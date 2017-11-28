#include "Teuchos_RCP.hpp"
#include "State.hh"
#include "PK.hh"


//
// Helper struct to store a run for the below function
// ================================================================================
struct Run {
  Run(const Teuchos::RCP<State>& S_,
      const Teuchos::RCP<PK>& pk_)
      : S(S_),
        pk(pk_) {}

  Teuchos::RCP<State> S;
  Teuchos::RCP<PK> pk;
};


//
// Helper function that mocks a coordinator to run the test.
// ================================================================================
std::pair<double,double> run_test(const Teuchos::RCP<State>& S,
        const Teuchos::RCP<PK>& pk) {

  pk->Setup();
  S->Setup();

  pk->Initialize();
  S->Initialize();

  int nsteps_good = 0;
  int nsteps_bad = 0;

  double t_final = 1.0;
  bool done = false;
  while (!done) {
    double dt = std::min(pk->get_dt(), 1.0);
    S->set_time("next", S->time()+dt);
    S->set_cycle("next", S->cycle()+1);

    bool fail = pk->AdvanceStep("", "next");

    if (fail) {
      pk->FailStep("", "next");
      nsteps_bad++;
    } else {
      pk->CommitStep("", "next");
      S->set_time("", S->time("next"));
      S->set_cycle("", S->cycle("next"));
      nsteps_good++;
    }

    done = S->time("") >= t_final-0.0001;
  }
  return std::make_pair(nsteps_good, nsteps_bad);
}



