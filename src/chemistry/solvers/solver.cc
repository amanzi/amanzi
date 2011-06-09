/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include "chemistry/includes/solver.hh"

Solver::Solver(void)
    : system_size_(0),
      A_(NULL) {
  right_hand_side_.clear();
  solution_.clear();
} // end Solver constructor

void Solver::Initialize(const int n) {
  set_system_size(n);
  right_hand_side_.resize(n);
  solution_.resize(n);
  if (A_) delete A_;
  A_ = new MatrixBlock(n);
} // end Init()

Solver::~Solver() {
  if (A_) delete A_;
  A_ = NULL;
} // end Solver destructor
