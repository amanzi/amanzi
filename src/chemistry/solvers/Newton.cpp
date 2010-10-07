#include "LU.hpp"
#include "Newton.hpp"

Newton::Newton(const int n) {
  static_cast<void>(n);
} // end Newton constructor

void Newton::solve() {
  cout << "Solved!\n";
} // end solve()


Newton::~Newton() {
  if (J) delete J;
  J = NULL;
} // end Newton destructor
