#include "chemistry/includes/direct_solver.hh"
#include "chemistry/includes/Block.hpp"

DirectSolver::DirectSolver(void)
    : Solver(),
      row_interchange_(0.),
      factored_(false) {
  pivoting_indices_.clear();
  row_scaling_.clear();
} // end DirectSolver constructor

void DirectSolver::Initialize(const int n) {
  set_system_size(n);
  right_hand_side_.resize(n);
  solution_.resize(n);
  A_ = new MatrixBlock(n);
  pivoting_indices_.resize(n);
  row_scaling_.resize(n);
} // end Init()

void DirectSolver::Solve(void) {
  LUDecomposition();
  LUBackSolve(solution_);
} // end Solve()

void DirectSolver::Solve(std::vector<double> &b) {
  LUDecomposition();
  LUBackSolve(b);
} // end Solve()

void DirectSolver::Solve(MatrixBlock *A, std::vector<double> &b) {
  A_->SetValues(A);
  LUDecomposition();
  LUBackSolve(b);
} // end Solve()

void DirectSolver::Solve(Block *A, std::vector<double> &b) {
  A_->SetValues(A);
  LUDecomposition();
  LUBackSolve(b);
} // end Solve()

void DirectSolver::LUDecomposition(void) {
  const double small_number = 1.e-20;
  int imax = 0;
  double temp, dum;

  double **a = A_->GetValues();
	row_interchange_ = 1.0;
	for (int i = 0; i < system_size(); i++) {
		double big = 0.0;
		for (int j = 0; j < system_size(); j++)
			if ((temp = fabs(a[i][j])) > big) big = temp;
    if (big == 0.0) std::cout << "Singular matrix in routine ludcmp";
		row_scaling_[i] = 1.0 / big;
	}
	for (int j = 0; j < system_size(); j++) {
		for (int i = 0; i < j; i++) {
			double sum = a[i][j];
			for (int k = 0; k < i; k++) 
        sum -= a[i][k] * a[k][j];
			a[i][j] = sum;
		}
		double big = 0.0;
		for (int i = j; i < system_size(); i++) {
			double sum =  a[i][j];
			for (int k = 0; k < j; k++) 
        sum -= a[i][k] * a[k][j];
			a[i][j] = sum;
			if ((dum = row_scaling_[i] * fabs(sum)) >= big) {
				big = dum;
				imax = i;
			}
		}
		if (j != imax) {
			for (int k = 0; k < system_size(); k++) {
				dum = a[imax][k];
				a[imax][k] = a[j][k];
				a[j][k] = dum;
			}
			row_interchange_ = -row_interchange_;
			row_scaling_[imax] = row_scaling_[j];
		}
		pivoting_indices_[j] = imax;
		if (a[j][j] == 0.0) a[j][j] = small_number;
		if (j != system_size() - 1) {
			dum = 1.0 / (a[j][j]);
			for (int i = j + 1; i < system_size(); i++) 
        a[i][j] *= dum;
		}
  }

  a = NULL;
  factored_ = true;

} // end LUDecomposition()


void DirectSolver::LUBackSolve(void) {
  // assumes that the right hand side has already been stored.
  for (int i = 0; i < system_size(); i++) {
    solution_[i] = right_hand_side_[i];
  }
  LUBackSolve(solution_);
} // end LUBackSolve()

void DirectSolver::LUBackSolve(std::vector<double> &b) {
  double **a = A_->GetValues();
  int ii = 0;
	for (int i = 0; i < system_size(); i++) {
		int ip=pivoting_indices_[i];
		double sum = b[ip];
		b[ip] = b[i];
		if (ii != 0) {
			for (int j = ii - 1; j < i; j++) 
        sum -= a[i][j] * b[j];
    }
		else if (sum != 0.0) {
			ii = i+1;
    }
		b[i] = sum;
	}
	for (int i = system_size() - 1; i >= 0; i--) {
		double sum = b[i];
		for (int j = i+1; j < system_size(); j++) 
      sum -= a[i][j] * b[j];
		b[i] = sum / a[i][i];
	}
  a = NULL;
} // end LUBackSolve()

DirectSolver::~DirectSolver() {
} // end DirectSolver destructor
