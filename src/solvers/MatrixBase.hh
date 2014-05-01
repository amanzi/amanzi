/*
  Example of a Matrix base class used in our templates. The 
  routines below are mandatory for any implementation of a 
  Matrix class to be compatible with Amanzi.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

class MatrixBase {
 public:
  // Space for the domain of the operator.
  const VectorSpace& DomainMap() const;

  // Space for the domain of the operator.
  const VectorSpace& RangeMap() const;

  // Apply matrix, b <-- Ax, returns ierr = 0 if success, !0 otherwise
  int Apply(const Epetra_Vector& x,
	    Epetra_Vector& b) const;

  // Apply the inverse, x <-- A^-1 b, returns ierr = 0 if success, !0 otherwise
  int ApplyInverse(const Epetra_Vector& b,
		   Epetra_Vector& x) const;
};
