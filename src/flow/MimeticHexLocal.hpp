#ifndef __MIMETICHEXLOCAL_H__
#define __MIMETICHEXLOCAL_H__

#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"

class MimeticHexLocal {

public:
  MimeticHexLocal() {}
  MimeticHexLocal(double x_[][3]);
  ~MimeticHexLocal() {};
  
  void update(double x_[][3]);

  void mass_matrix(Epetra_SerialDenseMatrix &matrix, bool invert = false)
      { mass_matrix(matrix, 1.0, invert); }
  void mass_matrix(Epetra_SerialDenseMatrix &matrix, double K, bool invert = false);
  void mass_matrix(Epetra_SerialDenseMatrix &matrix, const Epetra_SerialSymDenseMatrix &K, bool invert = false);
  
  void diff_op(double, const double&, const double[], double&, double[]);
  void diff_op(const Epetra_SerialSymDenseMatrix&, const double&, const double[], double&, double[]);
  void diff_op(double, const double&, const Epetra_SerialDenseVector&, double&, Epetra_SerialDenseVector&);
  void diff_op(const Epetra_SerialSymDenseMatrix&, const double&, const Epetra_SerialDenseVector&, double&, Epetra_SerialDenseVector&);

private:

  double hvol;
  double cwgt[8];
  Epetra_SerialDenseMatrix face_normal;
  
};

#endif
