#ifndef __MIMETICHEXLOCAL_H__
#define __MIMETICHEXLOCAL_H__

#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"

namespace Amanzi
{

class MimeticHexLocal {

public:
  MimeticHexLocal() {}
  MimeticHexLocal(double *x) { update(x); }
  MimeticHexLocal(double x[][3]) { update(x); }
  ~MimeticHexLocal() {};

  void update(const Epetra_SerialDenseMatrix &x);
  void update(double *x);
  void update(double x[][3]);

  void mass_matrix(Epetra_SerialDenseMatrix &matrix, bool invert = false) const
      { mass_matrix(matrix, 1.0, invert); }
  void mass_matrix(Epetra_SerialDenseMatrix &matrix, double K, bool invert = false) const;
  void mass_matrix(Epetra_SerialDenseMatrix &matrix, const Epetra_SerialSymDenseMatrix &K, bool invert = false) const;

  void diff_op(double, const double&, const double[], double&, double[]) const;
  void diff_op(const Epetra_SerialSymDenseMatrix&, const double&, const double[], double&, double[]) const;
  void diff_op(double, const double&, const Epetra_SerialDenseVector&, double&, Epetra_SerialDenseVector&) const;
  void diff_op(const Epetra_SerialSymDenseMatrix&, const double&, const Epetra_SerialDenseVector&, double&, Epetra_SerialDenseVector&) const;

  void GravityFlux(const double g[], double gflux[]) const;
  void CellFluxVector(double Fface[], double Fcell[]) const;

private:

  double hvol;
  double cwgt[8];
  Epetra_SerialDenseMatrix face_normal;

};

} // close namespace Amanzi

#endif
