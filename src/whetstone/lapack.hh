/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.
 
  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  LAPACK wrappers.
*/

#ifndef  AMANZI_LAPACK_HH_
#define  AMANZI_LAPACK_HH_

namespace Amanzi {
namespace WhetStone {

/* generic mashine */
#define PREFIX
#define F77_LAPACK_MANGLE(lcase,UCASE) lcase ## _

#define DSYEV_F77  F77_LAPACK_MANGLE(dsyev,DSYEV)
#define DGETRF_F77 F77_LAPACK_MANGLE(dgetrf,DGETRF)
#define DGETRS_F77 F77_LAPACK_MANGLE(dgetrs,DGETRS)
#define DGETRI_F77 F77_LAPACK_MANGLE(dgetri,DGETRI)
#define DGESVD_F77 F77_LAPACK_MANGLE(dgesvd,DGESVD)
#define DPOSV_F77  F77_LAPACK_MANGLE(dposv,DPOSV)
#define DGESV_F77  F77_LAPACK_MANGLE(dgesv,DGESV)
#define DGEEV_F77  F77_LAPACK_MANGLE(dgeev,DGEEV)
#define DGELS_F77  F77_LAPACK_MANGLE(dgels,DGELS)
#define DGEES_F77  F77_LAPACK_MANGLE(dgees,DGEES)
#define DGTSV_F77  F77_LAPACK_MANGLE(dgtsv,DGTSV)
#define DPOTRF_F77 F77_LAPACK_MANGLE(dpotrf,DPOTRF)
#define DPOTRI_F77 F77_LAPACK_MANGLE(dpotri,DPOTRI)

#ifdef __cplusplus
extern "C" {
#endif

void PREFIX DSYEV_F77(const char* jobz, const char* uplo, 
                      int* n, double* a, int* lda, 
                      double* w, double* work, int* lwork, int* info);

void PREFIX DGETRF_F77(int* nrow, int* ncol, double* a, int* lda, 
                       int* ipiv, int* info); 

void PREFIX DGETRS_F77(const char* uplo,
                       int* n, int* nrhs, double* a, int* lda, 
                       int* ipiv, double* b, int* ldb, int* info); 

void PREFIX DGETRI_F77(int* n, double* a, int* lda, 
                       int* ipiv, double* work, int* lwork, int* info); 

void PREFIX DGESVD_F77(const char* jobu, const char* jobvt, 
                       int* nrow, int* ncol, double* a, int* lda, 
                       double* s, double *u, int* ldu, double* vt, int* ldvt, 
                       double* work, int* lwork, int* info);

void PREFIX DPOSV_F77(const char* uplo, 
                      int* n, int* nrhs, double *a, int* lda, double *b, int* ldb, 
                      int* info);

void PREFIX DGESV_F77(int* n, int* nrhs, double* a, int* lda, int* ipiv, 
                      double* b, int* ldb, int* info);

void PREFIX DGEEV_F77(const char* jobu, const char* jobvt,
                      int* n, double* a, int* lda, double* wr, double* wi, 
                      double* vl, int* ldvl, double* vr, int* ldvr, 
                      double* work, int* lwork, int* info);

void PREFIX DGELS_F77(const char* trans,
                      int* m, int* n, int* nrhs, double* a, int* lda,
                      double* b, int* ldb, double* work, int* lwork, int* info);
 
void PREFIX DGEES_F77(const char* jobvs, const char* sort, int* select(double, double),
                      int* n, double* a, int* lda, int* sdim, 
                      double* wr, double* wi, double* vs, int* ldvs,
                      double* work, int* lwork, int* bwork, int* info);

void PREFIX DGTSV_F77(int* n, int* nrhs,
                      double* down, double* diag, double* up,
                      double* b, int* ldb, int* info);

void PREFIX DPOTRF_F77(const char* uplo, int* n, double* a, int* lda, int* info); 

void PREFIX DPOTRI_F77(const char* uplo, int* n, double* a, int* lda, int* info); 

#ifdef __cplusplus
}
#endif

}  // namespace WhetStone
}  // namespace Amanzi

#endif
