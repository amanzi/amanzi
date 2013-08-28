/*
 This is the mimetic discretization component of the Amanzi code. 
 
 Copyright 2010-20XX held jointly by LANS/LANL, LBNL, and PNNL. 
 Amanzi is released under the three-clause BSD License. 
 The terms of use and "as is" disclaimer for this license are 
 provided in the top-level COPYRIGHT file.
 
Version: 2.0
Release name: naka-to.
Author: Konstantin Lipnikov (lipnikov@lanl.gov)

Usage: 
*/

#ifndef  __LAPACK_HH__
#define  __LAPACK_HH__

namespace Amanzi {
namespace WhetStone {

/* generic mashine */
#define PREFIX
#define F77_LAPACK_MANGLE(lcase,UCASE) lcase ## _

#define DSYEV_F77  F77_LAPACK_MANGLE(dsyev,DSYEV)
#define DGETRF_F77 F77_LAPACK_MANGLE(dgetrf,DGETRF)
#define DGETRI_F77 F77_LAPACK_MANGLE(dgetri,DGETRI)
#define DGESVD_F77 F77_LAPACK_MANGLE(dgesvd,DGESVD)
#define DPOSV_F77  F77_LAPACK_MANGLE(dposv,DPOSV)

#ifdef __cplusplus
extern "C" {
#endif

void PREFIX DSYEV_F77(const char* jobz, const char* uplo, 
                      int* n, double* a, int* lda, 
                      double* w, double* work, int* lwork, int* info);

void PREFIX DGETRF_F77(int* nrow, int* ncol, double* a, int* lda, 
                       int* ipiv, int* info); 

void PREFIX DGETRI_F77(int* n, double* a, int* lda, 
                       int* ipiv, double* work, int* lwork, int* info); 

void PREFIX DGESVD_F77(const char* jobu, const char* jobvt, 
                       int* nrow, int* ncol, double* a, int* lda, 
                       double* s, double *u, int* ldu, double* vt, int* ldvt, 
                       double* work, int* lwork, int* info);

void PREFIX DPOSV_F77(const char* uplo, 
                      int* n, int* nrhs, double *a, int* lda, double *b, int* ldb, 
                      int* info);

#ifdef __cplusplus
}
#endif

}  // namespace WhetStone
}  // namespace Amanzi

#endif
