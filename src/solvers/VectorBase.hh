/*
  Solvers

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)

  Example of a Vector base class used in our templates. The 
  routines below are mandatory for any implementation of a 
  vector class to be compatible with Amanzi.

  Note that this class is NEVER used, it simply is documentation of
  the implied interface.
*/

class VectorBase {
 public:
  // copy constructor -- this should change to Clone()! --etc
  VectorBase(const VectorBase& u);

  // Underlying vector map.
  VectorSpace& Map();

  // (*this) = a * (*this) + b * u
  VectorBase& Update(double b, const VectorBase& u, double a);

  // Eucleadian norm of vector, returns error code 0 if success, !0 otherwise
  int Norm2(double* norm) const;

  // dot product: a = (*this) * u, returns error code 0 if success, !0 otherwise
  int Dot(const VectorBase& u, double* a);
};
