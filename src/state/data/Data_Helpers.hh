/*
  State

  Copyright 2010-202x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Helpers that know how to read/write/etc data.
*/

#ifndef AMANZI_STATE_DATA_HELPERS_HH_
#define AMANZI_STATE_DATA_HELPERS_HH_

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"

#include "errors.hh"

#include "StateDefs.hh"
#include "Operator_DataHelpers.hh"

namespace Amanzi {

namespace Functions {
class BoundaryFunction;
}
struct TensorVector;
class CompositeVector;
class TreeVector;
class TreeVectorSpace;

namespace Helpers {

// ======================================================================
// Visualization
// ======================================================================

// Default simply dispatches to Vis.  This fails to compile! if not either
// specifically implemented in Visualization class or specialized below.
template <typename T>
void
WriteVis(const Visualization& vis,
         const Key& fieldname,
         const std::vector<std::string>* subfieldnames,
         const T& t)
{
  UserWriteVis(vis, fieldname, subfieldnames, t);
}


// ======================================================================
// Checkpoint
// ======================================================================

// WriteCheckpoint writes data to file
//
// Default simply dispatches to Checkpoint.  This fails to compile! if not
// either specifically implemented in Checkpoint class or specialized below.
template <typename T>
void
WriteCheckpoint(const Checkpoint& chkp,
                const Key& fieldname,
                const std::vector<std::string>* subfieldnames,
                const T& t)
{
  UserWriteCheckpoint(chkp, fieldname, subfieldnames, t);
}

// ReadCheckpoint reads data from file
//
// Default simply dispatches to Checkpoint.  This fails to compile! if not
// either specifically implemented in Checkpoint class or specialized below.
template <typename T>
bool
ReadCheckpoint(const Checkpoint& chkp,
               const Key& fieldname,
               const std::vector<std::string>* subfieldnames,
               T& t)
{
  return UserReadCheckpoint(chkp, fieldname, subfieldnames, t);
}


// ======================================================================
// Initialize
// ======================================================================

template <typename T>
bool
Initialize(Teuchos::ParameterList& plist,
           T& t,
           const Key& fieldname,
           const std::vector<std::string>* subfieldnames)
{
  return UserInitialize(plist, t, fieldname, subfieldnames); // user imlementation is required
}

// ======================================================================
// Operator Assignment
// ======================================================================
template <typename T>
typename std::enable_if<!std::is_assignable<T&, const T&>::value, void>::type
Assign(T& dest, const T& source)
{
  // please C++17... then remove the enable_if junk above and use this instead
  // if constexpr(std::is_assignable<T&, const T&>::value) {
  //   dest = source;
  // } else {
  //   UserAssign(dest, source);
  // }
  UserAssign(dest, source);
}

template <typename T>
typename std::enable_if<std::is_assignable<T&, const T&>::value, void>::type
Assign(T& dest, const T& source)
{
  dest = source;
}


// ======================================================================
// Equivalency of factories
//
// Only needed for factories that cannot be default-constructed.
// ======================================================================
// another that is much easier in C++17, but until then we spell it out
template <typename F>
bool
Equivalent(const F& one, const F& two)
{
  return one == two;
}

template <>
bool
Equivalent(const Epetra_Map& one, const Epetra_Map& two);
template <>
bool
Equivalent(const Epetra_BlockMap& one, const Epetra_BlockMap& two);
template <>
bool
Equivalent(const TreeVectorSpace& one, const TreeVectorSpace& two);


// ======================================================================
// Specializations for simple data types
// ======================================================================
template <>
void
WriteVis<double>(const Visualization& vis,
                 const Key& fieldname,
                 const std::vector<std::string>* subfieldnames,
                 const double& t);

template <>
void
WriteCheckpoint<double>(const Checkpoint& chkp,
                        const Key& fieldname,
                        const std::vector<std::string>* subfieldnames,
                        const double& t);


template <>
bool
ReadCheckpoint<double>(const Checkpoint& chkp,
                       const Key& fieldname,
                       const std::vector<std::string>* subfieldnames,
                       double& t);

template <>
bool
Initialize<double>(Teuchos::ParameterList& plist,
                   double& t,
                   const Key& fieldname,
                   const std::vector<std::string>* subfieldnames);


template <>
void
WriteVis<int>(const Visualization& vis,
              const Key& fieldname,
              const std::vector<std::string>* subfieldnames,
              const int& t);

template <>
void
WriteCheckpoint<int>(const Checkpoint& chkp,
                     const Key& fieldname,
                     const std::vector<std::string>* subfieldnames,
                     const int& t);

template <>
bool
ReadCheckpoint<int>(const Checkpoint& chkp,
                    const Key& fieldname,
                    const std::vector<std::string>* subfieldnames,
                    int& t);

template <>
bool
Initialize<int>(Teuchos::ParameterList& plist,
                int& t,
                const Key& fieldname,
                const std::vector<std::string>* subfieldnames);


// ======================================================================
// Specializations for CompositeVector
// ======================================================================
template <>
void
WriteVis<CompositeVector>(const Visualization& vis,
                          const Key& fieldname,
                          const std::vector<std::string>* subfieldnames,
                          const CompositeVector& vec);

template <>
void
WriteCheckpoint<CompositeVector>(const Checkpoint& chkp,
                                 const Key& fieldname,
                                 const std::vector<std::string>* subfieldnames,
                                 const CompositeVector& vec);

template <>
bool
ReadCheckpoint<CompositeVector>(const Checkpoint& chkp,
                                const Key& fieldname,
                                const std::vector<std::string>* subfieldnames,
                                CompositeVector& vec);

template <>
bool
Initialize<CompositeVector>(Teuchos::ParameterList& plist,
                            CompositeVector& t,
                            const Key& fieldname,
                            const std::vector<std::string>* subfieldnames);


// ======================================================================
// Specializations for Epetra_Vector
// ======================================================================
template <>
void
WriteVis<Epetra_Vector>(const Visualization& vis,
                        const Key& fieldname,
                        const std::vector<std::string>* subfieldnames,
                        const Epetra_Vector& vec);

template <>
void
WriteCheckpoint<Epetra_Vector>(const Checkpoint& chkp,
                               const Key& fieldname,
                               const std::vector<std::string>* subfieldnames,
                               const Epetra_Vector& vec);

template <>
bool
ReadCheckpoint<Epetra_Vector>(const Checkpoint& chkp,
                              const Key& fieldname,
                              const std::vector<std::string>* subfieldnames,
                              Epetra_Vector& vec);

template <>
bool
Initialize<Epetra_Vector>(Teuchos::ParameterList& plist,
                          Epetra_Vector& t,
                          const Key& fieldname,
                          const std::vector<std::string>* subfieldnames);

// Specializations for TreeVector
// ======================================================================
template <>
void
WriteVis<TreeVector>(const Visualization& vis,
                     const Key& fieldname,
                     const std::vector<std::string>* subfieldnames,
                     const TreeVector& vec);

template <>
void
WriteCheckpoint<TreeVector>(const Checkpoint& chkp,
                            const Key& fieldname,
                            const std::vector<std::string>* subfieldnames,
                            const TreeVector& vec);

template <>
bool
ReadCheckpoint<TreeVector>(const Checkpoint& chkp,
                           const Key& fieldname,
                           const std::vector<std::string>* subfieldnames,
                           TreeVector& vec);

template <>
bool
Initialize<TreeVector>(Teuchos::ParameterList& plist,
                       TreeVector& t,
                       const Key& fieldname,
                       const std::vector<std::string>* subfieldnames);


// ======================================================================
// Specializations for geometric objects
// ======================================================================
template <>
void
WriteVis<AmanziGeometry::Point>(const Visualization& vis,
                                const Key& fieldname,
                                const std::vector<std::string>* subfieldnames,
                                const AmanziGeometry::Point& vec);

template <>
inline void
WriteCheckpoint<AmanziGeometry::Point>(const Checkpoint& chkp,
                                       const Key& fieldname,
                                       const std::vector<std::string>* subfieldnames,
                                       const AmanziGeometry::Point& p)
{}

template <>
inline bool
ReadCheckpoint<AmanziGeometry::Point>(const Checkpoint& chkp,
                                      const Key& fieldname,
                                      const std::vector<std::string>* subfieldnames,
                                      AmanziGeometry::Point& p)
{
  return true;
}

template <>
bool
Initialize<AmanziGeometry::Point>(Teuchos::ParameterList& plist,
                                  AmanziGeometry::Point& p,
                                  const Key& fieldname,
                                  const std::vector<std::string>* subfieldnames);


// ======================================================================
// Specializations for WhetStone::Tensor
// ======================================================================
template <>
inline void
WriteVis<TensorVector>(const Visualization& vis,
                       const Key& fieldname,
                       const std::vector<std::string>* subfieldnames,
                       const TensorVector& tensor)
{}

template <>
inline void
WriteCheckpoint<TensorVector>(const Checkpoint& chkp,
                              const Key& fieldname,
                              const std::vector<std::string>* subfieldnames,
                              const TensorVector& tensor)
{}

template <>
inline bool
ReadCheckpoint<TensorVector>(const Checkpoint& chkp,
                             const Key& fieldname,
                             const std::vector<std::string>* subfieldnames,
                             TensorVector& tensor)
{
  return true;
}

template <>
inline bool
Initialize<TensorVector>(Teuchos::ParameterList& plist,
                         TensorVector& tensor,
                         const Key& fieldname,
                         const std::vector<std::string>* subfieldnames)
{
  return true;
}


// ======================================================================
// Specializations for Functions::BoundaryFunction
// ======================================================================
template <>
inline void
WriteVis<Functions::BoundaryFunction>(const Visualization& vis,
                                      const Key& fieldname,
                                      const std::vector<std::string>* subfieldnames,
                                      const Functions::BoundaryFunction& bc)
{}

template <>
inline void
WriteCheckpoint<Functions::BoundaryFunction>(const Checkpoint& chkp,
                                             const Key& fieldname,
                                             const std::vector<std::string>* subfieldnames,
                                             const Functions::BoundaryFunction& bc)
{}

template <>
inline bool
ReadCheckpoint<Functions::BoundaryFunction>(const Checkpoint& chkp,
                                            const Key& fieldname,
                                            const std::vector<std::string>* subfieldnames,
                                            Functions::BoundaryFunction& bc)
{
  return true;
}

template <>
inline bool
Initialize<Functions::BoundaryFunction>(Teuchos::ParameterList& plist,
                                        Functions::BoundaryFunction& bc,
                                        const Key& fieldname,
                                        const std::vector<std::string>* subfieldnames)
{
  return true;
}

template <>
inline void
Assign<Functions::BoundaryFunction>(Functions::BoundaryFunction& dest,
                                    const Functions::BoundaryFunction& source)
{
  Errors::Message msg("Functions::BoundaryFunction: assignment operator not supported.");
  Exceptions::amanzi_throw(msg);
}


// ======================================================================
// Specializations for Teuchos::Array<double>
// ======================================================================
template <>
void
WriteVis<Teuchos::Array<double>>(const Visualization& vis,
                                 const Key& fieldname,
                                 const std::vector<std::string>* subfieldnames,
                                 const Teuchos::Array<double>& vec);

template <>
void
WriteCheckpoint<Teuchos::Array<double>>(const Checkpoint& chkp,
                                        const Key& fieldname,
                                        const std::vector<std::string>* subfieldnames,
                                        const Teuchos::Array<double>& vec);

template <>
bool
ReadCheckpoint<Teuchos::Array<double>>(const Checkpoint& chkp,
                                       const Key& fieldname,
                                       const std::vector<std::string>* subfieldnames,
                                       Teuchos::Array<double>& vec);

template <>
bool
Initialize<Teuchos::Array<double>>(Teuchos::ParameterList& plist,
                                   Teuchos::Array<double>& t,
                                   const Key& fieldname,
                                   const std::vector<std::string>* subfieldnames);


} // namespace Helpers
} // namespace Amanzi

#endif
