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

#include "errors.hh"

#include "CompositeVector.hh"
#include "Op.hh"
#include "TensorVector.hh"
#include "Operator.hh"
#include "BCs.hh"
#include "BoundaryFunction.hh"

#include "FunctionFactory.hh"

#include "Checkpoint.hh"
#include "StateDefs.hh"
#include "Visualization.hh"

namespace Amanzi {
namespace Helpers {

// ======================================================================
// Visualization
// ======================================================================

// Default simply dispatches to Vis.  This fails to compile! if not either
// specifically implemented in Visualization class or specialized below.
template <typename T>
void WriteVis(const Visualization& vis, const Key& fieldname,
              const std::vector<std::string>& subfieldnames, const T& t) {
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
void WriteCheckpoint(const Checkpoint& chkp, const Key& fieldname, const T& t) {
  UserWriteCheckpoint(chkp, fieldname, t);
}

// ReadCheckpoint reads data from file
//
// Default simply dispatches to Checkpoint.  This fails to compile! if not
// either specifically implemented in Checkpoint class or specialized below.
template <typename T>
void ReadCheckpoint(const Checkpoint& chkp, const Key& fieldname, T& t) {
  UserReadCheckpoint(chkp, fieldname, t);
}


// ======================================================================
// Initialize
// ======================================================================

template <typename T>
bool Initialize(Teuchos::ParameterList& plist, T& t, const Key& fieldname,
                const std::vector<std::string>& subfieldnames) {
  return UserInitialize(plist, t, fieldname, subfieldnames);  // user imlementation is required
}

// ======================================================================
// Operator Assignment
// ======================================================================
template <typename T>
void Assign(T& dest, const T& source) {
  dest = source;
}


// ======================================================================
// Specializations for simple data types
// ======================================================================
template <>
void WriteVis<double>(const Visualization& vis, const Key& fieldname,
                      const std::vector<std::string>& subfieldnames,
                      const double& t);

template <>
void WriteCheckpoint<double>(const Checkpoint& chkp, const Key& fieldname,
                             const double& t);

template <>
void ReadCheckpoint<double>(const Checkpoint& chkp, const Key& fieldname,
                            double& t);

template <>
bool Initialize<double>(Teuchos::ParameterList& plist,
                        double& t, const Key& fieldname,
                        const std::vector<std::string>& subfieldnames);

template <>
void WriteVis<int>(const Visualization& vis, const Key& fieldname,
                   const std::vector<std::string>& subfieldnames, const int& t);

template <>
void WriteCheckpoint<int>(const Checkpoint& chkp, const Key& fieldname,
                          const int& t);

template <>
void ReadCheckpoint<int>(const Checkpoint &chkp, const Key &fieldname,
                         int &t);

template <>
bool Initialize<int>(Teuchos::ParameterList& plist,
                     int& t, const Key& fieldname,
                     const std::vector<std::string>& subfieldnames);


// ======================================================================
// Specializations for CompositeVector
// ======================================================================
template <>
void WriteVis<CompositeVector>(const Visualization& vis, const Key& fieldname,
                               const std::vector<std::string>& subfieldnames,
                               const CompositeVector& vec);

template <>
void WriteCheckpoint<CompositeVector>(const Checkpoint& chkp,
                                      const Key& fieldname,
                                      const CompositeVector& vec);

template <>
void ReadCheckpoint<CompositeVector>(const Checkpoint& chkp,
                                     const Key& fieldname,
                                     CompositeVector& vec);

template <>
bool Initialize<CompositeVector>(Teuchos::ParameterList& plist,
                                 CompositeVector& t, const Key& fieldname,
                                 const std::vector<std::string>& subfieldnames);


// ======================================================================
// Specializations for geometric objects
// ======================================================================
template <>
void WriteVis<AmanziGeometry::Point>(const Visualization& vis, const Key& fieldname,
                                     const std::vector<std::string>& subfieldnames,
                                     const AmanziGeometry::Point& vec);

template <> inline
void WriteCheckpoint<AmanziGeometry::Point>(const Checkpoint& chkp, const Key& fieldname,
                                            const AmanziGeometry::Point& p) {}

template <> inline
void ReadCheckpoint<AmanziGeometry::Point>(const Checkpoint& chkp, const Key& fieldname,
                                           AmanziGeometry::Point& p) {}

template <>
bool Initialize<AmanziGeometry::Point>(Teuchos::ParameterList& plist,
                                       AmanziGeometry::Point& p, const Key& fieldname,
                                       const std::vector<std::string>& subfieldnames);


// ======================================================================
// Specializations for Op
// ======================================================================
template <> inline
void WriteVis<Operators::Op>(const Visualization& vis, const Key& fieldname,
                             const std::vector<std::string>& subfieldnames,
                             const Operators::Op& vec) {}

template <> inline
void WriteCheckpoint<Operators::Op>(const Checkpoint& chkp, const Key& fieldname,
                                    const Operators::Op& vec) {}

template <> inline
void ReadCheckpoint<Operators::Op>(const Checkpoint& chkp, const Key& fieldname,
                                   Operators::Op& vec) {}

template <> inline
bool Initialize<Operators::Op>(Teuchos::ParameterList& plist,
                               Operators::Op& t, const Key& fieldname,
                               const std::vector<std::string>& subfieldnames) {
  return true;
}


// ======================================================================
// Specializations for Operator
// ======================================================================
template <> inline
void WriteVis<Operators::Operator>(const Visualization& vis, const Key& fieldname,
                                   const std::vector<std::string>& subfieldnames,
                                   const Operators::Operator& global_operator) {}

template <> inline
void WriteCheckpoint<Operators::Operator>(const Checkpoint& chkp, const Key& fieldname,
                                          const Operators::Operator& global_operator) {}

template <> inline
void ReadCheckpoint<Operators::Operator>(const Checkpoint& chkp, const Key& fieldname,
                                         Operators::Operator& global_operator) {}

template <> inline
bool Initialize<Operators::Operator>(Teuchos::ParameterList& plist,
                                     Operators::Operator& t, const Key& fieldname,
                                     const std::vector<std::string>& subfieldnames) {
  return true;
}

template <> inline
void Assign<Operators::Operator>(Operators::Operator& dest,
        const Operators::Operator& source) {
  Errors::Message msg("Operators::Operator: assignment operator not supported.");
  Exceptions::amanzi_throw(msg);
}


// ======================================================================
// Specializations for WhetStone::Tensor
// ======================================================================
template <> inline
void WriteVis<TensorVector>(const Visualization& vis, const Key& fieldname,
                            const std::vector<std::string>& subfieldnames,
                            const TensorVector& tensor) {}

template <> inline
void WriteCheckpoint<TensorVector>(const Checkpoint& chkp, const Key& fieldname,
                                   const TensorVector& tensor) {}

template <> inline
void ReadCheckpoint<TensorVector>(const Checkpoint& chkp, const Key& fieldname,
                                  TensorVector& tensor) {}

template <> inline
bool Initialize<TensorVector>(Teuchos::ParameterList& plist,
                              TensorVector& tensor, const Key& fieldname,
                              const std::vector<std::string>& subfieldnames) {
  return true;
}


// ======================================================================
// Specializations for Operators::BCs
// ======================================================================
template <> inline
void
WriteVis<Operators::BCs>(const Visualization& vis, const Key& fieldname,
                         const std::vector<std::string>& subfieldnames,
                         const Operators::BCs& bc) {}

template <> inline
void WriteCheckpoint<Operators::BCs>(const Checkpoint& chkp, const Key& fieldname,
                                            const Operators::BCs& bc) {}

template <> inline
void ReadCheckpoint<Operators::BCs>(const Checkpoint& chkp, const Key& fieldname,
                                           Operators::BCs& bc) {}

template <> inline
bool
Initialize<Operators::BCs>(Teuchos::ParameterList& plist,
                           Operators::BCs& bc, const Key& fieldname,
                           const std::vector<std::string>& subfieldnames) {
  return true;
}


// ======================================================================
// Specializations for Functions::BoundaryFunction
// ======================================================================
template <> inline
void
WriteVis<Functions::BoundaryFunction>(const Visualization& vis, const Key& fieldname,
                                      const std::vector<std::string>& subfieldnames,
                                      const Functions::BoundaryFunction& bc) {}

template <> inline
void
WriteCheckpoint<Functions::BoundaryFunction>(const Checkpoint& chkp,
        const Key& fieldname, const Functions::BoundaryFunction& bc) {}

template <> inline
void
ReadCheckpoint<Functions::BoundaryFunction>(const Checkpoint& chkp,
        const Key& fieldname, Functions::BoundaryFunction& bc) {}

template <> inline
bool
Initialize<Functions::BoundaryFunction>(Teuchos::ParameterList& plist,
        Functions::BoundaryFunction& bc, const Key& fieldname,
        const std::vector<std::string>& subfieldnames) {
  return true;
}

} // namespace Helpers
} // namespace Amanzi

#endif
