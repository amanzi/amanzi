/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
  State

  Helpers that know how to read/write/etc data.
*/

/*!

.. _constants-scalar-spec:
.. admonition:: constants-scalar-spec

   * `"value`" ``[double]`` Value of a scalar constant

.. _constants-dense-vector-spec:
.. admonition:: constants-dense-vector-spec

   * `"value`" ``[Array(double)]`` Value of a dense, local vector.

.. _constants-point-spec:
.. admonition:: constants-point-spec

   * `"value`" ``[Array(double)]`` Array containing the values of the point.


.. _constants-composite-vector-spec:
.. admonition:: constants-composite-vector-spec

   * `"constant`" ``[double]`` **optional** Constant value.
   * `"value`" ``[double]`` **optional** Constant value, same as `"constant`" above.
   * `"function`" ``[composite-vector-function-spec-list]`` **optional**
     Initialize from a function, see CompositeVectorFunction_
   * `"restart file`" ``[string]`` **optional** Path to a checkpoint file from
     which to read the values.
   * `"cells from file`" ``[string]`` **optional** Same as `"restart file`",
     but only reads the cell component.
   * `"exodus file initialization`" ``[exodus-file-initialization-spec]``
     **optional** See `Exodus File Initialization`_.
   * `"initialize from 1D column`" ``[column-file-initialization-spec]``
     **optional** See `Column File Initialization`_.

*/

#ifndef AMANZI_STATE_DATA_HELPERS_HH_
#define AMANZI_STATE_DATA_HELPERS_HH_

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"

#include "errors.hh"

#include "CompositeVector.hh"
#include "StateDefs.hh"
#include "Operator_DataHelpers.hh"

namespace Amanzi {

// namespace Functions {
// class BoundaryFunction;
// }
struct TensorVector;
class TreeVector;
class TreeVectorSpace;
template <typename T>
struct MultiPatch;

namespace Helpers {

// ======================================================================
// Visualization
// ======================================================================

// Default simply dispatches to Vis.  This fails to compile! if not either
// specifically implemented in Visualization class or specialized below.
template <typename T>
void
WriteVis(const Visualization& vis, Teuchos::ParameterList& attrs, const T& t)
{
  vis.write(attrs, t);
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
WriteCheckpoint(const Checkpoint& chkp, Teuchos::ParameterList& attrs, const T& t)
{
  // hack to have dof indices in length 1 vectors for checkpointing/regression
  // tests to match master's behavior
  attrs.set("always write subfield dof", true);
  chkp.write(attrs, t);
}

// ReadCheckpoint reads data from file
//
// Default simply dispatches to Checkpoint.  This fails to compile! if not
// either specifically implemented in Checkpoint class or specialized below.
template <typename T>
void
ReadCheckpoint(const Checkpoint& chkp, Teuchos::ParameterList& attrs, T& t)
{
  chkp.read(attrs, t);
}


// ======================================================================
// Initialize
// ======================================================================

template <typename T>
bool
Initialize(Teuchos::ParameterList& plist, T& t)
{
  return false;
}

// ======================================================================
// Operator Assignment
// ======================================================================
template <typename T>
void
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
  if (&one == &two) return true;
  if constexpr (std::is_member_function_pointer<decltype(&F::isSameAs)>::value) {
    return one.isSameAs(two);
  } else {
    return one == two;
  }
}


// ======================================================================
// Specializations for simple data types
// ======================================================================
template <>
bool
Initialize<double>(Teuchos::ParameterList& plist, double& t);


template <>
bool
Initialize<int>(Teuchos::ParameterList& plist, int& t);

template <>
bool
Initialize<bool>(Teuchos::ParameterList& plist, bool& t);


// ======================================================================
// Specializations for CompositeVector
// ======================================================================
template <>
bool
Initialize<CompositeVector>(Teuchos::ParameterList& plist, CompositeVector& t);


// ======================================================================
// Specializations for Vector
// ======================================================================
template <>
bool
Initialize<Vector_type>(Teuchos::ParameterList& plist, Vector_type& t);

// ======================================================================
// Specializations for TreeVector
// ======================================================================
template <>
bool
Initialize<TreeVector>(Teuchos::ParameterList& plist, TreeVector& t);

// It isn't clear that this should work in general, because TreeVectors may
// combine CompositeVectors that exist on incompatible domain comms.  It hasn't
// crashed yet though...
template <>
void
WriteVis<TreeVector>(const Visualization& vis,
                     Teuchos::ParameterList& attrs,
                     const TreeVector& vec);

// It isn't clear that this should work in general, because TreeVectors may
// combine CompositeVectors that exist on incompatible domain comms.  In
// particular I don't think this should work in multi-file checkpointing.  But
// it hasn't crashed yet...
template <>
void
WriteCheckpoint<TreeVector>(const Checkpoint& chkp,
                            Teuchos::ParameterList& attrs,
                            const TreeVector& vec);

// It isn't clear that this should work in general, because TreeVectors may
// combine CompositeVectors that exist on incompatible domain comms.  In
// particular I don't think this should work in multi-file checkpointing.  But
// it hasn't crashed yet...
template <>
void
ReadCheckpoint<TreeVector>(const Checkpoint& chkp, Teuchos::ParameterList& attrs, TreeVector& vec);

// ======================================================================
// Specializations for geometric objects
// ======================================================================

template <>
bool
Initialize<AmanziGeometry::Point>(Teuchos::ParameterList& plist, AmanziGeometry::Point& p);


// ======================================================================
// Specializations for WhetStone::Tensor
// ======================================================================
// no support for checkpoint tensors yet... they should be functions only
// anyway at this point and therefore need not be checkpointed.
template <>
inline bool
Initialize<TensorVector>(Teuchos::ParameterList& plist, TensorVector& tensor)
{
  return false;
}

// no support for vis of tensors
template <>
inline void
WriteVis<TensorVector>(const Visualization& vis,
                       Teuchos::ParameterList& attrs,
                       const TensorVector& vec)
{}

// no support for checkpoint tensors yet... they should be functions only
// anyway at this point and therefore need not be checkpointed.
template <>
inline void
WriteCheckpoint<TensorVector>(const Checkpoint& chkp,
                              Teuchos::ParameterList& attrs,
                              const TensorVector& vec)
{}

// no support for checkpoint tensors yet... they should be functions only
// anyway at this point and therefore need not be checkpointed.
template <>
inline void
ReadCheckpoint<TensorVector>(const Checkpoint& chkp,
                             Teuchos::ParameterList& attrs,
                             TensorVector& vec)
{}


// ======================================================================
// Specializations for Teuchos::Array<double>
// ======================================================================
template <>
void
WriteVis<Teuchos::Array<double>>(const Visualization& vis,
                                 Teuchos::ParameterList& attrs,
                                 const Teuchos::Array<double>& vec);

template <>
void
WriteCheckpoint<Teuchos::Array<double>>(const Checkpoint& chkp,
                                        Teuchos::ParameterList& attrs,
                                        const Teuchos::Array<double>& vec);

template <>
void
ReadCheckpoint<Teuchos::Array<double>>(const Checkpoint& chkp,
                                       Teuchos::ParameterList& attrs,
                                       Teuchos::Array<double>& vec);

template <>
bool
Initialize<Teuchos::Array<double>>(Teuchos::ParameterList& plist, Teuchos::Array<double>& t);


// ======================================================================
// Specializations for Patch<double>, Patch<int>
// ======================================================================
template <>
inline void
WriteVis<MultiPatch<double>>(const Visualization& vis,
                             Teuchos::ParameterList& attrs,
                             const MultiPatch<double>& vec)
{}

template <>
inline void
WriteCheckpoint<MultiPatch<double>>(const Checkpoint& chkp,
                                    Teuchos::ParameterList& attrs,
                                    const MultiPatch<double>& vec)
{}

template <>
inline void
ReadCheckpoint<MultiPatch<double>>(const Checkpoint& chkp,
                                   Teuchos::ParameterList& attrs,
                                   MultiPatch<double>& vec)
{}

template <>
inline bool
Initialize<MultiPatch<double>>(Teuchos::ParameterList& plist, MultiPatch<double>& t)
{
  return true;
}


template <>
inline void
WriteVis<MultiPatch<int>>(const Visualization& vis,
                          Teuchos::ParameterList& attrs,
                          const MultiPatch<int>& vec)
{}

template <>
inline void
WriteCheckpoint<MultiPatch<int>>(const Checkpoint& chkp,
                                 Teuchos::ParameterList& attrs,
                                 const MultiPatch<int>& vec)
{}

template <>
inline void
ReadCheckpoint<MultiPatch<int>>(const Checkpoint& chkp,
                                Teuchos::ParameterList& attrs,
                                MultiPatch<int>& vec)
{}

template <>
inline bool
Initialize<MultiPatch<int>>(Teuchos::ParameterList& plist, MultiPatch<int>& t)
{
  return true;
}


template <>
inline void
WriteVis<Operators::BCs>(const Visualization& vis,
                         Teuchos::ParameterList& attrs,
                         const Operators::BCs& vec)
{}

template <>
inline void
WriteCheckpoint<Operators::BCs>(const Checkpoint& chkp,
                                Teuchos::ParameterList& attrs,
                                const Operators::BCs& vec)
{}

template <>
inline void
ReadCheckpoint<Operators::BCs>(const Checkpoint& chkp,
                               Teuchos::ParameterList& attrs,
                               Operators::BCs& vec)
{}

template <>
inline bool
Initialize<Operators::BCs>(Teuchos::ParameterList& plist, Operators::BCs& t)
{
  return true;
}


} // namespace Helpers
} // namespace Amanzi

#endif
