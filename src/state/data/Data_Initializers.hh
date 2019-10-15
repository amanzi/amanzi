/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! Data_Initializers are functions that know how to initialize data.

/*
  Initialize from parameter list.
*/

#ifndef AMANZI_DATA_INITIALIZERS_HH_
#define AMANZI_DATA_INITIALIZERS_HH_

#include "Teuchos_ParameterList.hpp"

#include "errors.hh"
#include "FunctionFactory.hh"

#include "CompositeVector.hh"
//#include "Op.hh"
//#include "TensorVector.hh"
//#include "Operator.hh"
//#include "BCs.hh"
//#include "BoundaryFunction.hh"

namespace Amanzi {
namespace Data_Initializers {

//
// Initialize
// ======================================================================

// This fails to compile if not provided by user.
template <typename T>
bool
Initialize(Teuchos::ParameterList& plist, const Teuchos::ParameterList& attrs,
           T& t)
{
  return UserInitialize(plist, attrs, t);
}


// Initialize from a single parameter value.
template <typename T>
bool
InitializePrimitiveByValue(Teuchos::ParameterList& plist,
                           const Teuchos::ParameterList& attrs, T& t)
{
  if (plist.isParameter("value")) {
    t = plist.template get<T>("value");
    return true;
  }
  return false;
}

// Initialize from a single parameter value.
template <template <typename> class V, typename T>
bool
InitializeVectorByValue(Teuchos::ParameterList& plist,
                        const Teuchos::ParameterList& attrs, V<T>& t)
{
  if (plist.isParameter("value")) {
    t.putScalar(plist.template get<T>("value"));
    return true;
  }
  return false;
}


//
// Specializations for simple data types
// ======================================================================
template <>
bool
Initialize<double>(Teuchos::ParameterList& plist,
                   const Teuchos::ParameterList& attrs, double& t);
template <>
bool
Initialize<int>(Teuchos::ParameterList& plist,
                const Teuchos::ParameterList& attrs, int& t);
template <>
bool
Initialize<std::string>(Teuchos::ParameterList& plist,
                        const Teuchos::ParameterList& attrs, std::string& t);

//
// Specializations for Vectors
// ======================================================================
template <typename Scalar>
bool
Initialize(Teuchos::ParameterList& plist, const Teuchos::ParameterList& attrs,
           Vector_type_<Scalar>& t)
{
  if (InitializeVectorByValue<Vector_type_, Scalar>(plist, attrs, t))
    return true;
  return false;
}

template <typename Scalar>
bool
Initialize(Teuchos::ParameterList& plist, const Teuchos::ParameterList& attrs,
           MultiVector_type_<Scalar>& t)
{
  if (InitializeVectorByValue<MultiVector_type_, Scalar>(plist, attrs, t))
    return true;
  return false;
}

template <typename Scalar>
bool
Initialize(Teuchos::ParameterList& plist, const Teuchos::ParameterList& attrs,
           CompositeVector_<Scalar>& t)
{
  if (InitializeVectorByValue<CompositeVector_, Scalar>(plist, attrs, t))
    return true;
  return false;
}


//
// Specializations for Vectors
// ======================================================================
// template <>
// bool Initialize<CompositeVector>(
//     Teuchos::ParameterList &plist, CompositeVector &t, const Key &fieldname,
//     const std::vector<std::string> &subfieldnames) {
//   bool fully_initialized = false;

// First try all initialization method which set the entire data structure.
// ------ Try to set values from a restart file -----
// FIXME EPETRA TO TPETRA: read restart
// if (plist.isParameter("restart file")) {
//   auto filename = plist.get<std::string>("restart file");
//   Checkpoint chkp(filename, t.Comm());
//   ReadCheckpoint(chkp, fieldname, t);
//   chkp.Finalize();
//   return true;
// }

// // ------ Try to set values from an file -----
// if (plist.isSublist("exodus file initialization")) {
//   // data must be pre-initialized to zero in case Exodus file does not
//   // provide all values.
//   t.putScalar(0.0);

//   Teuchos::ParameterList& file_list = plist.sublist("exodus file
//   initialization"); Functions::ReadExodusIIMeshFunction(file_list, t);
//   return true;
// }

// // ------ Set values using a constant -----
// if (plist.isParameter("constant")) {
//   double value = plist.get<double>("constant");
//   t.putScalar(value);
//   return true;
// }

// Next try all partial initialization methods -- typically cells.
// ------ Try to set cell values from a restart file -----
// FIXME EPETRA TO TPETRA: read file
// if (plist.isParameter("cells from file")) {
//   auto filename = plist.get<std::string>("cells from file");
//   Checkpoint chkp(filename, t.Comm());

//   // read just the cells
//   auto vec_c = t.ViewComponent("cell", false);
//   for (int i = 0; i != t.getNumVectors("cell"); ++i) {
//     std::stringstream name;
//     name << fieldname << ".cell." << i;
//     // FIXME EPETRA TO TPETRA
//     // chkp.Read(name.str(), *vec_c(i));
//   }
//   chkp.Finalize();
// }

// ------ Set values from 1D solution -----
// FIXME EPETRA TO TPETRA: read restart
// if (plist.isSublist("initialize from 1D column")) {
//   Teuchos::ParameterList &init_plist =
//       plist.sublist("initialize from 1D column");
//   if (!init_plist.isParameter("f header"))
//     init_plist.set("f header", std::string("/") + fieldname);
//   Functions::ReadColumnMeshFunction(init_plist, t);
//   fully_initialized = true;
// }

// ------ Set values using a function -----
// FIXME EPETRA TO TPETRA: mesh_functions doesn't compile
// if (plist.isSublist("function")) {
//   Teuchos::ParameterList func_plist = plist.sublist("function");

//   // -- potential use of a mapping operator first --
//   bool map_normal = plist.get<bool>("dot with normal", false);
//   if (map_normal) {
//     // FIXME EPETRA TO TPETRA
//     // map_normal take a vector and dots it with face normals
//     AMANZI_ASSERT(t.NumComponents() == 1);                 // one comp
//     AMANZI_ASSERT(t.HasComponent("face"));                 // is named face
//     AMANZI_ASSERT(t.Location("face") == AmanziMesh::FACE); // is on face
//     AMANZI_ASSERT(t.getNumVectors("face") == 1);              // and is
//     scalar

//     // create a vector on faces of the appropriate dimension
//     int dim = t.Mesh()->space_dimension();

//     CompositeVectorSpace cvs;
//     cvs.SetMesh(t.Mesh());
//     cvs.SetComponent("face", AmanziMesh::FACE, dim);
//     CompositeVector vel_vec(cvs);

//     // Evaluate the velocity function
//     auto func = Functions::CreateCompositeVectorFunction(func_plist,
//     vel_vec.getMap()); func->Compute(0.0, vel_vec);

//     // Dot the velocity with the normal
//     unsigned int nfaces_owned =
//         t.Mesh()->num_entities(AmanziMesh::FACE,
//         AmanziMesh::Parallel_type::OWNED);

//     auto dat_f = t.ViewComponent("face", false);
//     auto vel_f = vel_vec.ViewComponent("face", false);

//     AmanziGeometry::Point vel(dim);
//     for (unsigned int f = 0; f != nfaces_owned; ++f) {
//       AmanziGeometry::Point normal = t.Mesh()->face_normal(f);
//       if (dim == 2) {
//         vel.set(vel_f[0][f], vel_f[1][f]);
//       } else if (dim == 3) {
//         vel.set(vel_f[0][f], vel_f[1][f], vel_f[2][f]);
//       } else {
//         AMANZI_ASSERT(0);
//       }
//       dat_f[0][f] = vel * normal;
//     }
//     return true;

//   } else {
//     // FIXME EPETRA TO TPETRA
//     // // no map, just evaluate the function
//     // Teuchos::RCP<Functions::CompositeVectorFunction> func =
//     //     Functions::CreateCompositeVectorFunction(func_plist, t.getMap());
//     // func->Compute(0.0, t);
//     // fully_initialized = true;
//   }
// }

//   if (fully_initialized) {
//     if ((t.HasComponent("face") || t.HasComponent("boundary_face")) &&
//         t.HasComponent("cell") &&
//         plist.get<bool>("initialize faces from cells", false)) {
//       // FIXME EPETRA TO TPETRA
//       // DeriveFaceValuesFromCellValues(t);
//       return true;
//     }
//   }
//   return fully_initialized;
// }


//
// Specializations for Op
// ======================================================================
/*
template <>
inline void
WriteVis<Operators::Op>(const Visualization &vis, const Key &fieldname,
                        const std::vector<std::string> &subfieldnames,
                        const Operators::Op &vec) {}

template <>
inline void WriteCheckpoint<Operators::Op>(const Checkpoint &chkp,
                                           const Key &fieldname,
                                           const Operators::Op &vec) {}

template <>
inline void ReadCheckpoint<Operators::Op>(const Checkpoint &chkp,
                                          const Key &fieldname,
                                          Operators::Op &vec) {}

template <>
inline bool
Initialize<Operators::Op>(Teuchos::ParameterList &plist, Operators::Op &t,
                          const Key &fieldname,
                          const std::vector<std::string> &subfieldnames) {
  return true;
}

//
// Specializations for Operator
// ======================================================================

template <>
inline void
WriteVis<Operators::Operator>(const Visualization &vis, const Key &fieldname,
                              const std::vector<std::string> &subfieldnames,
                              const Operators::Operator &global_operator) {}

template <>
inline void WriteCheckpoint<Operators::Operator>(
    const Checkpoint &chkp, const Key &fieldname,
    const Operators::Operator &global_operator) {}

template <>
inline void
ReadCheckpoint<Operators::Operator>(const Checkpoint &chkp,
                                    const Key &fieldname,
                                    Operators::Operator &global_operator) {}

template <>
inline bool
Initialize<Operators::Operator>(Teuchos::ParameterList &plist,
                                Operators::Operator &t, const Key &fieldname,
                                const std::vector<std::string> &subfieldnames) {
  return true;
}

//
// Specializations for WhetStone::Tensor
// ======================================================================

template <>
inline void
WriteVis<TensorVector>(const Visualization &vis, const Key &fieldname,
                       const std::vector<std::string> &subfieldnames,
                       const TensorVector &tensor) {}

template <>
inline void WriteCheckpoint<TensorVector>(const Checkpoint &chkp,
                                          const Key &fieldname,
                                          const TensorVector &tensor) {}

template <>
inline void ReadCheckpoint<TensorVector>(const Checkpoint &chkp,
                                         const Key &fieldname,
                                         TensorVector &tensor) {}

template <>
inline bool
Initialize<TensorVector>(Teuchos::ParameterList &plist, TensorVector &tensor,
                         const Key &fieldname,
                         const std::vector<std::string> &subfieldnames) {
  return true;
}

//
// Specializations for Operators::BCs
// ======================================================================
*/
/*
template <>
inline void
WriteVis<Operators::BCs>(const Visualization &vis, const Key &fieldname,
                         const std::vector<std::string> &subfieldnames,
                         const Operators::BCs &bc) {}

template <>
inline void WriteCheckpoint<Operators::BCs>(const Checkpoint &chkp,
                                            const Key &fieldname,
                                            const Operators::BCs &bc) {}

template <>
inline void ReadCheckpoint<Operators::BCs>(const Checkpoint &chkp,
                                           const Key &fieldname,
                                           Operators::BCs &bc) {}

template <>
inline bool
Initialize<Operators::BCs>(Teuchos::ParameterList &plist,
                           Operators::BCs &bc, const Key &fieldname,
                           const std::vector<std::string> &subfieldnames) {
  return true;
}


//
// Specializations for Functions::BoundaryFunction
// ======================================================================

template <>
inline void
WriteVis<Functions::BoundaryFunction>(const Visualization &vis, const Key
&fieldname, const std::vector<std::string> &subfieldnames, const
Functions::BoundaryFunction &bc) {}

template <>
inline void WriteCheckpoint<Functions::BoundaryFunction>(const Checkpoint &chkp,
                                            const Key &fieldname,
                                            const Functions::BoundaryFunction
&bc) {}

template <>
inline void ReadCheckpoint<Functions::BoundaryFunction>(const Checkpoint &chkp,
                                           const Key &fieldname,
                                           Functions::BoundaryFunction &bc) {}

template <>
inline bool
Initialize<Functions::BoundaryFunction>(Teuchos::ParameterList &plist,
                           Functions::BoundaryFunction &bc, const Key
&fieldname, const std::vector<std::string> &subfieldnames) { return true;
}
*/
} // namespace Data_Initializers
} // namespace Amanzi

#endif
