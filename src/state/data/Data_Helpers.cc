/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
Amanzi

License: BSD, see $AMANZI_DIR/COPYRIGHT
Author: Ethan Coon

Helpers that know how to read/write/etc data.

------------------------------------------------------------------------- */

#include "Data_Helpers.hh"

#include "CompositeVectorFunction.hh"
#include "CompositeVectorFunctionFactory.hh"

namespace Amanzi {
namespace Helpers {

//
// Specializations for simple data types
// ======================================================================

template <>
void WriteVis<double>(const Visualization &vis, const Key &fieldname,
                      const std::vector<std::string> &subfieldnames,
                      const double &t) {
  vis.Write(fieldname, t);
}

template <>
void WriteCheckpoint<double>(const Checkpoint &chkp, const Key &fieldname,
                             const double &t) {
  chkp.Write(fieldname, t);
}

template <>
void ReadCheckpoint<double>(const Checkpoint &chkp, const Key &fieldname,
                            double &t) {
  chkp.Read(fieldname, t);
}

// specialization for double
template <>
bool Initialize<double>(Teuchos::ParameterList &plist, double &t,
                        const Key &fieldname,
                        const std::vector<std::string> &subfieldnames) {
  // initialize by list-provided value
  if (plist.isParameter("value")) {
    t = plist.get<double>("value");
    return true;
  }
  return false;
}

template <>
void WriteVis<int>(const Visualization &vis, const Key &fieldname,
                   const std::vector<std::string> &subfieldnames,
                   const int &t) {
  vis.Write(fieldname, t);
}

template <>
void WriteCheckpoint<int>(const Checkpoint &chkp, const Key &fieldname,
                          const int &t) {
  chkp.Write(fieldname, t);
}

template <>
void ReadCheckpoint<int>(const Checkpoint &chkp, const Key &fieldname, int &t) {
  chkp.Read(fieldname, t);
}

// specialization for int
template <>
bool Initialize<int>(Teuchos::ParameterList &plist, int &t,
                     const Key &fieldname,
                     const std::vector<std::string> &subfieldnames) {
  // initialize by list-provided value
  if (plist.isParameter("value")) {
    t = plist.get<int>("value");
    return true;
  }
  return false;
}

//
// Specializations for CompositeVector
// ======================================================================

template <>
void WriteVis<CompositeVector>(const Visualization &vis, const Key &fieldname,
                               const std::vector<std::string> &subfieldnames,
                               const CompositeVector &vec) {
  return;
  // if (vec.HasComponent("cell")) {
  //   // FIXME EPETRA TO TPETRA: need HostView of DeviceView vector
  //   auto vec_c = vec.ViewComponent("cell", false);
  //   if (subfieldnames.size() > 0) {
  //     if (vec.getNumVectors("cell") != subfieldnames.size()) {
  //       Errors::Message msg;
  //       msg << "While Visualizing \"" << fieldname << "\" a vector of lenth "
  //           << vec.getNumVectors("cell") << ", subfieldnames of length "
  //           << (int)subfieldnames.size() << " were provided.";
  //       throw(msg);
  //     }

  //     // I feel like this could be fixed with boost zip_iterator, but the
  //     // Epetra vectors need iterated over and they aren't, by default,
  //     // iterable.
  //     for (int i = 0; i != vec.getNumVectors("cell"); ++i) {
  //       // FIXME EPETRA TO TPETRA
  //       //vis.Write(subfieldnames[i], *vec_c(i));
  //     }

  //   } else {
  //     for (int i = 0; i != vec.getNumVectors("cell"); ++i) {
  //       std::stringstream name;
  //       name << fieldname << ".cell." << i;
  //       // FIXME EPETRA TO TPETRA: write HostView to disk
  //       // vis.Write(name.str(), *vec_c(i));
  //     }
  //   }
  // }
}

template <>
void WriteCheckpoint<CompositeVector>(const Checkpoint &chkp,
                                      const Key &fieldname,
                                      const CompositeVector &vec) {
  return;
  // for (const auto &cname : vec) {
  //   // FIXME EPETRA TO TPETRA: need HostView of DeviceView vector
  //   auto vec_c = vec.ViewComponent(cname, false);

  //   for (int i = 0; i != vec.getNumVectors(cname); ++i) {
  //     std::stringstream name;
  //     name << fieldname << "." << cname << "." << i;
  //     // FIXME EPETRA TO TPETRA: write HostView to disk
  //     // chkp.Write(name.str(), *vec_c(i));
  //   }
  // }
}

template <>
void ReadCheckpoint<CompositeVector>(const Checkpoint &chkp,
                                     const Key &fieldname,
                                     CompositeVector &vec) {
  return;
  // for (const auto &cname : vec) {
  //   // FIXME EPETRA TO TPETRA: need HostView of DeviceView vector
  //   auto vec_c = vec.ViewComponent(cname, false);

  //   for (int i = 0; i != vec.getNumVectors(cname); ++i) {
  //     std::stringstream name;
  //     name << fieldname << "." << cname << "." << i;
  //     // FIXME EPETRA TO TPETRA
  //     //chkp.Read(name.str(), *vec_c(i));
  //   }
  // }
}

template <>
bool Initialize<CompositeVector>(
    Teuchos::ParameterList &plist, CompositeVector &t, const Key &fieldname,
    const std::vector<std::string> &subfieldnames) {
  bool fully_initialized = false;

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

  // ------ Set values using a constant -----
  if (plist.isParameter("constant")) {
    double value = plist.get<double>("constant");
    t.putScalar(value);
    return true;
  }

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
  //     AMANZI_ASSERT(t.getNumVectors("face") == 1);              // and is scalar

  //     // create a vector on faces of the appropriate dimension
  //     int dim = t.Mesh()->space_dimension();

  //     CompositeVectorSpace cvs;
  //     cvs.SetMesh(t.Mesh());
  //     cvs.SetComponent("face", AmanziMesh::FACE, dim);
  //     CompositeVector vel_vec(cvs);

  //     // Evaluate the velocity function
  //     auto func = Functions::CreateCompositeVectorFunction(func_plist, vel_vec.getMap());
  //     func->Compute(0.0, vel_vec);

  //     // Dot the velocity with the normal
  //     unsigned int nfaces_owned =
  //         t.Mesh()->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);

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

  if (fully_initialized) {
    if ((t.HasComponent("face") || t.HasComponent("boundary_face")) &&
        t.HasComponent("cell") &&
        plist.get<bool>("initialize faces from cells", false)) {
      // FIXME EPETRA TO TPETRA
      // DeriveFaceValuesFromCellValues(t);
      return true;
    }
  }
  return fully_initialized;
}

} // namespace Helpers
} // namespace Amanzi
