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

#include "Data_Helpers.hh"

#include "CompositeVector.hh"
#include "TensorVector.hh"
#include "BoundaryFunction.hh"
#include "TreeVectorSpace.hh"
#include "TreeVector.hh"
#include "TensorVector.hh"

#include "ColumnMeshFunction.hh"
#include "CompositeVectorFunction.hh"
#include "CompositeVectorFunctionFactory.hh"
#include "IO.hh"

namespace Amanzi {
namespace Helpers {

// ======================================================================
// Specializations for simple data types
// ======================================================================
template <>
void
WriteVis<double>(const Visualization& vis,
                 const Key& fieldname,
                 const std::vector<std::string>* subfieldnames,
                 const double& t)
{
  vis.Write(fieldname, t);
}

template <>
void
WriteCheckpoint<double>(const Checkpoint& chkp,
                        const Key& fieldname,
                        const std::vector<std::string>* subfieldnames,
                        const double& t)
{
  chkp.Write(fieldname, t);
}

template <>
bool
ReadCheckpoint<double>(const Checkpoint& chkp,
                       const Key& fieldname,
                       const std::vector<std::string>* subfieldnames,
                       double& t)
{
  return chkp.Read(fieldname, t);
}

template <>
bool
Initialize<double>(Teuchos::ParameterList& plist,
                   double& t,
                   const Key& fieldname,
                   const std::vector<std::string>* subfieldnames)
{
  // initialize by list-provided value
  if (plist.isParameter("value")) {
    t = plist.get<double>("value");
    return true;
  }
  return false;
}


template <>
void
WriteVis<int>(const Visualization& vis,
              const Key& fieldname,
              const std::vector<std::string>* subfieldnames,
              const int& t)
{
  vis.Write(fieldname, t);
}

template <>
void
WriteCheckpoint<int>(const Checkpoint& chkp,
                     const Key& fieldname,
                     const std::vector<std::string>* subfieldnames,
                     const int& t)
{
  chkp.Write(fieldname, t);
}

template <>
bool
ReadCheckpoint<int>(const Checkpoint& chkp,
                    const Key& fieldname,
                    const std::vector<std::string>* subfieldnames,
                    int& t)
{
  return chkp.Read(fieldname, t);
}

template <>
bool
Initialize<int>(Teuchos::ParameterList& plist,
                int& t,
                const Key& fieldname,
                const std::vector<std::string>* subfieldnames)
{
  // initialize by list-provided value
  if (plist.isParameter("value")) {
    t = plist.get<int>("value");
    return true;
  }
  return false;
}


// ======================================================================
// Specializations for geometric objects
// ======================================================================
template <>
void
WriteVis<AmanziGeometry::Point>(const Visualization& vis,
                                const Key& fieldname,
                                const std::vector<std::string>* subfieldnames,
                                const AmanziGeometry::Point& vec)
{}

template <>
bool
Initialize<AmanziGeometry::Point>(Teuchos::ParameterList& plist,
                                  AmanziGeometry::Point& p,
                                  const Key& fieldname,
                                  const std::vector<std::string>* subfieldnames)
{
  if (plist.isParameter("value")) {
    auto tmp = plist.get<Teuchos::Array<double>>("value").toVector();
    if (tmp.size() < 4) {
      p.set(tmp.size(), tmp.data());
      return true;
    }
  }
  return false;
}


// ======================================================================
// Specializations for CompositeVector
// ======================================================================
template <>
void
WriteVis<CompositeVector>(const Visualization& vis,
                          const Key& fieldname,
                          const std::vector<std::string>* subfieldnames,
                          const CompositeVector& vec)
{
  if (vec.HasComponent("cell")) {
    const auto& vec_c = *vec.ViewComponent("cell");
    if (subfieldnames && subfieldnames->size() > 0) {
      if (vec_c.NumVectors() != subfieldnames->size()) {
        Errors::Message msg;
        msg << "While Visualizing \"" << fieldname << "\" a vector of lenth " << vec_c.NumVectors()
            << ", subfieldnames of length " << (int)subfieldnames->size() << " were provided.";
        throw(msg);
      }

      std::vector<Key> full_names;
      for (int i = 0; i != vec_c.NumVectors(); ++i) {
        Key full_name = fieldname + "." + (*subfieldnames)[i];
        full_names.emplace_back(full_name);
      }
      vis.WriteVector(vec_c, full_names);

    } else {
      std::vector<Key> full_names;
      if (vec_c.NumVectors() > 1) {
        for (int i = 0; i != vec_c.NumVectors(); ++i) {
          Key full_name = fieldname + "." + std::to_string(i);
          full_names.emplace_back(full_name);
        }
      } else {
        full_names.emplace_back(fieldname);
      }
      vis.WriteVector(vec_c, full_names);
    }
  }
}

template <>
void
WriteCheckpoint<CompositeVector>(const Checkpoint& chkp,
                                 const Key& fieldname,
                                 const std::vector<std::string>* subfieldnames,
                                 const CompositeVector& vec)
{
  for (const auto& cname : vec) {
    const auto& vec_c = *vec.ViewComponent(cname, false);

    if (subfieldnames && subfieldnames->size() > 0) {
      if (subfieldnames && vec_c.NumVectors() != subfieldnames->size()) {
        Errors::Message msg;
        msg << "While Visualizing \"" << fieldname << "\" a vector of lenth " << vec_c.NumVectors()
            << ", subfieldnames of length " << (int)subfieldnames->size() << " were provided.";
        throw(msg);
      }

      std::vector<std::string> comp_subfieldnames(vec_c.NumVectors());
      for (int i = 0; i != vec_c.NumVectors(); ++i) {
        comp_subfieldnames[i] = fieldname + '.' + cname + '.' + (*subfieldnames)[i];
      }
      chkp.WriteVector(vec_c, comp_subfieldnames);
    } else {
      std::vector<std::string> comp_subfieldnames(vec_c.NumVectors());
      for (int i = 0; i != vec_c.NumVectors(); ++i) {
        std::stringstream name;
        name << fieldname << "." << cname << "." << i;
        comp_subfieldnames[i] = name.str();
      }
      chkp.WriteVector(vec_c, comp_subfieldnames);
    }
  }
}

template <>
bool
ReadCheckpoint<CompositeVector>(const Checkpoint& chkp,
                                const Key& fieldname,
                                const std::vector<std::string>* subfieldnames,
                                CompositeVector& vec)
{
  bool flag(true);
  for (const auto& cname : vec) {
    auto& vec_c = *vec.ViewComponent(cname, false);

    if (subfieldnames && subfieldnames->size() > 0) {
      if (vec_c.NumVectors() != subfieldnames->size()) {
        Errors::Message msg;
        msg << "While Visualizing \"" << fieldname << "\" a vector of lenth " << vec_c.NumVectors()
            << ", subfieldnames of length " << (int)subfieldnames->size() << " were provided.";
        throw(msg);
      }

      for (int i = 0; i != vec_c.NumVectors(); ++i) {
        Key name = fieldname + '.' + cname + '.' + (*subfieldnames)[i];
        chkp.Read(name, *vec_c(i));
      }
    } else {
      for (int i = 0; i != vec_c.NumVectors(); ++i) {
        std::stringstream name;
        name << fieldname << "." << cname << "." << i;
        flag &= chkp.Read(name.str(), *vec_c(i));
      }
    }
  }
  return flag;
}

template <>
bool
Initialize<CompositeVector>(Teuchos::ParameterList& plist,
                            CompositeVector& t,
                            const Key& fieldname,
                            const std::vector<std::string>* subfieldnames)
{
  bool initialized = false;

  // First try all initialization method which set the entire data structure.
  // ------ Try to set values from a restart file -----
  if (plist.isParameter("restart file")) {
    auto filename = plist.get<std::string>("restart file");
    Checkpoint chkp(filename, t.Comm(), Keys::getDomain(fieldname));
    ReadCheckpoint(chkp, fieldname, subfieldnames, t);
    chkp.Finalize();
    initialized = true;
  }

  // // ------ Try to set values from an file -----
  if (plist.isSublist("exodus file initialization")) {
    // data must be pre-initialized to zero in case Exodus file does not
    // provide all values.
    t.PutScalar(0.0);

    Teuchos::ParameterList& file_list = plist.sublist("exodus file initialization");
    ReadVariableFromExodusII(file_list, t);
    initialized = true;
  }

  // ------ Set values using a constant -----
  if (plist.isParameter("constant")) {
    double value = plist.get<double>("constant");
    t.PutScalar(value);
    initialized = true;
  }
  if (plist.isParameter("value")) {
    double value = plist.get<double>("value");
    t.PutScalar(value);
    initialized = true;
  }

  // Next try all partial initialization methods -- typically cells.
  // ------ Try to set cell values from a restart file -----
  if (plist.isParameter("cells from file")) {
    AMANZI_ASSERT(false);
    auto filename = plist.get<std::string>("cells from file");
    Checkpoint chkp(filename, t.Comm());

    // read just the cells
    auto& vec_c = *t.ViewComponent("cell");
    for (int i = 0; i != vec_c.NumVectors(); ++i) {
      std::stringstream name;
      name << fieldname << ".cell." << i;
      chkp.Read(name.str(), *vec_c(i));
    }
    chkp.Finalize();
    initialized = true;
  }


  // ------ Set values from 1D solution -----
  if (plist.isSublist("initialize from 1D column")) {
    Teuchos::ParameterList& init_plist = plist.sublist("initialize from 1D column");
    if (!init_plist.isParameter("f header"))
      init_plist.set("f header", std::string("/") + fieldname);
    Functions::ReadColumnMeshFunction(init_plist, t);
    initialized = true;
  }

  // ------ Set values using a function -----
  if (plist.isSublist("function")) {
    double t_ini = plist.get<double>("time", 0.0);
    Teuchos::ParameterList func_plist = plist.sublist("function");

    std::vector<std::string> complist;

    // -- potential use of a mapping operator first --
    bool map_normal = plist.get<bool>("dot with normal", false);
    if (map_normal) {
      // map_normal take a vector and dots it with face normals
      AMANZI_ASSERT(t.NumComponents() == 1);                 // one comp
      AMANZI_ASSERT(t.HasComponent("face"));                 // is named face
      AMANZI_ASSERT(t.Location("face") == AmanziMesh::Entity_kind::FACE); // is on face
      AMANZI_ASSERT(t.NumVectors("face") == 1);              // and is scalar

      // create a vector on faces of the appropriate dimension
      int dim = t.Mesh()->getSpaceDimension();

      CompositeVectorSpace cvs;
      cvs.SetMesh(t.Mesh());
      cvs.SetComponent("face", AmanziMesh::Entity_kind::FACE, dim);
      Teuchos::RCP<CompositeVector> vel_vec = Teuchos::rcp(new CompositeVector(cvs));

      // Evaluate the velocity function
      auto func = Functions::CreateCompositeVectorFunction(func_plist, vel_vec->Map(), complist);
      func->Compute(t_ini, vel_vec.ptr());

      // CV's map may differ from the regular mesh map due to presense of fractures
      const auto& fmap = *t.Map().Map("face", true);
      int nfaces_owned = t.Mesh()->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);

      Epetra_MultiVector& dat_f = *t.ViewComponent("face");
      const Epetra_MultiVector& vel_f = *vel_vec->ViewComponent("face");

      // Dot the velocity with the normal
      // -- two branches: single flux per face, multiple fluxes
      int dir;
      AmanziGeometry::Point vel(dim);
      for (int f = 0; f != nfaces_owned; ++f) {
        for (int i = 0; i < dim; ++i) vel[i] = vel_f[i][f];

        int ndofs = fmap.ElementSize(f);
        int g = fmap.FirstPointInElement(f);
        if (ndofs == 1) {
          const AmanziGeometry::Point& normal = t.Mesh()->getFaceNormal(f);
          dat_f[0][g] = vel * normal;
        } else {
          AmanziMesh::Entity_ID_View cells;
          cells = t.Mesh()->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);

          for (int i = 0; i < ndofs; ++i) {
            const AmanziGeometry::Point& normal = t.Mesh()->getFaceNormal(f, cells[i], &dir);
            dat_f[0][g + i] = (vel * normal) * dir;
          }
        }
      }
      initialized = true;

    } else {
      auto t_ptr = Teuchos::rcpFromRef(t).ptr();

      // no map, just evaluate the function
      auto func = Functions::CreateCompositeVectorFunction(func_plist, t.Map(), complist);
      func->Compute(t_ini, t_ptr);
      initialized = true;
    }
  }

  if ((t.HasComponent("face") || t.HasComponent("boundary_face")) && t.HasComponent("cell") &&
      plist.get<bool>("initialize faces from cells", false)) {
    DeriveFaceValuesFromCellValues(t);
  }
  return initialized;
}


// ======================================================================
// Specializations for Epetra_Vector
// ======================================================================
template <>
void
WriteVis<Epetra_Vector>(const Visualization& vis,
                        const Key& fieldname,
                        const std::vector<std::string>* subfieldnames,
                        const Epetra_Vector& vec)
{
  vis.Write(fieldname, vec);
}

template <>
void
WriteCheckpoint<Epetra_Vector>(const Checkpoint& chkp,
                               const Key& fieldname,
                               const std::vector<std::string>* subfieldnames,
                               const Epetra_Vector& vec)
{
  chkp.Write(fieldname, vec);
}

template <>
bool
ReadCheckpoint<Epetra_Vector>(const Checkpoint& chkp,
                              const Key& fieldname,
                              const std::vector<std::string>* subfieldnames,
                              Epetra_Vector& vec)
{
  return chkp.Read(fieldname, vec);
}

template <>
bool
Initialize<Epetra_Vector>(Teuchos::ParameterList& plist,
                          Epetra_Vector& t,
                          const Key& fieldname,
                          const std::vector<std::string>* subfieldnames)
{
  bool initialized = false;

  // First try all initialization method which set the entire data structure.
  // ------ Try to set values from a restart file -----
  if (plist.isParameter("restart file")) {
    auto filename = plist.get<std::string>("restart file");
    auto comm = Teuchos::rcpFromRef(t.Comm());
    Checkpoint chkp(filename, comm);
    ReadCheckpoint(chkp, fieldname, subfieldnames, t);
    chkp.Finalize();
    initialized = true;
  }

  // ------ Set values using a constant -----
  if (plist.isParameter("constant")) {
    double value = plist.get<double>("constant");
    t.PutScalar(value);
    initialized = true;
  }
  if (plist.isParameter("value")) {
    double value = plist.get<double>("value");
    t.PutScalar(value);
    initialized = true;
  }
  return initialized;
}


template <>
bool
Equivalent(const Epetra_Map& one, const Epetra_Map& two)
{
  return one.SameAs(two);
}

template <>
bool
Equivalent(const Epetra_BlockMap& one, const Epetra_BlockMap& two)
{
  return one.SameAs(two);
}

template <>
bool
Equivalent(const TreeVectorSpace& one, const TreeVectorSpace& two)
{
  return one.SameAs(two);
}


// ======================================================================
// Specializations for TreeVector
// ======================================================================
template <>
void
WriteVis<TreeVector>(const Visualization& vis,
                     const Key& fieldname,
                     const std::vector<std::string>* subfieldnames,
                     const TreeVector& vec)
{
  int i = 0;
  for (const auto& subvec : vec) {
    std::string subvec_name = fieldname + "_" + std::to_string(i);
    WriteVis(vis, subvec_name, nullptr, *subvec);
    ++i;
  }
  if (vec.Data() != Teuchos::null) { WriteVis(vis, fieldname + "_data", nullptr, *vec.Data()); }
}

template <>
void
WriteCheckpoint<TreeVector>(const Checkpoint& chkp,
                            const Key& fieldname,
                            const std::vector<std::string>* subfieldnames,
                            const TreeVector& vec)
{
  int i = 0;
  for (const auto& subvec : vec) {
    std::string subvec_name = fieldname + "_" + std::to_string(i);
    WriteCheckpoint(chkp, subvec_name, nullptr, *subvec);
    ++i;
  }
  if (vec.Data() != Teuchos::null) {
    WriteCheckpoint(chkp, fieldname + "_data", nullptr, *vec.Data());
  }
}

template <>
bool
ReadCheckpoint<TreeVector>(const Checkpoint& chkp,
                           const Key& fieldname,
                           const std::vector<std::string>* subfieldnames,
                           TreeVector& vec)
{
  bool flag(true);
  int i = 0;
  for (const auto& subvec : vec) {
    std::string subvec_name = fieldname + "_" + std::to_string(i);
    flag &= ReadCheckpoint(chkp, subvec_name, nullptr, *subvec);
    ++i;
  }
  if (vec.Data() != Teuchos::null) {
    flag &= ReadCheckpoint(chkp, fieldname + "_data", nullptr, *vec.Data());
  }
  return flag;
}

template <>
bool
Initialize<TreeVector>(Teuchos::ParameterList& plist,
                       TreeVector& t,
                       const Key& fieldname,
                       const std::vector<std::string>* subfieldnames)
{
  return false;
}


// ======================================================================
// Specializations for Teuchos::Array<double>
// ======================================================================
template <>
void
WriteVis<Teuchos::Array<double>>(const Visualization& vis,
                                 const Key& fieldname,
                                 const std::vector<std::string>* subfieldnames,
                                 const Teuchos::Array<double>& vec)
{
  if (subfieldnames) {
    AMANZI_ASSERT(vec.size() != subfieldnames->size());
    for (int i = 0; i != vec.size(); ++i) {
      vis.Write(fieldname + "_" + (*subfieldnames)[i], vec[i]);
    }
  } else {
    for (int i = 0; i != vec.size(); ++i) {
      vis.Write(fieldname + "_" + std::to_string(i), vec[i]);
    }
  }
}

template <>
void
WriteCheckpoint<Teuchos::Array<double>>(const Checkpoint& chkp,
                                        const Key& fieldname,
                                        const std::vector<std::string>* subfieldnames,
                                        const Teuchos::Array<double>& vec)
{
  if (subfieldnames) {
    AMANZI_ASSERT(vec.size() != subfieldnames->size());
    for (int i = 0; i != vec.size(); ++i) {
      chkp.Write(fieldname + "_" + (*subfieldnames)[i], vec[i]);
    }
  } else {
    for (int i = 0; i != vec.size(); ++i) {
      chkp.Write(fieldname + "_" + std::to_string(i), vec[i]);
    }
  }
}

template <>
bool
ReadCheckpoint<Teuchos::Array<double>>(const Checkpoint& chkp,
                                       const Key& fieldname,
                                       const std::vector<std::string>* subfieldnames,
                                       Teuchos::Array<double>& vec)
{
  bool flag = true;
  if (subfieldnames) {
    AMANZI_ASSERT(vec.size() != subfieldnames->size());
    for (int i = 0; i != vec.size(); ++i) {
      flag &= chkp.Read(fieldname + "_" + (*subfieldnames)[i], vec[i]);
    }
  } else {
    for (int i = 0; i != vec.size(); ++i) {
      flag &= chkp.Read(fieldname + "_" + std::to_string(i), vec[i]);
    }
  }
  return flag;
}


template <>
bool
Initialize<Teuchos::Array<double>>(Teuchos::ParameterList& plist,
                                   Teuchos::Array<double>& t,
                                   const Key& fieldname,
                                   const std::vector<std::string>* subfieldnames)
{
  if (plist.isParameter("value")) {
    auto arr = plist.get<Teuchos::Array<double>>("value");
    if (t.size() != arr.size()) {
      Errors::Message msg;
      msg << "While initializing \"" << fieldname << "\" an array of lenth " << (int)t.size()
          << " was expected, but an array of length " << (int)arr.size() << " was provided.";
      throw(msg);
    }
    t = arr;
    return true;
  }
  return false;
}


} // namespace Helpers
} // namespace Amanzi
