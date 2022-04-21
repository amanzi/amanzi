/*
  State

  Copyright 2010-202x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Helpers that know how to read/write/etc data.
*/

#include "ColumnMeshFunction.hh"
#include "CompositeVectorFunction.hh"
#include "CompositeVectorFunctionFactory.hh"

#include "Data_Helpers.hh"
#include "IO.hh"

namespace Amanzi {
namespace Helpers {

// ======================================================================
// Specializations for simple data types
// ======================================================================
template <>
void WriteVis<double>(const Visualization& vis, const Key& fieldname,
                      const std::vector<std::string>& subfieldnames,
                      const double& t) {
  vis.Write(fieldname, t);
}

template <>
void WriteCheckpoint<double>(const Checkpoint& chkp, const Key& fieldname,
                             const std::vector<std::string>& subfieldnames,
                             const double& t) {
  chkp.Write(fieldname, t);
}

template <>
void ReadCheckpoint<double>(const Checkpoint& chkp, const Key& fieldname,
                            const std::vector<std::string>& subfieldnames,
                            double& t) {
  chkp.Read(fieldname, t);
}

template <>
bool Initialize<double>(Teuchos::ParameterList& plist, double& t,
                        const Key& fieldname,
                        const std::vector<std::string>& subfieldnames) {
  // initialize by list-provided value
  if (plist.isParameter("value")) {
    t = plist.get<double>("value");
    return true;
  }
  return false;
}


template <>
void WriteVis<int>(const Visualization& vis, const Key& fieldname,
                   const std::vector<std::string>& subfieldnames,
                   const int& t) {
  vis.Write(fieldname, t);
}

template <>
void WriteCheckpoint<int>(const Checkpoint& chkp, const Key& fieldname,
                          const std::vector<std::string>& subfieldnames,
                          const int& t) {
  chkp.Write(fieldname, t);
}

template <>
void ReadCheckpoint<int>(const Checkpoint& chkp, const Key& fieldname,
                         const std::vector<std::string>& subfieldnames,
                         int& t) {
  chkp.Read(fieldname, t);
}

template <>
bool Initialize<int>(Teuchos::ParameterList& plist, int& t,
                     const Key& fieldname,
                     const std::vector<std::string>& subfieldnames) {
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
void WriteVis<AmanziGeometry::Point>(const Visualization& vis, const Key& fieldname,
                                     const std::vector<std::string>& subfieldnames,
                                     const AmanziGeometry::Point& vec)
{
}

template <>
bool Initialize<AmanziGeometry::Point>(Teuchos::ParameterList& plist,
                                       AmanziGeometry::Point& p, const Key& fieldname,
                                       const std::vector<std::string>& subfieldnames)
{
  if (plist.isParameter("value")) {
    auto tmp = plist.get<Teuchos::Array<double> >("value").toVector();
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
void WriteVis<CompositeVector>(const Visualization& vis, const Key& fieldname,
                               const std::vector<std::string>& subfieldnames,
                               const CompositeVector& vec)
{
  if (vec.HasComponent("cell")) {
    const auto& vec_c = *vec.ViewComponent("cell");
    if (subfieldnames.size() > 0) {
      if (vec_c.NumVectors() != subfieldnames.size()) {
        Errors::Message msg;
        msg << "While Visualizing \"" << fieldname << "\" a vector of lenth "
            << vec_c.NumVectors() << ", subfieldnames of length "
            << (int)subfieldnames.size() << " were provided.";
        throw(msg);
      }

      // I feel like this could be fixed with boost zip_iterator, but the
      // Epetra vectors need iterated over and they aren't, by default,
      // iterable.
      for (int i = 0; i != vec_c.NumVectors(); ++i) {
        Key fullname = fieldname + std::string(".cell.") + subfieldnames[i];
        vis.Write(fullname, *vec_c(i));
      }

    } else {
      for (int i = 0; i != vec_c.NumVectors(); ++i) {
        std::stringstream name;
        name << fieldname << ".cell." << i;
        vis.Write(name.str(), *vec_c(i));
      }
    }
  }
}

template <>
void WriteCheckpoint<CompositeVector>(const Checkpoint& chkp,
                                      const Key& fieldname,
                                      const std::vector<std::string>& subfieldnames,
                                      const CompositeVector& vec) {
  for (const auto& cname : vec) {
    const auto& vec_c = *vec.ViewComponent(cname, false);

    if (subfieldnames.size() > 0) {
      if (vec_c.NumVectors() != subfieldnames.size()) {
        Errors::Message msg;
        msg << "While Visualizing \"" << fieldname << "\" a vector of lenth "
            << vec_c.NumVectors() << ", subfieldnames of length "
            << (int)subfieldnames.size() << " were provided.";
        throw(msg);
      }

      for (int i = 0; i != vec_c.NumVectors(); ++i) {
        Key name = fieldname + '.' + cname + '.' + subfieldnames[i];
        chkp.Write(name, *vec_c(i));
      }
    } else {
      for (int i = 0; i != vec_c.NumVectors(); ++i) {
        std::stringstream name;
        name << fieldname << "." << cname << "." << i;
        chkp.Write(name.str(), *vec_c(i));
      }
    }
  }
}

template <>
void ReadCheckpoint<CompositeVector>(const Checkpoint& chkp,
                                     const Key& fieldname,
                                     const std::vector<std::string>& subfieldnames,
                                     CompositeVector& vec) {
  for (const auto& cname : vec) {
    auto& vec_c = *vec.ViewComponent(cname, false);

    if (subfieldnames.size() > 0) {
      if (vec_c.NumVectors() != subfieldnames.size()) {
        Errors::Message msg;
        msg << "While Visualizing \"" << fieldname << "\" a vector of lenth "
            << vec_c.NumVectors() << ", subfieldnames of length "
            << (int)subfieldnames.size() << " were provided.";
        throw(msg);
      }

      for (int i = 0; i != vec_c.NumVectors(); ++i) {
        Key name = fieldname + '.' + cname + '.' + subfieldnames[i];
        chkp.Read(name, *vec_c(i));
      }
    } else {
      for (int i = 0; i != vec_c.NumVectors(); ++i) {
        std::stringstream name;
        name << fieldname << "." << cname << "." << i;
        chkp.Read(name.str(), *vec_c(i));
      }
    }
  }
}

template <>
bool Initialize<CompositeVector>(
    Teuchos::ParameterList& plist, CompositeVector& t, const Key& fieldname,
    const std::vector<std::string>& subfieldnames) {
  bool fully_initialized = false;

  // First try all initialization method which set the entire data structure.
  // ------ Try to set values from a restart file -----
  if (plist.isParameter("restart file")) {
    auto filename = plist.get<std::string>("restart file");
    Checkpoint chkp(filename, t.Comm());
    ReadCheckpoint(chkp, fieldname, subfieldnames, t);
    chkp.Finalize();
    return true;
  }

  // // ------ Try to set values from an file -----
  if (plist.isSublist("exodus file initialization")) {
    // data must be pre-initialized to zero in case Exodus file does not
    // provide all values.
    t.PutScalar(0.0);

    Teuchos::ParameterList& file_list = plist.sublist("exodus file initialization");
    ReadVariableFromExodusII(file_list, t);
    return true;
  }

  // ------ Set values using a constant -----
  if (plist.isParameter("constant")) {
    double value = plist.get<double>("constant");
    t.PutScalar(value);
    return true;
  }
  if (plist.isParameter("value")) {
    double value = plist.get<double>("value");
    t.PutScalar(value);
    return true;
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
    return true;
  }


  // ------ Set values from 1D solution -----
  if (plist.isSublist("initialize from 1D column")) {
    Teuchos::ParameterList& init_plist = plist.sublist("initialize from 1D column");
    if (!init_plist.isParameter("f header"))
      init_plist.set("f header", std::string("/") + fieldname);
    Functions::ReadColumnMeshFunction(init_plist, t);
    return true;
  }

  // ------ Set values using a function -----
  if (plist.isSublist("function")) {
    Teuchos::ParameterList func_plist = plist.sublist("function");

    std::vector<std::string> complist;

    // -- potential use of a mapping operator first --
    bool map_normal = plist.get<bool>("dot with normal", false);
    if (map_normal) {
      // map_normal take a vector and dots it with face normals
      AMANZI_ASSERT(t.NumComponents() == 1);                 // one comp
      AMANZI_ASSERT(t.HasComponent("face"));                 // is named face
      AMANZI_ASSERT(t.Location("face") == AmanziMesh::FACE); // is on face
      AMANZI_ASSERT(t.NumVectors("face") == 1);              // and is scalar

      // create a vector on faces of the appropriate dimension
      int dim = t.Mesh()->space_dimension();

      CompositeVectorSpace cvs;
      cvs.SetMesh(t.Mesh());
      cvs.SetComponent("face", AmanziMesh::FACE, dim);
      Teuchos::RCP<CompositeVector> vel_vec = Teuchos::rcp(new CompositeVector(cvs));

      // Evaluate the velocity function
      auto func = Functions::CreateCompositeVectorFunction(func_plist, vel_vec->Map(), complist);
      func->Compute(0.0, vel_vec.ptr());

      // CV's map may differ from the regular mesh map due to presense of fractures
      const auto& fmap = *t.Map().Map("face", true);
      int nfaces_owned = t.Mesh()->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);

      Epetra_MultiVector& dat_f = *t.ViewComponent("face");
      const Epetra_MultiVector& vel_f = *vel_vec->ViewComponent("face");

      // Dot the velocity with the normal
      AmanziGeometry::Point vel(dim);
      for (int f = 0; f != nfaces_owned; ++f) {
        const AmanziGeometry::Point& normal = t.Mesh()->face_normal(f);
        for (int i = 0; i < dim; ++i) vel[i] = vel_f[i][f]; 

        int g = fmap.FirstPointInElement(f);
        for (int i = 0; i < fmap.ElementSize(f); ++i) {
          dat_f[0][g + i] = vel * normal; 
        }
      }
      return true;

    } else {
      auto t_ptr = Teuchos::rcpFromRef(t).ptr();

      // no map, just evaluate the function
      auto func = Functions::CreateCompositeVectorFunction(func_plist, t.Map(), complist);
      func->Compute(0.0, t_ptr);
      fully_initialized = true;
    }
  }

  if (fully_initialized) {
    if ((t.HasComponent("face") || t.HasComponent("boundary_face")) &&
         t.HasComponent("cell") &&
         plist.get<bool>("initialize faces from cells", false)) {
      DeriveFaceValuesFromCellValues(t);
      return true;
    }
    if (t.HasComponent("boundary_face") &&
        t.HasComponent("cell") &&
        plist.get<bool>("initialize faces from cells", false)) {
      DeriveFaceValuesFromCellValues(t);
      return true;
    }
  }
  return fully_initialized;
}

}  // namespace Helpers
}  // namespace Amanzi

