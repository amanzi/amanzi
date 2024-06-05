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
//#include "BoundaryFunction.hh"
#include "TreeVectorSpace.hh"
#include "TreeVector.hh"

//#include "ColumnMeshFunction.hh"
#include "CompositeVectorFunction.hh"
//#include "CompositeVectorFunctionFactory.hh"
//#include "IO.hh"

namespace Amanzi {
namespace Helpers {

// ======================================================================
// Specializations for simple data types
// ======================================================================
template <>
bool
Initialize<double>(Teuchos::ParameterList& plist, double& t)
{
  // initialize by list-provided value
  if (plist.isParameter("value")) {
    t = plist.get<double>("value");
    return true;
  }
  return false;
}


template <>
bool
Initialize<int>(Teuchos::ParameterList& plist, int& t)
{
  // initialize by list-provided value
  if (plist.isParameter("value")) {
    t = plist.get<int>("value");
    return true;
  }
  return false;
}

template <>
bool
Initialize<bool>(Teuchos::ParameterList& plist, bool& t)
{
  // initialize by list-provided value
  if (plist.isParameter("value")) {
    t = plist.get<bool>("value");
    return true;
  }
  return false;
}


// ======================================================================
// Specializations for geometric objects
// ======================================================================
template <>
bool
Initialize<AmanziGeometry::Point>(Teuchos::ParameterList& plist, AmanziGeometry::Point& p)
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
bool
Initialize<CompositeVector>(Teuchos::ParameterList& plist, CompositeVector& t)
{
  bool initialized = false;
  std::string fieldname = plist.name();

  // First try all initialization method which set the entire data structure.
  // ------ Try to set values from a restart file -----
  if (plist.isParameter("restart file")) {
    auto filename = plist.get<std::string>("restart file");
    Checkpoint chkp(filename, t.getComm(), Keys::getDomain(fieldname));
    ReadCheckpoint(chkp, plist, t);
    initialized = true;
  }

  // // ------ Try to set values from an file -----
  // if (plist.isSublist("exodus file initialization")) {
  //   // data must be pre-initialized to zero in case Exodus file does not
  //   // provide all values.
  //   t.putScalar(0.0);

  //   Teuchos::ParameterList& file_list = plist.sublist("exodus file initialization");
  //   ReadVariableFromExodusII(file_list, t);
  //   initialized = true;
  // }

  // ------ Set values using a constant -----
  if (plist.isParameter("constant")) {
    double value = plist.get<double>("constant");
    t.putScalar(value);
    initialized = true;
  }
  if (plist.isParameter("value")) {
    double value = plist.get<double>("value");
    t.putScalar(value);
    initialized = true;
  }

  // Next try all partial initialization methods -- typically cells.
  // ------ Try to set cell values from a restart file -----
  if (plist.isParameter("cells from file")) {
    AMANZI_ASSERT(false);
    auto filename = plist.get<std::string>("cells from file");
    Checkpoint chkp(filename, t.getComm());

    // read just the cells
    auto& vec_c = *t.getComponent("cell", false);
    Teuchos::ParameterList attrs(plist);
    attrs.set<AmanziMesh::Entity_kind>("location", AmanziMesh::Entity_kind::CELL);
    chkp.read(attrs, vec_c);
    initialized = true;
  }


  // ------ Set values from 1D solution -----
  // if (plist.isSublist("initialize from 1D column")) {
  //   Teuchos::ParameterList& init_plist = plist.sublist("initialize from 1D column");
  //   if (!init_plist.isParameter("f header"))
  //     init_plist.set("f header", std::string("/") + fieldname);
  //   Functions::ReadColumnMeshFunction(init_plist, t);
  //   initialized = true;
  // }

  // ------ Set values using a function -----
  if (plist.isSublist("function")) {
    double t_ini = plist.get<double>("time", 0.0);
    Teuchos::ParameterList func_plist = plist.sublist("function");

    std::vector<std::string> complist;

    // -- potential use of a mapping operator first --
    bool map_normal = plist.get<bool>("dot with normal", false);
    if (map_normal) {
      // map_normal take a vector and dots it with face normals
      AMANZI_ASSERT(t.size() == 1);                                          // one comp
      AMANZI_ASSERT(t.hasComponent("face"));                                 // is named face
      AMANZI_ASSERT(t.getLocation("face") == AmanziMesh::Entity_kind::FACE); // is on face
      AMANZI_ASSERT(t.getNumVectors("face") == 1);                           // and is scalar

      // create a vector on faces of the appropriate dimension
      int dim = t.getMesh()->getSpaceDimension();

      CompositeVectorSpace cvs;
      cvs.SetMesh(t.getMesh());
      cvs.SetComponent("face", AmanziMesh::Entity_kind::FACE, dim);
      auto vel_vec = cvs.Create();

      // Evaluate the velocity function
      auto func = Teuchos::rcp(new Functions::CompositeVectorFunction(func_plist, t.getMesh()));
      func->Compute(t_ini, *vel_vec);

      // CV's map may differ from the regular mesh map due to presense of fractures
      //Map_ptr_type fmap = t.getMap()->getMap("face", true);

      // Dot the velocity with the normal
      // -- two branches: single flux per face, multiple fluxes
      {
        auto dat_f = t.viewComponent("face", false);
        const auto vel_f = vel_vec->viewComponent("face", false);
        const AmanziMesh::Mesh& mesh = *t.getMesh();
        Kokkos::parallel_for(
          "CV::Initialize, velocity dot-with-normal",
          mesh.getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED),
          KOKKOS_LAMBDA(const int f) {
            AmanziGeometry::Point vel(dim);
            for (int i = 0; i < dim; ++i) vel[i] = vel_f(f, i);

            // int ndofs = fmap.ElementSize(f);
            // int g = fmap.FirstPointInElement(f);
            // if (ndofs == 1) {
            auto normal = mesh.getFaceNormal(f);
            dat_f(f, 0) = vel * normal;
            // } else {
            //   AmanziMesh::Entity_ID_List cells;
            //   cells = t.getMesh()->getFaceCells(f);

            //   for (int i = 0; i < ndofs; ++i) {
            //     const AmanziGeometry::Point& normal = t.getMesh()->getFaceNormal(f, cells[i], &dir);
            //     dat_f[0][g + i] = (vel * normal) * dir;
            //   }
            // }
          });
        initialized = true;
      }

    } else {
      // no map, just evaluate the function
      auto func = Teuchos::rcp(new Functions::CompositeVectorFunction(func_plist, t.getMesh()));
      func->Compute(t_ini, t);
      initialized = true;
    }
  }

  if ((t.hasComponent("face") || t.hasComponent("boundary_face")) && t.hasComponent("cell") &&
      plist.get<bool>("initialize faces from cells", false)) {
    DeriveFaceValuesFromCellValues(t);
  }
  return initialized;
}


// ======================================================================
// Specializations for Epetra_Vector
// ======================================================================
template <>
bool
Initialize<Vector_type>(Teuchos::ParameterList& plist, Vector_type& t)
{
  bool initialized = false;

  // First try all initialization method which set the entire data structure.
  // ------ Try to set values from a restart file -----
  if (plist.isParameter("restart file")) {
    auto filename = plist.get<std::string>("restart file");
    auto comm = t.getMap()->getComm();
    Checkpoint chkp(filename, comm);
    ReadCheckpoint(chkp, plist, t);
    initialized = true;
  }

  // ------ Set values using a constant -----
  if (plist.isParameter("constant")) {
    double value = plist.get<double>("constant");
    t.putScalar(value);
    initialized = true;
  }
  if (plist.isParameter("value")) {
    double value = plist.get<double>("value");
    t.putScalar(value);
    initialized = true;
  }
  return initialized;
}


// ======================================================================
// Specializations for TreeVector
// ======================================================================
template <>
bool
Initialize<TreeVector>(Teuchos::ParameterList& plist, TreeVector& t)
{
  return false;
}


template <>
void
WriteVis<TreeVector>(const Visualization& vis, Teuchos::ParameterList& attrs, const TreeVector& vec)
{
  int i = 0;
  for (const auto& subvec : vec) {
    Teuchos::ParameterList sublist(attrs);
    sublist.setName(attrs.name() + "_" + std::to_string(i));
    WriteVis(vis, sublist, *subvec);
    ++i;
  }
  if (vec.getData() != Teuchos::null) {
    Teuchos::ParameterList datalist(attrs);
    datalist.setName(attrs.name() + "_data");
    WriteVis(vis, datalist, *vec.getData());
  }
}

template <>
void
WriteCheckpoint<TreeVector>(const Checkpoint& chkp,
                            Teuchos::ParameterList& attrs,
                            const TreeVector& vec)
{
  int i = 0;
  for (const auto& subvec : vec) {
    Teuchos::ParameterList sublist(attrs);
    sublist.setName(attrs.name() + "_" + std::to_string(i));
    WriteCheckpoint(chkp, sublist, *subvec);
    ++i;
  }
  if (vec.getData() != Teuchos::null) {
    Teuchos::ParameterList datalist(attrs);
    datalist.setName(attrs.name() + "_data");
    WriteCheckpoint(chkp, datalist, *vec.getData());
  }
}

template <>
void
ReadCheckpoint<TreeVector>(const Checkpoint& chkp, Teuchos::ParameterList& attrs, TreeVector& vec)
{
  int i = 0;
  for (const auto& subvec : vec) {
    Teuchos::ParameterList sublist(attrs);
    sublist.setName(attrs.name() + "_" + std::to_string(i));
    ReadCheckpoint(chkp, sublist, *subvec);
    ++i;
  }
  if (vec.getData() != Teuchos::null) {
    Teuchos::ParameterList datalist(attrs);
    datalist.setName(attrs.name() + "_data");
    ReadCheckpoint(chkp, datalist, *vec.getData());
  }
}


// ======================================================================
// Specializations for Teuchos::Array<double>
// ======================================================================
template <>
void
WriteVis<Teuchos::Array<double>>(const Visualization& vis,
                                 Teuchos::ParameterList& attrs,
                                 const Teuchos::Array<double>& vec)
{
  auto names = OutputUtils::getNames(attrs, vec.size());
  for (int i = 0; i != vec.size(); ++i) {
    Teuchos::ParameterList attrs_i(attrs);
    attrs_i.setName(names[i]);
    WriteVis(vis, attrs_i, vec[i]);
  }
}

template <>
void
WriteCheckpoint<Teuchos::Array<double>>(const Checkpoint& chkp,
                                        Teuchos::ParameterList& attrs,
                                        const Teuchos::Array<double>& vec)
{
  auto names = OutputUtils::getNames(attrs, vec.size());
  for (int i = 0; i != vec.size(); ++i) {
    Teuchos::ParameterList attrs_i(attrs);
    attrs_i.setName(names[i]);
    WriteCheckpoint(chkp, attrs_i, vec[i]);
  }
}

template <>
void
ReadCheckpoint<Teuchos::Array<double>>(const Checkpoint& chkp,
                                       Teuchos::ParameterList& attrs,
                                       Teuchos::Array<double>& vec)
{
  auto names = OutputUtils::getNames(attrs, vec.size());
  for (int i = 0; i != vec.size(); ++i) {
    Teuchos::ParameterList attrs_i(attrs);
    attrs_i.setName(names[i]);
    ReadCheckpoint(chkp, attrs_i, vec[i]);
  }
}


template <>
bool
Initialize<Teuchos::Array<double>>(Teuchos::ParameterList& plist, Teuchos::Array<double>& t)
{
  if (plist.isParameter("value")) {
    auto arr = plist.get<Teuchos::Array<double>>("value");
    if (t.size() != arr.size()) {
      Errors::Message msg;
      msg << "While initializing \"" << plist.name() << "\" an array of lenth " << (int)t.size()
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
