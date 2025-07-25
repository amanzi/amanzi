/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------

ATS and Amanzi

Interface for CompositeVector, an implementation of a slightly improved
Epetra_MultiVector which spans multiple simplices and knows how to
communicate itself.

CompositeVectors are a collection of vectors defined on a common mesh and
communicator.  Each vector, or component, has a name (used as a key), a mesh
Entity_kind (CELL, FACE, NODE, or BOUNDARY_FACE), and a number of degrees of
freedom (dofs).  This, along with the Epetra_BlockMap provided by the mesh
on a given Entity_kind, is enough to create an Epetra_MultiVector.

Note that construction of the CompositeVector does not allocate the
Epetra_MultiVectors.  CreateData() must be called before usage.

Access using operator() is slow, and should only be used for debugging.
Prefer to use the ViewComponent() accessors.

Ghost cell updates are managed by the CompositeVector.  The design of this pattern
is prompted by two things:

  -- The need for updated ghost cell information is typically known by the
     user just prior to being used, not just after the non-ghost values are
     updated.

  -- Occasionally multiple functions need ghost values, but no changes to
     owned data have been made between these functions.  However, it is not
     always possible for the second call to know, for certain, that the first
     call did the communication.  Versatility means many code paths may be
     followed.

Therefore we use the following pattern:

  -- Each time the values of the vector are changed, flags are marked to
     record that the ghost values are stale.

  -- Each time ghost cells are needed, that flag is checked and communication
     is done, IF NEEDED.

Keeping these flags correct is therefore critical.  To do this, access to
vectors must follow some patterns; prefer to change via one of the first three
methods.  The following modifications tag the flag:

  -- Any of the usual PutScalar(), Apply(), etc methods tag changed.

  -- Non-const calls to ViewComponent() tag changed.

  -- GatherMasterToGhosted() tags changed.

  -- The ChangedValues() call manually tags changed.

  -- Scatter() called in non-INSERT mode tags changd.

Known ways to break this paradigm.  If you do any of these, it is not my fault
that you get strange parallel bugs! :

  -- Store a non-const pointer to the underlying Epetra_MultiVector.

      *** FIX: NEVER store a pointer to the underlying data, just keep
          pointers to the CompositeVector itself. ***

  -- Grab a non-const pointer, call Scatter(), then change the values of the
     local data.  This is the nasty one, because it is both subtle and common.
     When you access a non-const pointer, the data is flagged as changed.
     Then you call Scatter(), and the data is flagged as unchanged.  Then you
     change the data from your old non-const pointer, and the data is changed,
     but not flagged.

     *** FIX: ALWAYS call ViewComponent() after Scatter() and before changing
         values! ***

     Note that any fix would require a RestoreValues() type call (PETSc's
     VecGetArray(), VecRestoreArray() make this convention explicit, but are
     not unbreakable either.)

  -- Using const_cast() and then changing the values.

     *** FIX: Const-correctness is your friend.  Keep your PKs const-correct,
         and you will never have this problem.

Note that non-INSERT modes of scatter are always done, and also always tag as
changed.  This is because subsequent calls with different modes would break
the code.


DOCUMENT VANDELAY HERE! FIX ME --etc
------------------------------------------------------------------------- */

#include "Epetra_Vector.h"

#include "dbc.hh"
#include "errors.hh"
#include "CompositeVector.hh"

#ifndef MANAGED_COMMUNICATION
#define MANAGED_COMMUNICATION 0
#endif

namespace Amanzi {

// Private constructor
CompositeVector::CompositeVector(const CompositeVectorSpace& space, bool ghosted, InitMode mode)
  : map_(Teuchos::rcp(new CompositeVectorSpace(space))),
    ghosted_(ghosted),
    indexmap_(space.indexmap_),
    names_(space.names_),
    ghost_are_current_(space.NumComponents(), false)
{
  InitMap_();
  if (mode > InitMode::NOALLOC) CreateData_();
  if (mode == InitMode::ZERO) PutScalarMasterAndGhosted(0.);
}


void
CompositeVector::InitMap_()
{
  // generate the master's maps
  std::vector<Teuchos::RCP<const Epetra_BlockMap>> mastermaps;
  for (const auto& name : *map_) {
    mastermaps.emplace_back(map_->Map(name, false));
  }

  // create the master BlockVector
  mastervec_ = Teuchos::rcp(new BlockVector(Comm(), names_, mastermaps, map_->num_dofs_));

  // do the same for the ghosted Vector, if necessary
  if (ghosted_) {
    // generate the ghost's maps
    std::vector<Teuchos::RCP<const Epetra_BlockMap>> ghostmaps;
    for (const auto& name : *map_) {
      ghostmaps.emplace_back(map_->Map(name, true));
    }
    // create the ghost BlockVector
    ghostvec_ = Teuchos::rcp(new BlockVector(Comm(), names_, ghostmaps, map_->num_dofs_));
  } else {
    ghostvec_ = mastervec_;
  }
};


// Sets sizes of vectors, instantiates Epetra_Vectors, and preps for lazy
// creation of everything else.
void
CompositeVector::CreateData_()
{
  if (!Mesh().get()) {
    Errors::Message message("CompositeVector: construction called with no mesh.");
    Exceptions::amanzi_throw(message);
  }
  if (NumComponents() == 0) {
    Errors::Message message("CompositeVector: construction called with no components.");
    Exceptions::amanzi_throw(message);
  }

  ghost_are_current_.resize(NumComponents(), false);

  // Create the ghost vector.  Note this is also the master vector if not ghosted.
  ghostvec_->CreateData();

  // If the vector is ghosted, create the master from views of the ghost.
  if (ghosted_) {
    for (name_iterator name = begin(); name != end(); ++name) {
      // get the ghost component's data
      Teuchos::RCP<Epetra_MultiVector> g_comp = ghostvec_->ViewComponent(*name);
      double** data;
      g_comp->ExtractView(&data);

      // create the master component
      Teuchos::RCP<Epetra_MultiVector> m_comp = Teuchos::rcp(new Epetra_MultiVector(
        View, *mastervec_->ComponentMap(*name), data, mastervec_->NumVectors(*name)));

      // push it back into the master vec
      mastervec_->SetComponent(*name, m_comp);
    }
  }
  PutScalarMasterAndGhosted(0.);
};


CompositeVector&
CompositeVector::operator=(const CompositeVector& other)
{
  if (this != &other) {
    AMANZI_ASSERT(Map().SubsetOf(other.Map()));

    if (Ghosted() && other.Ghosted()) {
      // If both are ghosted, copy the ghosted vector.
      // Exception: boundary_face component created on a fly is not ghosted
      for (name_iterator name = begin(); name != end(); ++name) {
        bool ghosted =
          (*name == "boundary_face" && other.HasImportedComponent("face")) ? false : true;
        Teuchos::RCP<Epetra_MultiVector> comp = ViewComponent(*name, ghosted);
        Teuchos::RCP<const Epetra_MultiVector> othercomp = other.ViewComponent(*name, ghosted);
        *comp = *othercomp;
      }

    } else {
      // Copy the non-ghosted data.  NOTE: any ghosted data is undefined!
      for (name_iterator name = begin(); name != end(); ++name) {
        Teuchos::RCP<Epetra_MultiVector> comp = ViewComponent(*name, false);
        Teuchos::RCP<const Epetra_MultiVector> othercomp = other.ViewComponent(*name, false);
        *comp = *othercomp;
      }
    }
  }

  ChangedValue();
  return *this;
};


// view data
// -- Access a view of a single component's data.
// Ghosted views are simply the vector itself, while non-ghosted views are
// lazily generated.
Teuchos::RCP<const Epetra_MultiVector>
CompositeVector::ViewComponent(const std::string& name, bool ghosted) const
{
  AMANZI_ASSERT(HasImportedComponent(name));
  if (name == "boundary_face") {
    if (!HasComponent("boundary_face") && HasComponent("face")) {
      ApplyVandelay_();
      return ghosted ? vandelay_vector_all_ : vandelay_vector_owned_;
    }
  }

  if (ghosted) {
    return ghostvec_->ViewComponent(name);
  } else {
    return mastervec_->ViewComponent(name);
  }
};


Teuchos::RCP<Epetra_MultiVector>
CompositeVector::ViewComponent(const std::string& name, bool ghosted)
{
  AMANZI_ASSERT(HasComponent(name));
  ChangedValue(name);
  if (ghosted) {
    return ghostvec_->ViewComponent(name);
  } else {
    return mastervec_->ViewComponent(name);
  }
};


// View a vector slice of a single DoF
Teuchos::RCP<const CompositeVector>
CompositeVector::GetVector(size_t i) const
{
  std::vector<int> num_dofs_i(NumComponents(), 1);

  std::vector<AmanziMesh::Entity_kind> locations;
  int check_same_ndofs = -1;
  for (const auto& name : names_) {
    if (check_same_ndofs < 0) check_same_ndofs = NumVectors(name);

    if (check_same_ndofs != NumVectors(name)) {
      Errors::Message msg(
        "Cannot slice a vector with differing numbers of DoFs across components.");
      Exceptions::amanzi_throw(msg);
    }
    locations.emplace_back(Location(name));
  }

  // create the space
  CompositeVectorSpace cvs;
  cvs.SetMesh(map_->Mesh())->SetGhosted(ghosted_)->SetComponents(names_, locations, num_dofs_i);

  // create the vector
  auto cv = Teuchos::rcp(new CompositeVector(cvs, ghosted_, InitMode::NOALLOC));

  // for each component, set the vector as the MultiVector's dof vector,
  // through a NON-OWNING (these are owned by this's MultiVector) RCP to the
  // Epetra_Vector.
  for (const auto& name : names_) {
    // Have to const-cast to construct the new vector, but will return the new
    // vector as const.  This is safe, despite appearances.
    Epetra_MultiVector const* comp_vec = (*ViewComponent(name, true))(i);
    Teuchos::RCP<Epetra_MultiVector> comp_vec_ptr =
      Teuchos::rcpFromRef<Epetra_MultiVector>(*const_cast<Epetra_MultiVector*>(comp_vec));
    cv->SetComponent(name, comp_vec_ptr);
  }
  return cv;
}


// View a vector slice of a single DoF
Teuchos::RCP<CompositeVector>
CompositeVector::GetVector(size_t i)
{
  std::vector<int> num_dofs_i(NumComponents(), 1);

  std::vector<AmanziMesh::Entity_kind> locations;
  int check_same_ndofs = -1;
  for (const auto& name : names_) {
    if (check_same_ndofs < 0) check_same_ndofs = NumVectors(name);

    if (check_same_ndofs != NumVectors(name)) {
      Errors::Message msg(
        "Cannot slice a vector with differing numbers of DoFs across components.");
      Exceptions::amanzi_throw(msg);
    }
    locations.emplace_back(Location(name));
  }

  // create the space
  CompositeVectorSpace cvs;
  cvs.SetMesh(map_->Mesh())->SetGhosted(ghosted_)->SetComponents(names_, locations, num_dofs_i);

  // create the vector
  auto cv = Teuchos::rcp(new CompositeVector(cvs, ghosted_, InitMode::NOALLOC));

  // for each component, set the vector as the MultiVector's dof vector,
  // through a NON-OWNING (these are owned by this's MultiVector) RCP to the
  // Epetra_Vector.
  for (const auto& name : names_) {
    cv->SetComponent(name, Teuchos::rcpFromRef(*(*ViewComponent(name, true))(i)));
  }
  return cv;
}


// Set data by pointer if possible, otherwise by copy.
void
CompositeVector::SetComponent(const std::string& name, const Teuchos::RCP<Epetra_MultiVector>& data)
{
  ChangedValue(name);

  if (ghostvec_->ComponentMap(name)->SameAs(data->Map())) {
    // setting the ghost vec -- drop in the data in the ghost
    ghostvec_->SetComponent(name, data);

    // and create a new view for the master
    double** vals;
    data->ExtractView(&vals);
    Teuchos::RCP<Epetra_MultiVector> m_comp = Teuchos::rcp(new Epetra_MultiVector(
      View, *mastervec_->ComponentMap(name), vals, mastervec_->NumVectors(name)));
    mastervec_->SetComponent(name, m_comp);

  } else if (mastervec_->ComponentMap(name)->SameAs(data->Map())) {
    *mastervec_->ViewComponent(name) = *data;
  } else {
    Errors::Message message("Attempted set of non-compatible Component.");
    Exceptions::amanzi_throw(message);
  }
};


// -- Scatter master values to ghosted values.
// Modes shown in Epetra_CombineMode.h, but the default is Insert, which
// overwrites the current ghost value with the (unique) new master value.
void
CompositeVector::ScatterMasterToGhosted(bool force) const
{
  for (name_iterator name = begin(); name != end(); ++name) {
    ScatterMasterToGhosted(*name, force);
  }
};


void
CompositeVector::ScatterMasterToGhosted(const std::string& name, bool force) const
{
  // NOTE: allowing const is a hack to allow non-owning PKs to nonetheless
  // update ghost cells, which may be necessary for their discretization
  AMANZI_ASSERT(ghosted_);
#ifdef HAVE_MPI
#if MANAGED_COMMUNICATION
  if (ghosted_ && ((!ghost_are_current_[Index_(name)]) || force)) {
#else
  if (ghosted_) {
#endif
    // communicate
    Teuchos::RCP<Epetra_MultiVector> g_comp = ghostvec_->ViewComponent(name);
    Teuchos::RCP<const Epetra_MultiVector> m_comp = mastervec_->ViewComponent(name);
    g_comp->Import(*m_comp, importer(name), Insert);

#if MANAGED_COMMUNICATION
    // mark as communicated
    ghost_are_current_[Index_(name)] = true;
#endif
  }
#endif
};


// Scatter master values to ghosted values, on all components, in a mode.
//
// Modes shown in Epetra_CombineMode.h, but the default is Insert, which
// overwrites the current ghost value with the (unique) new master value.
//
// Note that although scatter changes things, it doesn't change master
// data, so we allow it to work on const.  This is necessary for a
// non-owning PK to communicate a non-owned vector.
//
// This Scatter() is not managed, and is always done.  Tags changed.
void
CompositeVector::ScatterMasterToGhosted(Epetra_CombineMode mode) const
{
  for (name_iterator name = begin(); name != end(); ++name) {
    ScatterMasterToGhosted(*name, mode);
  }
}

// Scatter master values to ghosted values, on all components, in a mode.
//
// Modes shown in Epetra_CombineMode.h, but the default is Insert, which
// overwrites the current ghost value with the (unique) new master value.
//
// Note that although scatter changes things, it doesn't change master
// data, so we allow it to work on const.  This is necessary for a
// non-owning PK to communicate a non-owned vector.
//
// This Scatter() is not managed, and is always done.  Tags changed.
void
CompositeVector::ScatterMasterToGhosted(const std::string& name, Epetra_CombineMode mode) const
{
  ChangedValue(name);

#ifdef HAVE_MPI
  if (ghosted_) {
    // communicate
    Teuchos::RCP<Epetra_MultiVector> g_comp = ghostvec_->ViewComponent(name);
    Teuchos::RCP<const Epetra_MultiVector> m_comp = mastervec_->ViewComponent(name);
    g_comp->Import(*m_comp, importer(name), Insert);

    // mark as communicated
    ghost_are_current_[Index_(name)] = true;
  }
#endif
};


// -- Combine ghosted values back to master values.
// Modes shown in Epetra_CombineMode.h, but the default is Add,
// where off-process values are first summed into the on-process value.
void
CompositeVector::GatherGhostedToMaster(Epetra_CombineMode mode)
{
  for (name_iterator name = begin(); name != end(); ++name) {
    GatherGhostedToMaster(*name, mode);
  }
};


void
CompositeVector::GatherGhostedToMaster(const std::string& name, Epetra_CombineMode mode)
{
  ChangedValue(name);
#ifdef HAVE_MPI
  if (ghosted_) {
    // communicate
    Teuchos::RCP<const Epetra_MultiVector> g_comp = ghostvec_->ViewComponent(name);
    Teuchos::RCP<Epetra_MultiVector> m_comp = mastervec_->ViewComponent(name);
    m_comp->Export(*g_comp, importer(name), mode);
  }
#endif
};


void
CompositeVector::CreateVandelayVector_() const
{
  vandelay_vector_all_ = Teuchos::rcp(
    new Epetra_MultiVector(Mesh()->getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE, true),
                           mastervec_->NumVectors("face"),
                           false));
  if (ghosted_) {
    double** data;
    vandelay_vector_all_->ExtractView(&data);
    vandelay_vector_owned_ = Teuchos::rcp(
      new Epetra_MultiVector(View,
                             Mesh()->getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE, false),
                             data,
                             mastervec_->NumVectors("face")));
  }

  // create new importer from face-component if it has more than one DOF per face
  //
  // Note that vandelay operates on OWNED entities, not ALL
  const auto& map = *ComponentMap("face", false);
  int nfaces = Mesh()->getMap(AmanziMesh::Entity_kind::FACE, false).NumMyElements();
  if (nfaces != map.NumMyPoints()) {
    vandelay_import_ = Teuchos::rcp(new Epetra_Import(Mesh()->getBoundaryFaceImporter()));

    auto data = vandelay_import_->PermuteFromLIDs();
    int n = vandelay_import_->NumPermuteIDs();

    for (int i = 0; i < n; ++i) {
      data[i] = map.FirstPointInElement(data[i]);
    }
  }
}

void
CompositeVector::ApplyVandelay_() const
{
  if (vandelay_vector_owned_ == Teuchos::null) {
    CreateVandelayVector_();
  }
  if (vandelay_import_ == Teuchos::null)
    vandelay_vector_owned_->Import(
      *ViewComponent("face", false), Mesh()->getBoundaryFaceImporter(), Insert);
  else vandelay_vector_owned_->Import(*ViewComponent("face", false), *vandelay_import_, Insert);
}


// return non-empty importer
const Epetra_Import&
CompositeVector::importer(const std::string& name) const
{
  return Mesh()->getImporter(Location(name));
}


// Mathematical operations
// -- result <- other \dot this
int
CompositeVector::Dot(const CompositeVector& other, double* result) const
{
  double tmp_result = 0.0;
  for (name_iterator lcv = begin(); lcv != end(); ++lcv) {
    if (other.HasImportedComponent(*lcv)) {
      std::vector<double> intermediate_result(ViewComponent(*lcv, false)->NumVectors(), 0.0);
      int ierr =
        ViewComponent(*lcv, false)->Dot(*other.ViewComponent(*lcv, false), &intermediate_result[0]);
      if (ierr) return ierr;

      for (int lcv_vector = 0; lcv_vector != NumVectors(*lcv); ++lcv_vector) {
        tmp_result += intermediate_result[lcv_vector];
      }
    }
  }
  *result = tmp_result;
  return 0;
};


// -- this <- scalarA*A + scalarThis*this
CompositeVector&
CompositeVector::Update(double scalarA, const CompositeVector& A, double scalarThis)
{
  //  AMANZI_ASSERT(map_->SubsetOf(*A.map_));
  ChangedValue();
  for (name_iterator lcv = begin(); lcv != end(); ++lcv) {
    if (A.HasImportedComponent(*lcv))
      ViewComponent(*lcv, false)->Update(scalarA, *A.ViewComponent(*lcv, false), scalarThis);
  }
  return *this;
};


// -- this <- scalarA*A + scalarB*B + scalarThis*this
CompositeVector&
CompositeVector::Update(double scalarA,
                        const CompositeVector& A,
                        double scalarB,
                        const CompositeVector& B,
                        double scalarThis)
{
  ChangedValue();
  for (name_iterator lcv = begin(); lcv != end(); ++lcv) {
    if (A.HasImportedComponent(*lcv) && B.HasImportedComponent(*lcv))
      ViewComponent(*lcv, false)
        ->Update(scalarA,
                 *A.ViewComponent(*lcv, false),
                 scalarB,
                 *B.ViewComponent(*lcv, false),
                 scalarThis);
  }
  return *this;
};


// -- this <- scalarAB * A@B + scalarThis*this  (@ is the elementwise product
int
CompositeVector::Multiply(double scalarAB,
                          const CompositeVector& A,
                          const CompositeVector& B,
                          double scalarThis)
{
  //  AMANZI_ASSERT(map_->SubsetOf(*A.map_));
  //  AMANZI_ASSERT(map_->SubsetOf(*B.map_));
  ChangedValue();
  int ierr = 0;
  for (name_iterator lcv = begin(); lcv != end(); ++lcv) {
    if (A.HasImportedComponent(*lcv) && B.HasImportedComponent(*lcv))
      ierr |=
        ViewComponent(*lcv, false)
          ->Multiply(
            scalarAB, *A.ViewComponent(*lcv, false), *B.ViewComponent(*lcv, false), scalarThis);
  }
  return ierr;
};


// -- this <- scalarAB * B / A + scalarThis*this  (/ is the elementwise division
int
CompositeVector::ReciprocalMultiply(double scalarAB,
                                    const CompositeVector& A,
                                    const CompositeVector& B,
                                    double scalarThis)
{
  AMANZI_ASSERT(map_->SubsetOf(*A.map_));
  AMANZI_ASSERT(map_->SubsetOf(*B.map_));
  ChangedValue();
  int ierr = 0;
  for (name_iterator lcv = begin(); lcv != end(); ++lcv) {
    if (A.HasImportedComponent(*lcv) && B.HasImportedComponent(*lcv))
      ierr |=
        ViewComponent(*lcv, false)
          ->ReciprocalMultiply(
            scalarAB, *A.ViewComponent(*lcv, false), *B.ViewComponent(*lcv, false), scalarThis);
  }
  return ierr;
};


// Mathematical operations
// -- minimum value by component
void
CompositeVector::MinValue(std::map<std::string, double>& value) const
{
  value.clear();

  for (int n = 0; n != names_.size(); ++n) {
    double tmp(1e+50), value_loc[1];
    const Epetra_MultiVector& comp = *ViewComponent(names_[n]);

    for (int i = 0; i != comp.NumVectors(); ++i) {
      comp(i)->MinValue(value_loc);
      tmp = std::min(tmp, value_loc[0]);
    }
    value[names_[n]] = tmp;
  }
};


// -- maximum value by component
void
CompositeVector::MaxValue(std::map<std::string, double>& value) const
{
  value.clear();

  for (int n = 0; n != names_.size(); ++n) {
    double tmp(-1e+50), value_loc[1];
    const Epetra_MultiVector& comp = *ViewComponent(names_[n]);

    for (int i = 0; i != comp.NumVectors(); ++i) {
      comp(i)->MaxValue(value_loc);
      tmp = std::max(tmp, value_loc[0]);
    }
    value[names_[n]] = tmp;
  }
};


// -- mean value by component
void
CompositeVector::MeanValue(std::map<std::string, double>& value) const
{
  value.clear();

  for (int n = 0; n != names_.size(); ++n) {
    int ni, nt(0);
    double tmp(0.0), value_loc[1];
    const Epetra_MultiVector& comp = *ViewComponent(names_[n]);

    for (int i = 0; i != comp.NumVectors(); ++i) {
      ni = comp(i)->GlobalLength();
      comp(i)->MeanValue(value_loc);
      tmp += value_loc[0] * ni;
      nt += ni;
    }
    value[names_[n]] = tmp / nt;
  }
};


// -----------------------------------------------------------------------------
// Interpolate pressure ICs on cells to ICs for lambda (faces).
// -----------------------------------------------------------------------------
void
DeriveFaceValuesFromCellValues(CompositeVector& cv)
{
  if (cv.HasComponent("face")) {
    cv.ScatterMasterToGhosted("cell");
    const Epetra_MultiVector& cv_c = *cv.ViewComponent("cell", true);
    Epetra_MultiVector& cv_f = *cv.ViewComponent("face", false);

    int f_owned = cv_f.MyLength();
    for (int f = 0; f != f_owned; ++f) {
      auto cells = cv.Mesh()->getFaceCells(f);
      int ncells = cells.size();

      double face_value = 0.0;
      for (int n = 0; n != ncells; ++n) {
        face_value += cv_c[0][cells[n]];
      }
      cv_f[0][f] = face_value / ncells;
    }
  } else if (cv.HasComponent("boundary_face")) {
    const Epetra_MultiVector& cv_c = *cv.ViewComponent("cell", true);
    Epetra_MultiVector& cv_f = *cv.ViewComponent("boundary_face", false);

    const Epetra_BlockMap& fb_map =
      cv.Mesh()->getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE, false);
    const Epetra_BlockMap& f_map = cv.Mesh()->getMap(AmanziMesh::Entity_kind::FACE, false);

    int fb_owned = cv_f.MyLength();
    for (int fb = 0; fb != fb_owned; ++fb) {
      int f_gid = fb_map.GID(fb);
      int f_lid = f_map.LID(f_gid);

      auto cells = cv.Mesh()->getFaceCells(f_lid);
      int ncells = cells.size();

      AMANZI_ASSERT((ncells == 1));

      double face_value = cv_c[0][cells[0]];
      cv_f[0][fb] = face_value;
    }
  }
}


// -----------------------------------------------------------------------------
// Non-member function: extension of a vector via data copy.
// -----------------------------------------------------------------------------
void
AddComponent(Teuchos::RCP<CompositeVector>& cv,
             const std::string& name,
             AmanziMesh::Entity_kind kind,
             int dim)
{
  // copy construct the CVS making it not owned and add the new component
  CompositeVectorSpace new_space(cv->Map());
  new_space.SetOwned(false);
  new_space.AddComponent(name, kind, dim);

  // create the new vector and copy data
  Teuchos::RCP<CompositeVector> new_cv = Teuchos::rcp(new CompositeVector(new_space));
  bool ghost = new_space.Ghosted();
  for (auto it = cv->Map().begin(); it != cv->Map().end(); ++it) {
    auto data1 = cv->ViewComponent(*it, ghost);
    auto data2 = new_cv->ViewComponent(*it, ghost);
    *data2 = *data1;
  }

  // replace the vector
  cv = new_cv;
}

} // namespace Amanzi
