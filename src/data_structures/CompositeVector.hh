/*
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@ornl.gov)
*/

/*!

Interface for CompositeVector, an implementation of a slightly improved
Epetra_MultiVector which spans multiple simplices and knows how to
communicate itself.

CompositeVectors are a collection of vectors defined on a common mesh and
communicator.  Each vector, or component, has a name (used as a key), a mesh
Entity_kind (CELL, FACE, NODE, or BOUNDARY_FACE), and a number of degrees of
freedom (dofs).  This, along with the Epetra_BlockMap provided from the mesh
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
     local data.  This is the nasty one, because it is both subtle and
     reasonable usage.  When you access a non-const pointer, the data is
     flagged as changed.  Then you call Scatter(), and the data is flagged as
     unchanged.  Then you change the data from your old non-const pointer, and
     the data is changed, but not flagged.

     *** FIX: ALWAYS call ViewComponent() after Scatter() and before changing
         values! ***

     *** FIX2: One way to avoid protect yourself within this convention if you
         need to use the pattern "change owned values, scatter, change owned
         values again" is to put non-const references in their own scope.
         For instance, the following practice is encourage:

         // ============  begin snip  ======================================
         CompositeVector my_cv;
         ...
         { // unnamed scope for MultiVector my_vec
           Epetra_MultiVector& my_vec = *my_cv.ViewComponent("cell",false);
           my_vec[0][0] = 12;
         } // close scope of MultiVector my_vec

         my_cv.ScatterMasterToGhosted()

         // Ref to MultiVector my_vec is now gone, so we cannot use it and
         // screw things up!

         { // unnamed scope for MultiVector my_vec
           // This is now safe!
           Epetra_MultiVector& my_vec = *my_cv.ViewComponent("cell",true);
           my_vec[0][0] = my_vec[0][ghost_index] + ...
         } // close scope of MultiVector my_vec
         // ============  end snip  ========================================

     Note that any fix would require a RestoreValues() type call (PETSc's
     VecGetArray(), VecRestoreArray() make this convention explicit, but are
     not unbreakable either.)

  -- Using const_cast() and then changing the values.

     *** FIX: Const-correctness is your friend.  Keep your PKs const-correct,
         and you will never have this problem.

Note that non-INSERT modes of scatter are always done, and also always tag as
changed.  This is because subsequent calls with different modes would break
the code.

This vector provides the duck-type interface Vec and may be used with time
integrators/nonlinear solvers.


DOCUMENT VANDELAY HERE! FIX ME --etc
------------------------------------------------------------------------- */

#ifndef AMANZI_COMPOSITEVECTOR_HH_
#define AMANZI_COMPOSITEVECTOR_HH_

#define CV_ENABLE_SET_FROM_OPERATOR 0

#include <vector>
#include "Teuchos_RCP.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_CombineMode.h"
#include "Epetra_Import.h"

#include "dbc.hh"
#include "Mesh.hh"
#include "data_structures_types.hh"
#include "BlockVector.hh"
#include "CompositeVectorSpace.hh"

namespace Amanzi {

class CompositeVector {
 public:
  using VectorSpace_t = CompositeVectorSpace;
  // -- Constructors --

  // Constructor from a CompositeVectorSpace (which is like a VectorSpace).
  CompositeVector(const CompositeVectorSpace& space);
  CompositeVector(const CompositeVectorSpace& space, bool ghosted);

  // Copy constructor.
  CompositeVector(const CompositeVector& other, InitMode mode = INIT_MODE_COPY);
  CompositeVector(const CompositeVector& other, bool ghosted, InitMode mode = INIT_MODE_COPY);

  // Assignment operator.
  CompositeVector& operator=(const CompositeVector& other);

  // -- Accessors to meta-data --

  // Space/VectorSpace/Map accessor.
  const CompositeVectorSpace& Map() const { return *map_; }

  // CompositeVector maintains its own ghosted value.
  bool Ghosted() const { return ghosted_; }

  // Much accessor functionality is delegated to the VectorSpace.
  typedef std::vector<std::string>::const_iterator name_iterator;
  name_iterator begin() const { return map_->begin(); }
  name_iterator end() const { return map_->end(); }
  unsigned int size() const { return map_->size(); }

  Comm_ptr_type Comm() const { return map_->Comm(); }
  Teuchos::RCP<const AmanziMesh::Mesh> Mesh() const { return map_->Mesh(); }
  bool HasComponent(const std::string& name) const { return map_->HasComponent(name); }
  int NumComponents() const { return size(); }
  AmanziMesh::Entity_kind Location(std::string name) const { return map_->Location(name); }

  int NumVectors(const std::string& name) const { return map_->NumVectors(name); }
  int GlobalLength() const { return mastervec_->GlobalLength(); }
  long int GetLocalElementCount() const { return mastervec_->GetLocalElementCount(); }

  // Provides the size of each component's vector, either ghosted or non-ghosted.
  unsigned int size(std::string name, bool ghosted = false) const
  {
    return ghosted ? ghostvec_->size(name) : mastervec_->size(name);
  }

  // Access the VectorSpace for each component.
  Teuchos::RCP<const Epetra_BlockMap> ComponentMap(std::string name, bool ghosted = false) const
  {
    return ghosted ? ghostvec_->ComponentMap(name) : mastervec_->ComponentMap(name);
  }

  // -- View data. --

  // Access a view of a single component's data.
  //
  // Const access -- this does not tag as changed.
  Teuchos::RCP<const Epetra_MultiVector>
  ViewComponent(std::string name, bool ghosted = false) const;

  // View entries in the vectors
  //
  // Return-by-value, this does not tag as changed.
  double operator()(std::string name, int i, int j) const { return (*ghostvec_)(name, i, j); }
  double operator()(std::string name, int j) const { return (*ghostvec_)(name, 0, j); }

  // -- Set data. --

  // Access a view of a single component's data.
  //
  // Non-const access -- tags changed.
  Teuchos::RCP<Epetra_MultiVector> ViewComponent(std::string name, bool ghosted = false);

#if CV_ENABLE_SET_FROM_OPERATOR
  // Set entries in the vectors.
  //
  // Using these is VERY STRONGLY DISCOURAGED.  Instead, call ViewComponent()
  // and set entries in the MultiVector.  THESE ARE VERY VERY SLOW (But they
  // can be handy in debugging.)  Tags changed.
  double& operator()(std::string name, int i, int j)
  {
    ChangedValue(name);
    return (*ghostvec_)(name, i, j);
  }

  // Set entries in the 0th vector.
  //
  // Using these is VERY STRONGLY DISCOURAGED.  Instead, call ViewComponent()
  // and set entries in the MultiVector.  THESE ARE VERY VERY SLOW (But they
  // can be handy in debugging.)  Tags changed.
  double& operator()(std::string name, int j)
  {
    ChangedValue(name);
    return (*ghostvec_)(name, 0, j);
  }
#endif


  // Set block by pointer if possible, copy if not?????? FIX ME --etc
  void SetComponent(std::string name, const Teuchos::RCP<Epetra_MultiVector>& data);

  // -- Communicate data and communication management. --

  // Mark all components as changed.
  void ChangedValue() const;

  // Mark a single component as changed.
  void ChangedValue(std::string name) const;

  // Scatter master values to ghosted values, on all components (INSERT mode).
  //
  // Note that although scatter changes things, it doesn't change master
  // data, so we allow it to work on const.  This is necessary for a
  // non-owning PK to communicate a non-owned vector.
  //
  // Note that the scatter is ifneeded, unless force=true.
  void ScatterMasterToGhosted(bool force = false) const;

  // Scatter master values to ghosted values, on one components (INSERT mode).
  //
  // Modes shown in Epetra_CombineMode.h, but the default is Insert, which
  // overwrites the current ghost value with the (unique) new master value.
  //
  // Note that although scatter changes things, it doesn't change master
  // data, so we allow it to work on const.  This is necessary for a
  // non-owning PK to communicate a non-owned vector.
  //
  // Note that the scatter is ifneeded, unless force=true.
  void ScatterMasterToGhosted(std::string name, bool force = false) const;
  void ScatterMasterToGhosted(const char* name, bool force = false) const
  {
    ScatterMasterToGhosted(std::string(name), force);
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
  void ScatterMasterToGhosted(Epetra_CombineMode mode) const;

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
  void ScatterMasterToGhosted(std::string name, Epetra_CombineMode mode) const;

  // Combine ghosted values back to master values.
  //
  // Modes shown in Epetra_CombineMode.h, but the default is Add,
  // where off-process values are first summed into the on-process value.
  //
  // This Scatter() is not managed, and is always done.  Tags changed.
  void GatherGhostedToMaster(Epetra_CombineMode mode = Add);
  void GatherGhostedToMaster(std::string name, Epetra_CombineMode mode = Add);

  // returns non-empty importer
  const Epetra_Import& importer(std::string name) const;

  // -- Assorted vector operations, this implements a Vec --

  // Sets all vectors to value.
  int PutScalar(double scalar);

  // Sets all vectors to value including ghosted elements.
  // Different name is given so it cannot be used in a templated code.
  int PutScalarMasterAndGhosted(double scalar);

  // Sets ghost elements to value.
  // Different name is given so it cannot be used in a templated code.
  int PutScalarGhosted(double scalar);

  // v(name,:,:) = scalar
  int PutScalar(std::string name, double scalar);

  // v(name,i,:) = scalar[i]
  int PutScalar(std::string name, std::vector<double> scalar);

  // this <- scalar*this
  int Scale(double scalar);
  int ScaleMasterAndGhosted(double scalar);

  // this <- abs(this)
  int Abs(const CompositeVector& other);

  // this(name,:,:) <- scalar*this(name,:,:)
  int Scale(std::string name, double scalar);

  // this <- this + scalar
  int Shift(double scalar);

  // this(name,:,:) <- scalar + this(name,:,:)
  int Shift(std::string name, double scalar);

  // this <- element wise reciprocal(this)
  int Reciprocal(const CompositeVector& other);

  // result <- other \dot this
  int Dot(const CompositeVector& other, double* result) const;

  // this <- scalarA*A + scalarThis*this
  CompositeVector& Update(double scalarA, const CompositeVector& A, double scalarThis);

  // this <- scalarA*A + scalarB*B + scalarThis*this
  CompositeVector& Update(double scalarA,
                          const CompositeVector& A,
                          double scalarB,
                          const CompositeVector& B,
                          double scalarThis);

  // this <- scalarAB * A@B + scalarThis*this  (@ is the elementwise product
  int
  Multiply(double scalarAB, const CompositeVector& A, const CompositeVector& B, double scalarThis);

  // this <- scalarAB * A^-1@B + scalarThis*this  (@ is the elementwise product
  int ReciprocalMultiply(double scalarAB,
                         const CompositeVector& A,
                         const CompositeVector& B,
                         double scalarThis);

  // -- norms --
  int NormInf(double* norm) const;
  int Norm1(double* norm) const;
  int Norm2(double* norm) const;

  int MinValue(double* value) const;
  int MaxValue(double* value) const;
  int MeanValue(double* value) const;

  void MinValue(std::map<std::string, double>& value) const;
  void MaxValue(std::map<std::string, double>& value) const;
  void MeanValue(std::map<std::string, double>& value) const;

  // -- Utilities --

  // Write components to outstream.
  void Print(std::ostream& os, bool data_io = true) const;

  // Populate by random numbers between -1 and 1.
  int Random();

 private:
  void InitMap_(const CompositeVectorSpace& space);
  void InitData_(const CompositeVector& other, InitMode mode);
  void CreateData_();

  int Index_(std::string name) const
  {
    std::map<std::string, int>::const_iterator item = indexmap_.find(name);
    AMANZI_ASSERT(item != indexmap_.end());
    return item->second;
  }


  // The Vandelay is an Importer/Exporter which allows face unknowns
  // to be spoofed as boundary face unknowns.
  void CreateVandelayVector_() const;
  void ApplyVandelay_() const;


 private:
  Teuchos::RCP<const CompositeVectorSpace> map_;
  bool ghosted_;

  // data enumerating the blocks
  std::map<std::string, int> indexmap_;
  std::vector<std::string> names_;
  mutable std::vector<bool> ghost_are_current_;

  // data containing the blocks
  mutable Teuchos::RCP<BlockVector> ghostvec_;
  Teuchos::RCP<BlockVector> mastervec_;

  // vector for boundary data
  mutable Teuchos::RCP<Epetra_MultiVector> vandelay_vector_;
};


inline void
CompositeVector::ChangedValue() const
{
  for (std::vector<bool>::iterator lcv = ghost_are_current_.begin();
       lcv != ghost_are_current_.end();
       ++lcv)
    *lcv = false;
}

inline void
CompositeVector::ChangedValue(std::string name) const
{
  ghost_are_current_[Index_(name)] = false;
}


inline int
CompositeVector::PutScalar(double scalar)
{
  ChangedValue();
  return mastervec_->PutScalar(scalar);
}

inline int
CompositeVector::PutScalarMasterAndGhosted(double scalar)
{
  ChangedValue();
  return ghostvec_->PutScalar(scalar);
}

inline int
CompositeVector::PutScalarGhosted(double scalar)
{
  ChangedValue();
  for (int lcv_comp = 0; lcv_comp != NumComponents(); ++lcv_comp) {
    int size_owned = mastervec_->size(names_[lcv_comp]);
    int size_ghosted = ghostvec_->size(names_[lcv_comp]);

    Epetra_MultiVector& vec = *ghostvec_->ViewComponent(names_[lcv_comp]);
    for (int j = 0; j != vec.NumVectors(); ++j) {
      for (int i = size_owned; i != size_ghosted; ++i) {
        int first = vec.Map().FirstPointInElement(i);
        int ndofs = vec.Map().ElementSize(i);
        for (int k = 0; k < ndofs; ++k) { vec[j][first + k] = scalar; }
      }
    }
  }
  return 0;
}

inline int
CompositeVector::PutScalar(std::string name, double scalar)
{
  ChangedValue(name);
  return mastervec_->PutScalar(name, scalar);
}

inline int
CompositeVector::PutScalar(std::string name, std::vector<double> scalar)
{
  ChangedValue(name);
  return mastervec_->PutScalar(name, scalar);
}

inline int
CompositeVector::Abs(const CompositeVector& other)
{
  ChangedValue();
  return mastervec_->Abs(*other.mastervec_);
}

inline int
CompositeVector::Scale(double scalar)
{
  ChangedValue();
  return mastervec_->Scale(scalar);
}

inline int
CompositeVector::ScaleMasterAndGhosted(double scalar)
{
  ChangedValue();
  return ghostvec_->Scale(scalar);
}

inline int
CompositeVector::Scale(std::string name, double scalar)
{
  ChangedValue(name);
  return mastervec_->Scale(name, scalar);
}

inline int
CompositeVector::Shift(double scalar)
{
  ChangedValue();
  return mastervec_->Shift(scalar);
}

inline int
CompositeVector::Shift(std::string name, double scalar)
{
  ChangedValue(name);
  return mastervec_->Shift(name, scalar);
}

inline int
CompositeVector::Reciprocal(const CompositeVector& other)
{
  ChangedValue();
  return mastervec_->Reciprocal(*other.mastervec_);
}

inline int
CompositeVector::NormInf(double* norm) const
{
  return mastervec_->NormInf(norm);
}

inline int
CompositeVector::Norm1(double* norm) const
{
  return mastervec_->Norm1(norm);
}

inline int
CompositeVector::Norm2(double* norm) const
{
  return mastervec_->Norm2(norm);
}

inline int
CompositeVector::MinValue(double* value) const
{
  return mastervec_->MinValue(value);
}

inline int
CompositeVector::MaxValue(double* value) const
{
  return mastervec_->MaxValue(value);
}

inline int
CompositeVector::MeanValue(double* value) const
{
  return mastervec_->MeanValue(value);
}

inline void
CompositeVector::Print(std::ostream& os, bool data_io) const
{
  return mastervec_->Print(os, data_io);
}

inline int
CompositeVector::Random()
{
  ChangedValue();
  return mastervec_->Random();
}


// -----------------------------------------------------------------------------
// Non-member functions.
// -----------------------------------------------------------------------------
void
DeriveFaceValuesFromCellValues(CompositeVector&);

void
AddComponent(Teuchos::RCP<CompositeVector>& cv,
             const std::string& name,
             AmanziMesh::Entity_kind kind,
             int dim);

} // namespace Amanzi

#endif
