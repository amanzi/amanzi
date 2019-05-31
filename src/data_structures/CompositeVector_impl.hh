/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*

Implementation for CompositeVector, an implementation of a block Tpetra
vector which spans multiple simplices and knows how to communicate
itself.

CompositeVectors are a collection of vectors defined on a common mesh and
communicator.  Each vector, or component, has a name (used as a key), a mesh
Entity_kind (CELL, FACE, NODE, or BOUNDARY_FACE), and a number of degrees of
freedom (dofs).  This, along with the Map_type provided from the mesh on a
given Entity_kind, is enough to create a Vector

Note that construction of the CompositeVector does not allocate the
Tpetra_Vector.  CreateData() must be called before usage.

Access using operator() is slow, and should only be used for debugging.
Prefer to use the ViewComponent() accessors.

DOCUMENT VANDELAY HERE! FIX ME --etc

*/

#ifndef AMANZI_COMPOSITEVECTOR_IMPL_HH_
#define AMANZI_COMPOSITEVECTOR_IMPL_HH_


#include "CompositeMap.hh"

namespace Amanzi {


// Constructor
template<typename Scalar>
CompositeVector_<Scalar>::CompositeVector_(const Teuchos::RCP<const CompositeMap>& map)
    : map_(map),
      mastervec_(Teuchos::rcp(new CompVector<Scalar>(map->Map(false)))),
      ghostvec_(Teuchos::rcp(new CompVector<Scalar>(map->Map(true)))),
      ghosted_(true)
{
  CreateData_();
}


template<typename Scalar>
CompositeVector_<Scalar>::CompositeVector_(const Teuchos::RCP<const CompositeMap>& map, bool ghosted)
    : map_(map),
      mastervec_(Teuchos::rcp(new CompVector<Scalar>(map->Map(false)))),
      ghostvec_(Teuchos::rcp(new CompVector<Scalar>(map->Map(true)))),
      ghosted_(ghosted)
{
  CreateData_();
}


template<typename Scalar>
CompositeVector_<Scalar>::CompositeVector_(const CompositeVector_<Scalar>& other, InitMode mode)
    : map_(other.Map()),
      mastervec_(Teuchos::rcp(new CompVector<Scalar>(other.Map()->Map(false)))),
      ghostvec_(Teuchos::rcp(new CompVector<Scalar>(other.Map()->Map(true)))),
      ghosted_(other.ghosted_)
{
  CreateData_();
  InitData_(other, mode);
}

template<typename Scalar>
CompositeVector_<Scalar>::CompositeVector_(const CompositeVector_<Scalar>& other, bool ghosted, InitMode mode)
    : map_(other->Map()),
      mastervec_(Teuchos::rcp(new CompVector<Scalar>(other.Map()->Map(false)))),
      ghostvec_(Teuchos::rcp(new CompVector<Scalar>(other.Map()->Map(true)))),
      ghosted_(ghosted)
{
  CreateData_();
  InitData_(other, mode);
}


// Initialize data
template<typename Scalar>
void
CompositeVector_<Scalar>::InitData_(const CompositeVector_<Scalar>& other, InitMode mode)
{
  // Trilinos inits to 0
  //  if (mode == INIT_MODE_ZERO) {
  //    PutScalar(0.);
  if (mode == INIT_MODE_COPY) {
    *this = other;
  }
}

// Sets sizes of vectors, instantiates Vectors, and preps for lazy
// creation of everything else.
template<typename Scalar>
void
CompositeVector_<Scalar>::CreateData_() {
  if (ghosted_) {
    // Create the ghost vector.  Note this is also the master vector if not ghosted.
    ghostvec_->CreateData();

    // If the vector is ghosted, create the master from views of the ghost.
    for (const auto& name : *this) {
      auto g_comp = ghostvec_->GetComponent(name);
      auto m_map = mastervec_->ComponentMap(name);
      if (m_map != Teuchos::null) {
        // get the ghost component's data
        auto m_comp = g_comp->offsetViewNonConst(m_map, 0);
        mastervec_->SetComponent(name, m_comp);
      }
    }
  } else {
    mastervec_->CreateData();
    ghostvec_ = mastervec_;
  }  
  PutScalarMasterAndGhosted(0.);
};


template<typename Scalar>
CompositeVector_<Scalar>& 
CompositeVector_<Scalar>::operator=(const CompositeVector_<Scalar>& other)
{
  if (this != &other) {
    if (ghosted_ && other.ghosted_) {
      // If both are ghosted, DEEP copy the ghosted vector.
      for (const auto& name : *this) {
        auto& comp = GetComponent_(name, true);
        const auto& othercomp = other.GetComponent_(name, true);
        comp.assign(othercomp); // note operator= is shallow copy, this is deep
      }

    } else {
      // DEEP Copy the non-ghosted data.  NOTE: any ghosted data is undefined!
      for (const auto& name : *this) {
        auto& comp = GetComponent_(name, false);
        const auto& othercomp = other.GetComponent_(name, false);
        comp.assign(othercomp); // note operator= is shallow copy, this is deep
      }
    }
  }
  return *this;
};



template<typename Scalar>
typename CompositeVector_<Scalar>::name_iterator
CompositeVector_<Scalar>::begin() const { return Map()->begin(); }

template<typename Scalar>
typename CompositeVector_<Scalar>::name_iterator
CompositeVector_<Scalar>::end() const { return Map()->end(); }

template<typename Scalar>
std::size_t
CompositeVector_<Scalar>::size() const { return Map()->size(); }

template<typename Scalar>
BlockMap_ptr_type
CompositeVector_<Scalar>::ComponentMap(const std::string& name, bool ghosted) const { return Map()->ComponentMap(name, ghosted); }


template<typename Scalar>
Comm_ptr_type
CompositeVector_<Scalar>::Comm() const { return Map()->Comm(); }

template<typename Scalar>
bool
CompositeVector_<Scalar>::HasComponent(const std::string& name) const
{
  return Map()->HasComponent(name);
}

template<typename Scalar>
std::size_t
CompositeVector_<Scalar>::NumVectors(const std::string& name) const
{
  return Map()->NumVectors(name);
}


template<typename Scalar>
LO
CompositeVector_<Scalar>::MyLength(const std::string& name, bool ghosted) const
{
  return ghosted ? ghostvec_->MyLength(name) : mastervec_->MyLength(name);
}


template<typename Scalar>
GO
CompositeVector_<Scalar>::GlobalLength(bool ghosted) const
{
  return Map()->GlobalLength(ghosted);
}

template<typename Scalar>
int
CompositeVector_<Scalar>::PutScalar(Scalar scalar) {
  return mastervec_->PutScalar(scalar);
}

template<typename Scalar>
int
CompositeVector_<Scalar>::PutScalarMasterAndGhosted(Scalar scalar) {
  return ghostvec_->PutScalar(scalar);
}

// template<typename Scalar>
// int
// CompositeVector_<Scalar>::PutScalarGhosted(Scalar scalar) {
//   for (int lcv_comp = 0; lcv_comp != NumComponents(); ++lcv_comp) {
//     int size_owned = mastervec_->size(names_[lcv_comp]);
//     int size_ghosted = ghostvec_->size(names_[lcv_comp]);

//     using Range_type = Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, int>;
//     /*
//     auto vec = ViewComponent(names_[lcv_comp], true);
//     Kokkos::parallel_for("PutScalarGhosted", Range_type(size_owned, size_ghosted), 
// 			 [=] (const int& i) { vec(i) = scalar; })
//     */
//     auto vec = ViewComponent(names_[lcv_comp], true);
//     auto vec_ghost_view = Kokkos::subview(vec, std::make_pair(size_owned, size_ghosted), Kokkos::ALL());
//     Kokkos::parallel_for("PutScalarGhosted", Range_type(0, size_ghosted - size_owned), 
    

//   }
//   return 0;
// }

template<typename Scalar>
int
CompositeVector_<Scalar>::PutScalar(const std::string& name, Scalar scalar) {
  return mastervec_->PutScalar(name, scalar);
}

// template<typename Scalar>
// int
// CompositeVector_<Scalar>::PutScalar(const std::string& name, std::vector<Scalar> scalar) {
//   return mastervec_->PutScalar(name, scalar);
// }

template<typename Scalar>
int
CompositeVector_<Scalar>::Abs(const CompositeVector_<Scalar>& other) {
  return mastervec_->Abs(*other.mastervec_);
}

template<typename Scalar>
int
CompositeVector_<Scalar>::Scale(Scalar scalar) {
  return mastervec_->Scale(scalar);
}

template<typename Scalar>
int
CompositeVector_<Scalar>::ScaleMasterAndGhosted(Scalar scalar) {
  return ghostvec_->Scale(scalar);
}

template<typename Scalar>
int
CompositeVector_<Scalar>::Scale(const std::string& name, Scalar scalar) {
  return mastervec_->Scale(name, scalar);
}

// template<typename Scalar>
// int
// CompositeVector_<Scalar>::Shift(Scalar scalar) {
//   return mastervec_->Shift(scalar);
// }

// template<typename Scalar>
// int
// CompositeVector_<Scalar>::Shift(const std::string& name, Scalar scalar) {
//   return mastervec_->Shift(name, scalar);
// }

template<typename Scalar>
int
CompositeVector_<Scalar>::Reciprocal(const CompositeVector_<Scalar>& other) {
  return mastervec_->Reciprocal(*other.mastervec_);
}

template<typename Scalar>
int
CompositeVector_<Scalar>::NormInf(Scalar* norm) const {
  return mastervec_->NormInf(norm);
}

template<typename Scalar>
int
CompositeVector_<Scalar>::Norm1(Scalar* norm) const {
  return mastervec_->Norm1(norm);
}

template<typename Scalar>
int
CompositeVector_<Scalar>::Norm2(Scalar* norm) const {
  return mastervec_->Norm2(norm);
}

// template<typename Scalar>
// int
// CompositeVector_<Scalar>::MinValue(Scalar* value) const {
//   return mastervec_->MinValue(value);
// }

// template<typename Scalar>
// int
// CompositeVector_<Scalar>::MaxValue(Scalar* value) const {
//   return mastervec_->MaxValue(value);
// }

// template<typename Scalar>
// int
// CompositeVector_<Scalar>::MeanValue(Scalar* value) const {
//   return mastervec_->MeanValue(value);
// }

template<typename Scalar>
void
CompositeVector_<Scalar>::Print(std::ostream &os, bool data_io) const {
  return mastervec_->Print(os, data_io);
}

template<typename Scalar>
int
CompositeVector_<Scalar>::Random() {
  return mastervec_->Random();
}


// -----------------------------------------------------------------------------
// Non-member functions.
// -----------------------------------------------------------------------------
//void DeriveFaceValuesFromCellValues(CompositeVector&);

// void AddComponent(Teuchos::RCP<CompositeVector> cv,
//                   const std::string& name, AmanziMesh::Entity_kind kind, int dim);



// view data
// -- Access a view of a single component's data.
// Ghosted views are simply the vector itself, while non-ghosted views are
// lazily generated.
template<typename Scalar>
const MultiVector_type_<Scalar>&
CompositeVector_<Scalar>::GetComponent_(const std::string& name, bool ghosted) const {
  if (name == std::string("boundary_face")) {
    if (!mastervec_->HasComponent("boundary_face") &&
        mastervec_->HasComponent("face")) {
      ApplyVandelay_();
      return *vandelay_vector_;
    }
  }
  if (ghosted) {
    return *ghostvec_->GetComponent(name);
  } else {
    return *mastervec_->GetComponent(name);
  }
};


template<typename Scalar>
MultiVector_type_<Scalar>&
CompositeVector_<Scalar>::GetComponent_(const std::string& name, bool ghosted) {
  if (name == std::string("boundary_face")) {
    if (!mastervec_->HasComponent("boundary_face") &&
        mastervec_->HasComponent("face")) {
      ApplyVandelay_();
      return *vandelay_vector_;
    }
  }
  if (ghosted) {
    return *ghostvec_->GetComponent(name);
  } else {
    return *mastervec_->GetComponent(name);
  }
};


// -- Scatter master values to ghosted values.
// Modes shown in Epetra_CombineMode.h, but the default is Insert, which
// overwrites the current ghost value with the (unique) new master value.
template<typename Scalar>
void
CompositeVector_<Scalar>::ScatterMasterToGhosted(Tpetra::CombineMode mode) const
{
  for (const auto& name : *this) {
    ScatterMasterToGhosted(name, mode);
  }
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
template<typename Scalar>
void
CompositeVector_<Scalar>::ScatterMasterToGhosted(const std::string& name, Tpetra::CombineMode mode) const
{
#ifdef HAVE_MPI
  if (ghosted_) {
    // check for and create the importer if needed
    auto import = Map()->Importer(name);
    
    // communicate
    auto g_comp = ghostvec_->GetComponent(name);
    auto m_comp = mastervec_->GetComponent(name);
    g_comp->doImport(*m_comp, *import, mode);
  }
#endif
};




// -- Combine ghosted values back to master values.
// Modes shown in Epetra_CombineMode.h, but the default is Add,
// where off-process values are first summed into the on-process value.
template<typename Scalar>
void
CompositeVector_<Scalar>::GatherGhostedToMaster(Tpetra::CombineMode mode)
{
  for (const auto& name=begin(); name!=end(); ++name) {
    GatherGhostedToMaster(*name, mode);
  }
};


template<typename Scalar>
void
CompositeVector_<Scalar>::GatherGhostedToMaster(const std::string& name, Tpetra::CombineMode mode)
{
#ifdef HAVE_MPI
  if (ghosted_) {
    auto import = Map()->Importer(name);
    // communicate
    auto g_comp = ghostvec_->GetComponent(name);
    auto m_comp = mastervec_->GetComponent(name);
    m_comp->doExport(*g_comp, *import, mode);
  }
#endif
};


// Mathematical operations
// -- result <- other \dot this
template<typename Scalar>
int
CompositeVector_<Scalar>::Dot(const CompositeVector_<Scalar>& other, Scalar* result) const {
  /*
  Scalar tmp_result = 0.0;
  for (const auto& lcv=begin(); lcv!=end(); ++lcv) {
    if (other.HasComponent(*lcv)) {
      std::vector<Scalar> intermediate_result(ViewComponent(*lcv,false)->NumVectors(),0.0);
      int ierr = ViewComponent(*lcv, false)->Dot(*other.ViewComponent(*lcv,false),
                                                 &intermediate_result[0]);
      if (ierr) return ierr;
      
      for (int lcv_vector = 0; lcv_vector != NumVectors(*lcv); ++lcv_vector) {
        tmp_result += intermediate_result[lcv_vector];
      }
    }
  }
  *result = tmp_result;
  */
  mastervec_->Dot(*other.mastervec_, result); /// FIXME Why is this in both Block and Composite?  Should be only here, as this one knows that the names match?  No need for block to be more than a container, other than potentially PutScalar()
  return 0;
};


// -- this <- scalarA*A + scalarThis*this
template<typename Scalar>
CompositeVector_<Scalar>&
CompositeVector_<Scalar>::Update(Scalar scalarA, const CompositeVector_<Scalar>& A,
                                         Scalar scalarThis) {
  //  AMANZI_ASSERT(Map()->SubsetOf(*A.map_));
  for (const auto& name : *this) {
    if (A.HasComponent(name))
      GetComponent_(name, false).update(scalarA, A.GetComponent_(name,false), scalarThis);
  }
  return *this;
};


// -- this <- scalarA*A + scalarB*B + scalarThis*this
template<typename Scalar>
CompositeVector_<Scalar>&
CompositeVector_<Scalar>::Update(Scalar scalarA, const CompositeVector_<Scalar>& A,
                 Scalar scalarB, const CompositeVector_<Scalar>& B, Scalar scalarThis) {
  //  AMANZI_ASSERT(Map()->SubsetOf(*A.map_));
  //  AMANZI_ASSERT(Map()->SubsetOf(*B.map_));
  for (const auto& name : *this) {
    if (A.HasComponent(name) && B.HasComponent(name))
      GetComponent_(name, false).update(scalarA, A.GetComponent_(name,false),
                                         scalarB, B.GetComponent_(name,false), scalarThis);
  }
  return *this;
};


// -- this <- scalarAB * A@B + scalarThis*this  (@ is the elementwise product
template<typename Scalar>
int
CompositeVector_<Scalar>::Multiply(Scalar scalarAB, const CompositeVector_<Scalar>& A,
        const CompositeVector_<Scalar>& B, Scalar scalarThis) {
  //  AMANZI_ASSERT(Map()->SubsetOf(*A.map_));
  //  AMANZI_ASSERT(Map()->SubsetOf(*B.map_));
  int ierr = 0;
  for (const auto& name : *this) {
    if (A.HasComponent(name) && B.HasComponent(name)) {
      if (A.GetComponent_(name, false).getNumVectors() > 1) {
        Errors::Message message("Not implemented multiply: Tpetra does not provide elementwise multiply.");
        Exceptions::amanzi_throw(message);
      }
      GetComponent_(name, false).elementWiseMultiply(scalarAB, *A.GetComponent_(name,false).getVector(0),
                                          B.GetComponent_(name,false), scalarThis);
    }
  }
  return ierr;
};


// -- this <- scalarAB * B / A + scalarThis*this  (/ is the elementwise division
// template<typename Scalar>
// int
// CompositeVector_<Scalar>::ReciprocalMultiply(Scalar scalarAB, const CompositeVector_<Scalar>& A,
//                                         const CompositeVector_<Scalar>& B, Scalar scalarThis) {
//   AMANZI_ASSERT(Map()->SubsetOf(*A.map_));
//   AMANZI_ASSERT(Map()->SubsetOf(*B.map_));
//   int ierr = 0;
//   for (const auto& name : *this) {
//     if (A.HasComponent(name) && B.HasComponent(name))
//       ierr |= ViewComponent(name, false)->ReciprocalMultiply(scalarAB,
//               *A.ViewComponent(name,false), *B.ViewComponent(name,false), scalarThis);
//   }
//   return ierr;
// };


// // Mathematical operations
// // -- minimum value by component
// template<typename Scalar>
// void
// CompositeVector_<Scalar>::MinValue(std::map<std::string, Scalar>& value) const {
//   value.clear();

//   for (int n = 0; n != names_.size(); ++n) {
//     Scalar tmp(1e+50), value_loc[1];
//     const MultiVector_type& comp = *ViewComponent(names_[n]);

//     for (int i = 0; i != comp.NumVectors(); ++i) {
//       comp(i)->MinValue(value_loc);
//       tmp = std::min(tmp, value_loc[0]);
//     }
//     value[names_[n]] = tmp;
//   }
// };


// // -- maximum value by component
// template<typename Scalar>
// void
// CompositeVector_<Scalar>::MaxValue(std::map<std::string, Scalar>& value) const {
//   value.clear();

//   for (int n = 0; n != names_.size(); ++n) {
//     Scalar tmp(-1e+50), value_loc[1];
//     const MultiVector_type& comp = *ViewComponent(names_[n]);

//     for (int i = 0; i != comp.NumVectors(); ++i) {
//       comp(i)->MaxValue(value_loc);
//       tmp = std::max(tmp, value_loc[0]);
//     }
//     value[names_[n]] = tmp;
//   }
// };


// // -- mean value by component
// template<typename Scalar>
// void
// CompositeVector_<Scalar>::MeanValue(std::map<std::string, Scalar>& value) const {
//   value.clear();

//   for (int n = 0; n != names_.size(); ++n) {
//     int ni, nt(0);
//     Scalar tmp(0.0), value_loc[1];
//     const MultiVector_type& comp = *ViewComponent(names_[n]);

//     for (int i = 0; i != comp.NumVectors(); ++i) {
//       ni = comp(i)->GlobalLength(); 
//       comp(i)->MeanValue(value_loc);
//       tmp += value_loc[0] * ni;
//       nt += ni;
//     }
//     value[names_[n]] = tmp / nt;
//   }
// };


// -----------------------------------------------------------------------------
// Interpolate pressure ICs on cells to ICs for lambda (faces).
// -----------------------------------------------------------------------------
// void DeriveFaceValuesFromCellValues(CompositeVector& cv) {
//   if (cv.HasComponent("face")) {
//     cv.ScatterMasterToGhosted("cell");
//     const MultiVector_type& cv_c = *cv.ViewComponent("cell",true);
//     MultiVector_type& cv_f = *cv.ViewComponent("face",false);

//     int f_owned = cv_f.MyLength();
//     for (int f=0; f!=f_owned; ++f) {
//       AmanziMesh::Entity_ID_List cells;
//       cv.Mesh()->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
//       int ncells = cells.size();

//       Scalar face_value = 0.0;
//       for (int n=0; n!=ncells; ++n) {
//         face_value += cv_c[0][cells[n]];
//       }
//       cv_f[0][f] = face_value / ncells;
//     }
//   }
//   else if (cv.HasComponent("boundary_face")) {
//     const MultiVector_type& cv_c = *cv.ViewComponent("cell",true);
//     MultiVector_type& cv_f = *cv.ViewComponent("boundary_face",false);

//     const Map_type& fb_map = cv.Mesh()->exterior_face_map(false);
//     const Map_type& f_map = cv.Mesh()->face_map(false);

//     int fb_owned = cv_f.MyLength();
//     for (int fb=0; fb!=fb_owned; ++fb) {
//       AmanziMesh::Entity_ID_List cells;

//       int f_gid = fb_map.GID(fb);
//       int f_lid = f_map.LID(f_gid);
      
//       cv.Mesh()->face_get_cells(f_lid, AmanziMesh::Parallel_type::ALL, &cells);
//       int ncells = cells.size();

//       AMANZI_ASSERT((ncells==1));

//       Scalar face_value = cv_c[0][cells[0]];
//       cv_f[0][fb] = face_value;
//     }
//   }
// }


// // -----------------------------------------------------------------------------
// // Non-member function: extension of a vector via data copy.
// // -----------------------------------------------------------------------------
// void AddComponent(Teuchos::RCP<CompositeVector> cv,
//                   const std::string& name, AmanziMesh::Entity_kind kind, int dim) {
//   // copy construct the CVS making it not owned and add the new component
//   CompositeMap new_space(cv->Map());
//   new_space.SetOwned(false);
//   new_space.AddComponent(name, kind, dim);

//   // create the new vector and copy data
//   Teuchos::RCP<CompositeVector> new_cv = Teuchos::rcp(new CompositeVector(new_space));
//   bool ghost = new_space.Ghosted();
//   std::vector<std::string>::const_iterator it;
//   for (it = cv->Map().begin(); it != cv->Map().end(); ++it) {
//     auto data1 = cv->ViewComponent(*it, ghost);
//     auto data2 = new_cv->ViewComponent(*it, ghost);
//     *data2 = *data1;
//   }

//   // replace the vector
//   cv = new_cv;
// }

} // namespace

#endif
