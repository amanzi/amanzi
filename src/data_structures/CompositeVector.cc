/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

ATS and Amanzi

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for CompositeVector, an implementation of a block Tpetra
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
------------------------------------------------------------------------- */

#include "Tpetra_Vector.hpp"

#include "dbc.hh"
#include "errors.hh"
#include "CompositeVector.hh"

namespace Amanzi {

// Constructor
CompositeVector::CompositeVector(const CompositeVectorSpace& space) :
    map_(Teuchos::rcp(new CompositeVectorSpace(space))),
    names_(space.names_),
    indexmap_(space.indexmap_),
    ghosted_(space.ghosted_)
{
  InitMap_(*map_);
  CreateData_();
}


CompositeVector::CompositeVector(const CompositeVectorSpace& space, bool ghosted) :
    map_(Teuchos::rcp(new CompositeVectorSpace(space,ghosted))),
    names_(space.names_),
    indexmap_(space.indexmap_),
    ghosted_(ghosted)
{
  InitMap_(*map_);
  CreateData_();
}


CompositeVector::CompositeVector(const CompositeVector& other,
                                 InitMode mode) :
    map_(Teuchos::rcp(new CompositeVectorSpace(*other.map_))),
    names_(other.names_),
    indexmap_(other.indexmap_),
    ghosted_(other.ghosted_)
{
  InitMap_(*map_);
  CreateData_();
  InitData_(other, mode);
}

CompositeVector::CompositeVector(const CompositeVector& other, bool ghosted,
                                 InitMode mode) :
    map_(Teuchos::rcp(new CompositeVectorSpace(*other.map_,ghosted))),
    names_(other.names_),
    indexmap_(other.indexmap_),
    ghosted_(ghosted)
{
  InitMap_(*map_);
  CreateData_();
  InitData_(other, mode);
}



void CompositeVector::InitMap_(const CompositeVectorSpace& space) {
  // generate the master's maps
  std::vector<Teuchos::RCP<const Map_type> > mastermaps;
  for (CompositeVectorSpace::name_iterator name=space.begin(); name!=space.end(); ++name) {
    if (space.Location(*name) == AmanziMesh::CELL) {
      mastermaps.push_back((Mesh()->cell_map(false)));
    } else if (space.Location(*name) == AmanziMesh::FACE) {
      mastermaps.push_back((Mesh()->face_map(false)));
    } else if (space.Location(*name) == AmanziMesh::NODE) {
      mastermaps.push_back((Mesh()->node_map(false)));
    } else if (space.Location(*name) == AmanziMesh::EDGE) {
      mastermaps.push_back((Mesh()->edge_map(false)));
    } else if (space.Location(*name) == AmanziMesh::BOUNDARY_FACE) {
      mastermaps.push_back((Mesh()->exterior_face_map(false)));
    }
  }

  // create the master BlockVector
  mastervec_ = Teuchos::rcp(new BlockVector(Comm(), names_,
          mastermaps, space.num_dofs_));

  // do the same for the ghosted Vector, if necessary
  if (space.Ghosted()) {
    // generate the ghost's maps
    std::vector<Teuchos::RCP<const Map_type> > ghostmaps;
    for (CompositeVectorSpace::name_iterator name=space.begin(); name!=space.end(); ++name) {
      if (space.Location(*name) == AmanziMesh::CELL) {
        ghostmaps.push_back((Mesh()->cell_map(true)));
      } else if (space.Location(*name) == AmanziMesh::FACE) {
        ghostmaps.push_back((Mesh()->face_map(true)));
      } else if (space.Location(*name) == AmanziMesh::NODE) {
        ghostmaps.push_back((Mesh()->node_map(true)));
      } else if (space.Location(*name) == AmanziMesh::EDGE) {
        ghostmaps.push_back((Mesh()->edge_map(true)));
      } else if (space.Location(*name) == AmanziMesh::BOUNDARY_FACE) {
        ghostmaps.push_back((Mesh()->exterior_face_map(true)));
      }
    }

    // create the ghost BlockVector
    ghostvec_ = Teuchos::rcp(new BlockVector(Comm(), names_,
            ghostmaps, space.num_dofs_));
  } else {
    ghostvec_ = mastervec_;
  }
};


// Initialize data
void CompositeVector::InitData_(const CompositeVector& other, InitMode mode) {
  // Trilinos inits to 0
  //  if (mode == INIT_MODE_ZERO) {
  //    PutScalar(0.);
  if (mode == INIT_MODE_COPY) {
    *this = other;
  }
}

// Sets sizes of vectors, instantiates Vectors, and preps for lazy
// creation of everything else.
void CompositeVector::CreateData_() {
  if (!Mesh().get()) {
    Errors::Message message("CompositeVector: construction called with no mesh.");
    Exceptions::amanzi_throw(message);
  }
  if (NumComponents() == 0) {
    Errors::Message message("CompositeVector: construction called with no components.");
    Exceptions::amanzi_throw(message);
  }

  if (importers_.size() == 0) {
    importers_.resize(NumComponents(), Teuchos::null);
  }

  // Create the ghost vector.  Note this is also the master vector if not ghosted.
  ghostvec_->CreateData();

  // If the vector is ghosted, create the master from views of the ghost.
  if (ghosted_) {
    for (name_iterator name=begin(); name!=end(); ++name) {
      // get the ghost component's data
      auto g_comp = ghostvec_->GetComponent(*name);
      auto m_comp = g_comp->offsetViewNonConst(mastervec_->ComponentMap(*name),0);
      mastervec_->SetComponent(*name, m_comp);
    }
  }
  PutScalarMasterAndGhosted(0.);
};


CompositeVector& CompositeVector::operator=(const CompositeVector& other) {
  if (this != &other) {
    AMANZI_ASSERT(Map().SubsetOf(other.Map()));

    if (Ghosted() && other.Ghosted()) {
      // If both are ghosted, copy the ghosted vector.
      for (name_iterator name=begin(); name!=end(); ++name) {
        auto& comp = GetComponent_(*name, true);
        const auto& othercomp = other.GetComponent_(*name, true);
        comp = othercomp;
      }

    } else {
      // Copy the non-ghosted data.  NOTE: any ghosted data is undefined!
      for (name_iterator name=begin(); name!=end(); ++name) {
        auto& comp = GetComponent_(*name, false);
        const auto& othercomp = other.GetComponent_(*name, false);
        comp = othercomp;
      }
    }
  }
  return *this;
};


// view data
// -- Access a view of a single component's data.
// Ghosted views are simply the vector itself, while non-ghosted views are
// lazily generated.
const MultiVector_type&
CompositeVector::GetComponent_(const std::string& name, bool ghosted) const {

/* -- FIXME EPETRA --> TPETRA
  if (name == std::string("boundary_face")) {
    if (!mastervec_->HasComponent("boundary_face") &&
        mastervec_->HasComponent("face")) {
      ApplyVandelay_();
      return vandelay_vector_;
    }
  }
  -- END FIXME EPETRA --> TPETRA */
  if (ghosted) {
    return *ghostvec_->GetComponent(name);
  } else {
    return *mastervec_->GetComponent(name);
  }
};


MultiVector_type&
CompositeVector::GetComponent_(const std::string& name, bool ghosted) {
/* -- FIXME EPETRA --> TPETRA
  if (name == std::string("boundary_face")) {
    if (!mastervec_->HasComponent("boundary_face") &&
        mastervec_->HasComponent("face")) {
      ApplyVandelay_();
      return vandelay_vector_;
    }
  }
  -- END FIXME EPETRA --> TPETRA */
  if (ghosted) {
    return *ghostvec_->GetComponent(name);
  } else {
    return *mastervec_->GetComponent(name);
  }
};


// Set data by pointer if possible, otherwise by copy.
void CompositeVector::SetComponent(const std::string& name,
        MultiVector_ptr_type data) {
  if (ghostvec_->ComponentMap(name)->isSameAs(*data->getMap())) {
    // setting the ghost vec -- drop in the data in the ghost
    ghostvec_->SetComponent(name, data);

    // and create a new view for the master
    auto m_comp = data->offsetViewNonConst(mastervec_->ComponentMap(name),0);
    mastervec_->SetComponent(name, m_comp);

  } else if (mastervec_->ComponentMap(name)->isSameAs(*data->getMap())) {
    *mastervec_->GetComponent(name) = *data;
  } else {
    Errors::Message message("Attempted set of non-compatible Component.");
    Exceptions::amanzi_throw(message);
  }
};


// -- Scatter master values to ghosted values.
// Modes shown in Epetra_CombineMode.h, but the default is Insert, which
// overwrites the current ghost value with the (unique) new master value.
void CompositeVector::ScatterMasterToGhosted(bool force) const {
  for (name_iterator name=begin(); name!=end(); ++name) {
    ScatterMasterToGhosted(*name, force);
  }
};


void
CompositeVector::ScatterMasterToGhosted(const std::string& name, bool force) const {
  // NOTE: allowing const is a hack to allow non-owning PKs to nonetheless
  // update ghost cells, which may be necessary for their discretization
#ifdef HAVE_MPI
  if (ghosted_) {
    // check for and create the importer if needed
    if (importers_[Index_(name)] == Teuchos::null) {
      Teuchos::RCP<const Map_type> target_map = ghostvec_->ComponentMap(name);
      Teuchos::RCP<const Map_type> source_map = mastervec_->ComponentMap(name);
      importers_[Index_(name)] =
        Teuchos::rcp(new Tpetra::Import<>(source_map, target_map));
    }

    // communicate
    auto g_comp = ghostvec_->GetComponent(name);
    auto m_comp = mastervec_->GetComponent(name);
    g_comp->doImport(*m_comp, *importers_[Index_(name)], Tpetra::INSERT);
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
CompositeVector::ScatterMasterToGhosted(Tpetra::CombineMode mode) const {
  for (name_iterator name=begin(); name!=end(); ++name) {
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
CompositeVector::ScatterMasterToGhosted(const std::string& name,
        Tpetra::CombineMode mode) const {

#ifdef HAVE_MPI
  if (ghosted_) {
    // check for and create the importer if needed
    if (importers_[Index_(name)] == Teuchos::null) {
      Teuchos::RCP<const Map_type> target_map = ghostvec_->ComponentMap(name);
      Teuchos::RCP<const Map_type> source_map = mastervec_->ComponentMap(name);
      importers_[Index_(name)] =
        Teuchos::rcp(new Tpetra::Import<>(source_map, target_map));
    }

    // communicate
    auto g_comp = ghostvec_->GetComponent(name);
    auto m_comp = mastervec_->GetComponent(name);
    g_comp->doImport(*m_comp, *importers_[Index_(name)], Tpetra::INSERT);
  }
#endif
};




// -- Combine ghosted values back to master values.
// Modes shown in Epetra_CombineMode.h, but the default is Add,
// where off-process values are first summed into the on-process value.
void CompositeVector::GatherGhostedToMaster(Tpetra::CombineMode mode) {
  for (name_iterator name=begin(); name!=end(); ++name) {
    GatherGhostedToMaster(*name, mode);
  }
};


void CompositeVector::GatherGhostedToMaster(const std::string& name,
        Tpetra::CombineMode mode) {
#ifdef HAVE_MPI
  if (ghosted_) {
    // check for and create the importer if needed
    if (importers_[Index_(name)] == Teuchos::null) {
      auto target_map = ghostvec_->ComponentMap(name);
      auto source_map = mastervec_->ComponentMap(name);
      importers_[Index_(name)] =
        Teuchos::rcp(new Tpetra::Import<>(source_map, target_map));
    }

    // communicate
    auto g_comp = ghostvec_->GetComponent(name);
    auto m_comp = mastervec_->GetComponent(name);
    m_comp->doExport(*g_comp, *importers_[Index_(name)], mode);
  }
#endif
};


// Vandelay operations
void CompositeVector::CreateVandelay_() const {
  vandelay_importer_ = Mesh()->exterior_face_importer();
  vandelay_vector_ = Teuchos::rcp(new MultiVector_type(Mesh()->exterior_face_map(false),
          mastervec_->NumVectors("face"), false));
}

void CompositeVector::ApplyVandelay_() const {
  if (vandelay_importer_ == Teuchos::null) {
    CreateVandelay_();
  }
  vandelay_vector_->doImport(GetComponent_("face",false), *vandelay_importer_, Tpetra::INSERT);
}


// return non-empty importer
Import_ptr_type
CompositeVector::importer(const std::string& name) {
  if (importers_[Index_(name)] == Teuchos::null) {
    Teuchos::RCP<const Map_type> target_map = ghostvec_->ComponentMap(name);
    Teuchos::RCP<const Map_type> source_map = mastervec_->ComponentMap(name);
    importers_[Index_(name)] =
      Teuchos::rcp(new Tpetra::Import<>(source_map, target_map));
  }
  return importers_[Index_(name)];
}


// Mathematical operations
// -- result <- other \dot this
int CompositeVector::Dot(const CompositeVector& other, double* result) const {
  /*
  double tmp_result = 0.0;
  for (name_iterator lcv=begin(); lcv!=end(); ++lcv) {
    if (other.HasComponent(*lcv)) {
      std::vector<double> intermediate_result(ViewComponent(*lcv,false)->NumVectors(),0.0);
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
CompositeVector& CompositeVector::Update(double scalarA, const CompositeVector& A,
                                         double scalarThis) {
  //  AMANZI_ASSERT(map_->SubsetOf(*A.map_));
  for (name_iterator lcv=begin(); lcv!=end(); ++lcv) {
    if (A.HasComponent(*lcv))
      GetComponent_(*lcv, false).update(scalarA, A.GetComponent_(*lcv,false), scalarThis);
  }
  return *this;
};


// -- this <- scalarA*A + scalarB*B + scalarThis*this
CompositeVector& CompositeVector::Update(double scalarA, const CompositeVector& A,
                 double scalarB, const CompositeVector& B, double scalarThis) {
  //  AMANZI_ASSERT(map_->SubsetOf(*A.map_));
  //  AMANZI_ASSERT(map_->SubsetOf(*B.map_));
  for (name_iterator lcv=begin(); lcv!=end(); ++lcv) {
    if (A.HasComponent(*lcv) && B.HasComponent(*lcv))
      GetComponent_(*lcv, false).update(scalarA, A.GetComponent_(*lcv,false),
                                         scalarB, B.GetComponent_(*lcv,false), scalarThis);
  }
  return *this;
};


// -- this <- scalarAB * A@B + scalarThis*this  (@ is the elementwise product
int CompositeVector::Multiply(double scalarAB, const CompositeVector& A,
        const CompositeVector& B, double scalarThis) {
  //  AMANZI_ASSERT(map_->SubsetOf(*A.map_));
  //  AMANZI_ASSERT(map_->SubsetOf(*B.map_));
  int ierr = 0;
  for (name_iterator lcv=begin(); lcv!=end(); ++lcv) {
    if (A.HasComponent(*lcv) && B.HasComponent(*lcv)) {
      if (A.GetComponent_(*lcv, false).getNumVectors() > 1) {
        Errors::Message message("Not implemented multiply: Tpetra does not provide elementwise multiply.");
        Exceptions::amanzi_throw(message);
      }
      GetComponent_(*lcv, false).elementWiseMultiply(scalarAB, *A.GetComponent_(*lcv,false).getVector(0),
                                          B.GetComponent_(*lcv,false), scalarThis);
    }
  }
  return ierr;
};


// -- this <- scalarAB * B / A + scalarThis*this  (/ is the elementwise division
// int CompositeVector::ReciprocalMultiply(double scalarAB, const CompositeVector& A,
//                                         const CompositeVector& B, double scalarThis) {
//   AMANZI_ASSERT(map_->SubsetOf(*A.map_));
//   AMANZI_ASSERT(map_->SubsetOf(*B.map_));
//   int ierr = 0;
//   for (name_iterator lcv=begin(); lcv!=end(); ++lcv) {
//     if (A.HasComponent(*lcv) && B.HasComponent(*lcv))
//       ierr |= ViewComponent(*lcv, false)->ReciprocalMultiply(scalarAB,
//               *A.ViewComponent(*lcv,false), *B.ViewComponent(*lcv,false), scalarThis);
//   }
//   return ierr;
// };


// // Mathematical operations
// // -- minimum value by component
// void CompositeVector::MinValue(std::map<std::string, double>& value) const {
//   value.clear();

//   for (int n = 0; n != names_.size(); ++n) {
//     double tmp(1e+50), value_loc[1];
//     const MultiVector_type& comp = *ViewComponent(names_[n]);

//     for (int i = 0; i != comp.NumVectors(); ++i) {
//       comp(i)->MinValue(value_loc);
//       tmp = std::min(tmp, value_loc[0]);
//     }
//     value[names_[n]] = tmp;
//   }
// };


// // -- maximum value by component
// void CompositeVector::MaxValue(std::map<std::string, double>& value) const {
//   value.clear();

//   for (int n = 0; n != names_.size(); ++n) {
//     double tmp(-1e+50), value_loc[1];
//     const MultiVector_type& comp = *ViewComponent(names_[n]);

//     for (int i = 0; i != comp.NumVectors(); ++i) {
//       comp(i)->MaxValue(value_loc);
//       tmp = std::max(tmp, value_loc[0]);
//     }
//     value[names_[n]] = tmp;
//   }
// };


// // -- mean value by component
// void CompositeVector::MeanValue(std::map<std::string, double>& value) const {
//   value.clear();

//   for (int n = 0; n != names_.size(); ++n) {
//     int ni, nt(0);
//     double tmp(0.0), value_loc[1];
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

//       double face_value = 0.0;
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

//       double face_value = cv_c[0][cells[0]];
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
//   CompositeVectorSpace new_space(cv->Map());
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

