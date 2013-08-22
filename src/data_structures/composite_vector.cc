/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Ethan Coon

   Interface for CompositeVector, an implementation of a slightly improved
   Epetra_MultiVector which spans multiple simplices and knows how to
   communicate itself.

   NOTE: All CompositeVector data is NOT initialized to zero!
   ------------------------------------------------------------------------- */

#include "dbc.hh"
#include "errors.hh"
#include "composite_vector.hh"

namespace Amanzi {

// Constructor
CompositeVector::CompositeVector(
        const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
        std::vector<std::string> names,
        std::vector<AmanziMesh::Entity_kind> locations,
        std::vector<int> num_dofs,
        bool ghosted) :
    mesh_(mesh),
    names_(names),
    locations_(locations),
    ghosted_(ghosted),
    created_(false) {

  num_components_ = locations_.size();

  // generate the reverse map of name->index
  for (unsigned int i=0; i!=num_components_; ++i) {
    indexmap_[names_[i]] = i;
  }

  // generate the master's maps
  std::vector<Teuchos::RCP<const Epetra_Map> > mastermaps(num_components_);
  for (unsigned int i=0; i!=num_components_; ++i) {
    if (locations_[i] == AmanziMesh::CELL) {
      mastermaps[i] = Teuchos::rcpFromRef(mesh_->cell_map(false));
    } else if (locations_[i] == AmanziMesh::FACE) {
      mastermaps[i] = Teuchos::rcpFromRef(mesh_->face_map(false));
    } else if (locations_[i] == AmanziMesh::NODE) {
      mastermaps[i] = Teuchos::rcpFromRef(mesh_->node_map(false));
    }
  }

  // create the master BlockVector
  mastervec_ = Teuchos::rcp(new BlockVector(comm(), names_, mastermaps, num_dofs));

  // do the same for the ghosted Vector, if necessary
  if (ghosted) {
    // generate the ghost's maps
    std::vector<Teuchos::RCP<const Epetra_Map> > ghostmaps(num_components_);
    for (unsigned int i=0; i!=num_components_; ++i) {
      if (locations_[i] == AmanziMesh::CELL) {
        ghostmaps[i] = Teuchos::rcpFromRef(mesh_->cell_map(true));
      } else if (locations_[i] == AmanziMesh::FACE) {
        ghostmaps[i] = Teuchos::rcpFromRef(mesh_->face_map(true));
      } else if (locations_[i] == AmanziMesh::NODE) {
        ghostmaps[i] = Teuchos::rcpFromRef(mesh_->node_map(true));
      }
    }

    // create the ghost BlockVector
    ghostvec_ =
      Teuchos::rcp(new BlockVector(comm(), names_, ghostmaps, num_dofs));
  } else {
    ghostvec_ = mastervec_;
  }
};


// copy constructor
CompositeVector::CompositeVector(const CompositeVector& other,
        ConstructMode mode) :
    mesh_(other.mesh_),
    names_(other.names_),
    locations_(other.locations_),
    importers_(other.importers_),
    ghosted_(other.ghosted_),
    num_components_(other.num_components_),
    indexmap_(other.indexmap_) {

  if (mode == CONSTRUCT_WITH_OLD_DATA) {
    // construct with the same old data, just make a new copy of pointers
    mastervec_ = Teuchos::rcp(new BlockVector(*other.get_master_vector(),
            CONSTRUCT_WITH_OLD_DATA));
    if (ghosted_) {
      ghostvec_ = Teuchos::rcp(new BlockVector(*other.get_ghost_vector(),
            CONSTRUCT_WITH_OLD_DATA));
    } else {
      ghostvec_ = mastervec_;
    }
    created_ = true;

  } else {
    // copy construct the master block vector, without data
    mastervec_ = Teuchos::rcp(new BlockVector(*other.get_master_vector(),
            CONSTRUCT_WITHOUT_DATA));
    // same for the ghost vector, if needed
    if (ghosted_) {
      ghostvec_ = Teuchos::rcp(new BlockVector(*other.get_ghost_vector(),
              CONSTRUCT_WITHOUT_DATA));
    } else {
      ghostvec_ = mastervec_;
    }

    // if new data is requested, create the vector
    if (mode == CONSTRUCT_WITH_NEW_DATA) CreateData();
  }
};


// Sets sizes of vectors, instantiates Epetra_Vectors, and preps for lazy
// creation of everything else.
void CompositeVector::CreateData() {
  if (importers_.size() == 0) {
    importers_.resize(num_components_, Teuchos::null);
  }

  // Create the ghost vector.  Note this is also the master vector if not ghosted.
  ghostvec_->CreateData();

  // If the vector is ghosted, create the master from views of the ghost.
  if (ghosted_) {
    for (name_iterator name=begin(); name!=end(); ++name) {
      // get the ghost component's data
      Teuchos::RCP<Epetra_MultiVector> g_comp = ghostvec_->ViewComponent(*name);
      double** data;
      g_comp->ExtractView(&data);

      // create the master component
      Teuchos::RCP<Epetra_MultiVector> m_comp =
        Teuchos::rcp(new Epetra_MultiVector(View, *mastervec_->map(*name), data,
                mastervec_->num_dofs(*name)));

      // push it back into the master vec
      mastervec_->SetComponent(*name, m_comp);
    }
  }
  created_ = true;
};


CompositeVector& CompositeVector::operator=(const CompositeVector& other) {
  if (this != &other) {

#ifdef ENABLE_DBC
    // check the meta-data is compatible
    if ((num_components_ > other.num_components_) || (mesh_ != other.mesh_)) {
      Errors::Message message("Attempted assignment of non-compatible CompositeVectors.");
      Exceptions::amanzi_throw(message);
    }
    for (name_iterator name=begin(); name!=end(); ++name) {
      // (Implicitly) check that each of this's names is in other's names, and
      // check that the num_dofs and locations of each name matches.
      ASSERT(num_dofs(*name) == other.num_dofs(*name));
      ASSERT(location(*name) == other.location(*name));
    }
#endif

    if (ghosted() && other.ghosted()) {
      // If both are ghosted, copy the ghosted vector.
      for (name_iterator name=begin(); name!=end(); ++name) {
        Teuchos::RCP<Epetra_MultiVector> comp = ViewComponent(*name, true);
        Teuchos::RCP<const Epetra_MultiVector> othercomp = other.ViewComponent(*name, true);
        *comp = *othercomp;
      }

    } else {
      // Copy the non-ghosted data.  NOTE: any ghosted data is undefined!
      for (name_iterator name=begin(); name!=end(); ++name) {
        Teuchos::RCP<Epetra_MultiVector> comp = ViewComponent(*name, false);
        Teuchos::RCP<const Epetra_MultiVector> othercomp = other.ViewComponent(*name, false);
        *comp = *othercomp;
      }
    }
  }
  return *this;
};


void CompositeVector::AssertCreatedOrDie_() const {
#ifdef ENABLE_DBC
  if (!created_) {
    Errors::Message message("CreateData has not been called on this CompositeVector.");
    Exceptions::amanzi_throw(message);
  }
#endif
};


// view data
// -- Access a view of a single component's data.
// Ghosted views are simply the vector itself, while non-ghosted views are
// lazily generated.
Teuchos::RCP<const Epetra_MultiVector>
CompositeVector::ViewComponent(std::string name, bool ghosted) const {
  AssertCreatedOrDie_();
  if (ghosted) {
    return ghostvec_->ViewComponent(name);
  } else {
    return mastervec_->ViewComponent(name);
  }
};


Teuchos::RCP<Epetra_MultiVector>
CompositeVector::ViewComponent(std::string name, bool ghosted) {
  AssertCreatedOrDie_();
  if (ghosted) {
    return ghostvec_->ViewComponent(name);
  } else {
    return mastervec_->ViewComponent(name);
  }
};


// Set data by pointer if possible, otherwise by copy.
void CompositeVector::SetComponent(std::string name,
        const Teuchos::RCP<Epetra_MultiVector>& data) {
  if (ghostvec_->map(name)->SameAs(data->Map())) {
    // setting the ghost vec -- drop in the data in the ghost
    ghostvec_->SetComponent(name, data);

    // and create a new view for the master
    double** vals;
    data->ExtractView(&vals);
    Teuchos::RCP<Epetra_MultiVector> m_comp =
      Teuchos::rcp(new Epetra_MultiVector(View, *mastervec_->map(name), vals,
              mastervec_->num_dofs(name)));
    mastervec_->SetComponent(name, m_comp);

  } else if (mastervec_->map(name)->SameAs(data->Map())) {
    *mastervec_->ViewComponent(name) = *data;
  } else {
    Errors::Message message("Attempted set of non-compatible Component.");
    Exceptions::amanzi_throw(message);
  }
};


// communicate
// -- Scatter master values to ghosted values.
// Modes shown in Epetra_CombineMode.h, but the default is Insert, which
// overwrites the current ghost value with the (unique) new master value.
void CompositeVector::ScatterMasterToGhosted(Epetra_CombineMode mode) const {
  // NOTE: allowing const is a hack to allow non-owning PKs to nonetheless
  // update ghost cells, which may be necessary for their discretization

  for (name_iterator name=begin(); name!=end(); ++name) {
    ScatterMasterToGhosted(*name, mode);
  }
};


void CompositeVector::ScatterMasterToGhosted(std::string name,
        Epetra_CombineMode mode) const {
  // NOTE: allowing const is a hack to allow non-owning PKs to nonetheless
  // update ghost cells, which may be necessary for their discretization

#ifdef HAVE_MPI
  if (ghosted_) {
    // check for and create the importer if needed
    if (importers_[index_(name)] == Teuchos::null) {
      Teuchos::RCP<const Epetra_Map> target_map = ghostvec_->map(name);
      Teuchos::RCP<const Epetra_Map> source_map = mastervec_->map(name);
      importers_[index_(name)] =
        Teuchos::rcp(new Epetra_Import(*target_map, *source_map));
    }

    // communicate
    Teuchos::RCP<Epetra_MultiVector> g_comp =
      ghostvec_->ViewComponent(name);
    Teuchos::RCP<const Epetra_MultiVector> m_comp =
      mastervec_->ViewComponent(name);
    g_comp->Import(*m_comp, *importers_[index_(name)], mode);
  }
#endif
};


// -- Combine ghosted values back to master values.
// Modes shown in Epetra_CombineMode.h, but the default is Add,
// where off-process values are first summed into the on-process value.
void CompositeVector::GatherGhostedToMaster(Epetra_CombineMode mode) {
  for (name_iterator name=begin(); name!=end(); ++name) {
    GatherGhostedToMaster(*name, mode);
  }
};


void CompositeVector::GatherGhostedToMaster(std::string name,
        Epetra_CombineMode mode) {
#ifdef HAVE_MPI
  if (ghosted_) {
    // check for and create the importer if needed
    if (importers_[index_(name)] == Teuchos::null) {
      Teuchos::RCP<const Epetra_Map> target_map = ghostvec_->map(name);
      Teuchos::RCP<const Epetra_Map> source_map = mastervec_->map(name);
      importers_[index_(name)] =
        Teuchos::rcp(new Epetra_Import(*target_map, *source_map));
    }

    // communicate
    Teuchos::RCP<const Epetra_MultiVector> g_comp =
      ghostvec_->ViewComponent(name);
    Teuchos::RCP<Epetra_MultiVector> m_comp =
      mastervec_->ViewComponent(name);
    m_comp->Export(*g_comp, *importers_[index_(name)], mode);
  }
#endif
};

} // namespace

