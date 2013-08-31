/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Implementation for a Field.

Field also stores some basic metadata for Vis, checkpointing, etc.
------------------------------------------------------------------------- */

#include <string>

#include "errors.hh"
#include "composite_vector.hh"
#include "composite_vector_factory.hh"
#include "composite_vector_function.hh"
#include "composite_vector_function_factory.hh"

#include "Field.hh"
#include "Field_CompositeVector.hh"

namespace Amanzi {

Field_CompositeVector::Field_CompositeVector(std::string fieldname, std::string owner) :
    Field::Field(fieldname, owner), data_() {
  type_ = COMPOSITE_VECTOR_FIELD;
};

Field_CompositeVector::Field_CompositeVector(std::string fieldname, std::string owner,
        const std::vector<std::vector<std::string> >& subfield_names) :
    Field::Field(fieldname, owner),
    data_(),
    subfield_names_(subfield_names) {
  type_ = COMPOSITE_VECTOR_FIELD;
};

Field_CompositeVector::Field_CompositeVector(std::string fieldname, std::string owner,
                                           Teuchos::RCP<CompositeVector>& data) :
    Field::Field(fieldname, owner), data_(data) {
  type_ = COMPOSITE_VECTOR_FIELD;
};

// copy constructor:
Field_CompositeVector::Field_CompositeVector(const Field_CompositeVector& other) :
    Field::Field(other),
    subfield_names_(other.subfield_names_) {
  data_ = Teuchos::rcp(new CompositeVector(*other.data_));
};

// Virtual copy constructor
Teuchos::RCP<Field> Field_CompositeVector::Clone() const {
  return Teuchos::rcp(new Field_CompositeVector(*this));
}

// Virtual copy constructor with non-empty name
Teuchos::RCP<Field> Field_CompositeVector::Clone(std::string fieldname) const {
  Teuchos::RCP<Field_CompositeVector> other = Teuchos::rcp(new Field_CompositeVector(*this));
  other->fieldname_ = fieldname;
  return other;
};

// Virtual copy constructor with non-empty name
Teuchos::RCP<Field> Field_CompositeVector::Clone(std::string fieldname, std::string owner) const {
  Teuchos::RCP<Field_CompositeVector> other = Teuchos::rcp(new Field_CompositeVector(*this));
  other->fieldname_ = fieldname;
  other->owner_ = owner;
  return other;
};

// Create the data
void Field_CompositeVector::CreateData() {
  data_->CreateData();
}

// write-access to the data
Teuchos::RCP<CompositeVector> Field_CompositeVector::GetFieldData() {
  return data_;
};

// Overwrite data by pointer, not copy
void Field_CompositeVector::SetData(const Teuchos::RCP<CompositeVector>& data) {
  data_ = data;
};

// Set data by copy.
void Field_CompositeVector::SetData(const CompositeVector& data) {
  *data_ = data;
};

// Initialize the data using a series of possible options from the input file.
void Field_CompositeVector::Initialize(Teuchos::ParameterList& plist) {
  // ------ protect against unset names -----
  EnsureSubfieldNames_();

  // ------ Try to set values from a restart file -----
  if (plist.isParameter("restart file")) {
    std::string filename = plist.get<string>("restart file");
    ReadCheckpoint_(filename);
    set_initialized();
    return;
  }

  // ------ Try to set cell values from a restart file -----
  if (plist.isParameter("cells from file")) {
    std::string filename = plist.get<string>("cells from file");
    ReadCellsFromCheckpoint_(filename);
    return;
  }

  // ------ Set values using a constant -----
  if (plist.isParameter("constant")) {
    double value = plist.get<double>("constant");
    data_->PutScalar(value);
    set_initialized();
    return;
  }

  // ------ Set values using a function -----
  if (plist.isSublist("function")) {
    Teuchos::ParameterList func_plist = plist.sublist("function");

    // -- potential use of a mapping operator first --
    bool map_normal = plist.get<bool>("dot with normal", false);
    if (map_normal) {
      // map_normal take a vector and dots it with face normals
      ASSERT(data_->num_components() == 1); // one comp
      ASSERT(data_->has_component("face")); // is named face
      ASSERT(data_->location("face") == AmanziMesh::FACE); // is on face
      ASSERT(data_->num_dofs("face") == 1);  // and is scalar

      // create a vector on faces of the appropriate dimension
      CompositeVectorFactory fac;
      fac.SetMesh(data_->mesh())
          ->SetComponent("face",AmanziMesh::FACE,data_->mesh()->space_dimension());
      Teuchos::RCP<CompositeVector> vel_vec = fac.CreateVector(false);
      vel_vec->CreateData();

      // Evaluate the velocity function
      Teuchos::RCP<Functions::CompositeVectorFunction> func =
          Functions::CreateCompositeVectorFunction(func_plist, *vel_vec);
      func->Compute(0., vel_vec.ptr());

      // Dot the velocity with the normal
      unsigned int nfaces_owned = data_->mesh()
          ->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);

      Epetra_MultiVector& dat_f = *data_->ViewComponent("face",false);
      const Epetra_MultiVector& vel_f = *vel_vec->ViewComponent("face",false);
      int dim = data_->mesh()->space_dimension();


      AmanziGeometry::Point vel(dim);
      for (unsigned int f=0; f!=nfaces_owned; ++f) {
        AmanziGeometry::Point normal = data_->mesh()->face_normal(f);
        if (dim == 2) {
          vel.set(vel_f[0][f], vel_f[1][f]);
        } else if (dim == 3) {
          vel.set(vel_f[0][f], vel_f[1][f], vel_f[2][f]);
        } else {
          ASSERT(0);
        }
        dat_f[0][f] = vel * normal;
      }
      set_initialized();
      return;

    } else {
      // no map, just evaluate the function
      Teuchos::RCP<Functions::CompositeVectorFunction> func =
          Functions::CreateCompositeVectorFunction(func_plist, *data_);
      func->Compute(0.0, data_.ptr());
      set_initialized();
      return;
    }
  }
};


void Field_CompositeVector::WriteVis(const Teuchos::Ptr<Visualization>& vis) {
  if (io_vis_ && (vis->mesh() == data_->mesh())) {
    EnsureSubfieldNames_();

    // loop over the components and dump them to the vis file if possible
    int i = 0;
    for (CompositeVector::name_iterator compname=data_->begin();
         compname!=data_->end(); ++compname) {
      // check that this vector is a cell vector (currently this is the only
      // type of vector we can visualize
      if (data_->location(*compname) == AmanziMesh::CELL) {
        // get the MultiVector that should be dumped
        Teuchos::RCP<Epetra_MultiVector> v = data_->ViewComponent(*compname, false);

        // construct the name for vis
        std::vector< std::string > vis_names(subfield_names_[i]);
        for (unsigned int j = 0; j!=subfield_names_[i].size(); ++j) {
          vis_names[j] = fieldname_ + std::string(".") + *compname
            + std::string(".") + subfield_names_[i][j];
        }
        vis->WriteVector(*v, vis_names);
      }
      i++;
    }
  }
};


void Field_CompositeVector::WriteCheckpoint(const Teuchos::Ptr<Checkpoint>& chk) {
  if (io_checkpoint_) {
    EnsureSubfieldNames_();

    // loop over the components and dump them to the checkpoint file if possible
    int i = 0;
    for (CompositeVector::name_iterator compname=data_->begin();
         compname!=data_->end(); ++compname) {
      // get the MultiVector that should be dumped
      Teuchos::RCP<Epetra_MultiVector> v = data_->ViewComponent(*compname, false);

      // construct name for the field in the checkpoint
      std::vector<std::string> chkp_names(subfield_names_[i]);
      for (unsigned int j = 0; j!=subfield_names_[i].size(); ++j) {
        chkp_names[j] = fieldname_ + "." + *compname + "." + subfield_names_[i][j];
      }
      chk->WriteVector(*v, chkp_names);
      i++;
    }
  }
};


void Field_CompositeVector::ReadCellsFromCheckpoint_(std::string filename) {
  Teuchos::RCP<Amanzi::HDF5_MPI> file_input =
      Teuchos::rcp(new Amanzi::HDF5_MPI(*data_->comm(), filename));
  EnsureSubfieldNames_();

  int i = 0;
  for (CompositeVector::name_iterator compname=data_->begin();
       compname!=data_->end(); ++compname) {

    if (*compname == std::string("cell")) {
      // get the MultiVector that should be read
      Teuchos::RCP<Epetra_MultiVector> vec = data_->ViewComponent(*compname, false);

      // construct name for the field in the checkpoint file
      std::vector<std::string> chkp_names(subfield_names_[i]);
      for (unsigned int j = 0; j!=subfield_names_[i].size(); ++j) {
        chkp_names[j] = fieldname_ + "." + *compname + "." + subfield_names_[i][j];
      }
      for (unsigned int j = 0; j!=subfield_names_[i].size(); ++j) {
        file_input->readData(*(*vec)(j), chkp_names[j]);
      }
    }
    ++i;
  }
}

void Field_CompositeVector::ReadCheckpoint_(std::string filename) {
  Teuchos::RCP<Amanzi::HDF5_MPI> file_input =
      Teuchos::rcp(new Amanzi::HDF5_MPI(*data_->comm(), filename));
  ReadCheckpoint(file_input.ptr());
}


// modify methods
// -- set data from file
void Field_CompositeVector::ReadCheckpoint(const Teuchos::Ptr<HDF5_MPI>& file_input) {
  EnsureSubfieldNames_();

  // loop over the components and dump them to the checkpoint file if possible
  int i = 0;
  for (CompositeVector::name_iterator compname=data_->begin();
       compname!=data_->end(); ++compname) {
    // get the MultiVector that should be read
    Teuchos::RCP<Epetra_MultiVector> vec = data_->ViewComponent(*compname, false);

    // construct name for the field in the checkpoint file
    std::vector<std::string> chkp_names(subfield_names_[i]);
    for (unsigned int j = 0; j!=subfield_names_[i].size(); ++j) {
      chkp_names[j] = fieldname_ + "." + *compname + "." + subfield_names_[i][j];
    }
    for (unsigned int j = 0; j!=subfield_names_[i].size(); ++j) {
      file_input->readData(*(*vec)(j), chkp_names[j]);
    }
    i++;
  }
}


void Field_CompositeVector::EnsureSubfieldNames_() {
  // set default values for subfield names, ensuring they are unique
  if (subfield_names_.size() == 0) {
    subfield_names_.resize(data_->num_components());

    unsigned int i = 0;
    for (CompositeVector::name_iterator compname=data_->begin();
         compname!=data_->end(); ++compname) {
      subfield_names_[i].resize(data_->num_dofs(*compname));

      for (unsigned int j=0; j!=subfield_names_[i].size(); ++j) {
        std::stringstream s;
        s << j;
        subfield_names_[i][j] = s.str();
      }
      ++i;
    }
  } else {
    for (int i=0; i!=subfield_names_.size(); ++i) {
      for (int j=0; j!=subfield_names_[i].size(); ++j) {
        if (subfield_names_[i][j].length() == 0) {
          std::stringstream s;
          s << j;
          subfield_names_[i][j] = s.str();
        }
      }
    }
  }
};

} // namespace
