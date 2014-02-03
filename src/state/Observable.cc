 /* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
Amanzi

License:
Author: Markus Berndt
        Ethan Coon (ecoon@lanl.gov)

Observable data object

------------------------------------------------------------------------- */

#include <string>
#include <algorithm>

#include <boost/filesystem/operations.hpp>

#include "errors.hh"
#include "Mesh.hh"
#include "State.hh"
#include "Field.hh"

#include "Observable.hh"


namespace Amanzi {

Observable::Observable(Teuchos::ParameterList& plist, Epetra_MpiComm *comm) :
    IOEvent(plist, comm),
    count_(0) {
  // process the spec
  name_ = plist.name();
  variable_ = plist.get<std::string>("Variable");
  region_ = plist.get<std::string>("Region");
  functional_ = plist.get<std::string>("Functional");

  // entity of region
  location_ = plist.get<std::string>("Location Name", "cell");

  // write mode
  interval_ = plist.get<int>("Write Interval", 0);
  write_ = interval_ > 0;

  if (write_) {
    filenamebase_ = plist.get<std::string>("Observation Output Filename");

    // open file only on process 0
    if (!comm->MyPID()) {
      std::string safename(name_);
      std::replace(safename.begin(), safename.end(), ' ', '_');
      std::replace(safename.begin(), safename.end(), ':', '_');
      std::stringstream filename;
      filename << filenamebase_ << "_" << safename;
      ASSERT(boost::filesystem::portable_file_name(filenamebase_));
      out_ = Teuchos::rcp(new std::ofstream(filenamebase_.c_str()));
    }
  }
}

void Observable::Update(const State& S,
                        Amanzi::ObservationData::DataTriple& data) {
  if (count_ == 0) WriteHeader_();

  ++count_;
  Update_(S, data);

  if (out_.get()) {
    if (data.is_valid) {
      *out_ << data.time << ", " << data.value << std::endl;
    } else {
      *out_ << data.time << ", " << "NaN" << std::endl;
    }

    if (count_ % interval_ == 0) out_->flush();
  }
}

void Observable::WriteHeader_() {
  if (out_.get()) {
    *out_ << "# Observation Name: " << name_ << std::endl;
    *out_ << "# Region: " << region_ << std::endl;
    *out_ << "# Functional: " << functional_ << std::endl;
    *out_ << "# Variable: " << variable_ << std::endl;
    *out_ << "# ==========================================================="
          << std::endl;
    *out_ << "#" << std::endl;
  }
}

void Observable::Update_(const State& S,
                          Amanzi::ObservationData::DataTriple& data) {
  data.time = S.time();

  Teuchos::RCP<const Field> field = S.GetField(variable_);

  if (field->type() == CONSTANT_SCALAR) {
    // scalars, just return the value
    data.value = *field->GetScalarData();
    data.is_valid = true;

  } else if (field->type() == COMPOSITE_VECTOR_FIELD) {
    // vector field
    Teuchos::RCP<const CompositeVector> vec = field->GetFieldData();
    ASSERT(vec->HasComponent(location_));

    AmanziMesh::Entity_kind entity = vec->Location(location_);
    AmanziMesh::Entity_ID_List ids;
    vec->Mesh()->get_set_entities(region_, entity, AmanziMesh::OWNED, &ids);

    double value(0.);
    double volume(0.);
    const Epetra_MultiVector& subvec = *vec->ViewComponent(location_, false);

    if (entity == AmanziMesh::CELL) {
      for (AmanziMesh::Entity_ID_List::const_iterator id=ids.begin();
           id!=ids.end(); ++id) {
        double vol = vec->Mesh()->cell_volume(*id);
        value += subvec[0][*id] * vol;
        volume += vol;
      }
    } else if (entity == AmanziMesh::FACE) {
      for (AmanziMesh::Entity_ID_List::const_iterator id=ids.begin();
           id!=ids.end(); ++id) {
        double vol = vec->Mesh()->face_area(*id);
        value += subvec[0][*id] * vol;
        volume += vol;
      }
    } else if (entity == AmanziMesh::NODE) {
      for (AmanziMesh::Entity_ID_List::const_iterator id=ids.begin();
           id!=ids.end(); ++id) {
        double vol = 1.0;
        value += subvec[0][*id] * vol;
        volume += vol;
      }
    }

    if (volume > 0) {
      if (functional_ == "Observation Data: Point") {
        data.value = value / volume;
        data.is_valid = true;
      } else if (functional_ == "Observation Data: Integral") {
        data.value = value;
        data.is_valid = true;
      } else {
        std::stringstream message;
        message << "Observable: unrecognized functional " << functional_;
        Errors::Message m(message.str());
        Exceptions::amanzi_throw(m);
      }
    } else {
      data.value = 0.;
      data.is_valid = false;
    }
  } else {
    data.value = 0.;
    data.is_valid = false;
  }
}


} // namespace


