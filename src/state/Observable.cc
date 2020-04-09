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

double ObservableExtensiveSum(double a, double b, double vol) { return a + b; }
double ObservableIntensiveSum(double a, double b, double vol) { return a + b*vol; }
double ObservableMin(double a, double b, double vol) { return std::min(a,b); }
double ObservableMax(double a, double b, double vol) { return std::max(a,b); }


Observable::Observable(Teuchos::ParameterList& plist) :
    IOEvent(plist),
    count_(0),
    write_(false)
{
  // process the spec
  name_ = plist.name();
  variable_ = plist.get<std::string>("variable");
  region_ = plist.get<std::string>("region");
  delimiter_ = plist.get<std::string>("delimiter", ",");

  functional_ = plist.get<std::string>("functional");
  if (functional_ == "observation data: point" ||
      functional_ == "observation data: integral" ||
      functional_ == "observation data: average") {
    function_ = &ObservableIntensiveSum;
  } else if (functional_ == "observation data: extensive integral") {
    function_ = &ObservableExtensiveSum;
  } else if (functional_ == "observation data: minimum") {
    function_ = &ObservableMin;
  } else if (functional_ == "observation data: maximum") {
    function_ = &ObservableMax;
  } else {
    Errors::Message msg;
    msg << "Observable: unrecognized functional " << functional_;
    Exceptions::amanzi_throw(msg);
  }
  
  // entity of region
  location_ = plist.get<std::string>("location name", "cell");

  // hack to orient flux to outward-normal along a boundary only
  flux_normalize_ = plist.get<bool>("direction normalized flux", false);
  if (flux_normalize_ && plist.isParameter("direction normalized flux direction")) {
    Teuchos::Array<double> direction =
        plist.get<Teuchos::Array<double> >("direction normalized flux direction");
    if (direction.size() == 2) {
      double norm = std::sqrt(std::pow(direction[0],2) + std::pow(direction[1],2));
      direction_ = Teuchos::rcp(new AmanziGeometry::Point(direction[0]/norm, direction[1]/norm));
    } else if (direction.size() == 3) {
      double norm = std::sqrt(std::pow(direction[0],2) + std::pow(direction[1],2)
              + std::pow(direction[2],2));
      direction_ = Teuchos::rcp(new AmanziGeometry::Point(direction[0]/norm,
              direction[1]/norm, direction[2]/norm));
    } else {
      Errors::Message msg;
      msg << "Observable: \"direction normalized flux direction\" cannot have dimension " << (int) direction.size() << ", must be 2 or 3.";
      Exceptions::amanzi_throw(msg);
    }
  }

  // write mode
  interval_ = plist.get<int>("write interval", 0);
  filenamebase_ = plist.get<std::string>("observation output filename");
}

void Observable::Update(const State& S,
                        Amanzi::ObservationData::DataQuadruple& data) {
  Update_(S, data);

  // open the file if I am the writing process and it isn't already open
  if (write_ && !out_.get()) {
    std::string safename(name_);
    std::replace(safename.begin(), safename.end(), ' ', '_');
    std::replace(safename.begin(), safename.end(), ':', '_');
    std::stringstream filename;
    filename << filenamebase_ << "_" << safename;
    AMANZI_ASSERT(boost::filesystem::portable_file_name(filenamebase_));
    out_ = Teuchos::rcp(new std::ofstream(filenamebase_.c_str()));
    WriteHeader_();
  }

  if (out_.get()) {
    if (data.is_valid) {
      *out_ << data.time << delimiter_ << " " << data.value << std::endl;
    } else {
      *out_ << data.time << delimiter_ << " " << "NaN" << std::endl;
    }

    if (count_ % interval_ == 0) out_->flush();
  }
  ++count_;
}

void Observable::Flush() {
  if (out_.get()) out_->flush();
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
    out_->precision(16);
    *out_ << std::scientific;
  }
}

void Observable::Update_(const State& S,
                         Amanzi::ObservationData::DataQuadruple& data) {
  data.time = S.time();

  Teuchos::RCP<const Field> field = S.GetField(variable_);

  if (field->type() == CONSTANT_SCALAR) {
    // only write on MPI_COMM_WORLD rank 0
    auto comm = Amanzi::getDefaultComm();
    write_ = comm->MyPID() == 0;

    // scalars, just return the value
    data.value = *field->GetScalarData();
    data.is_valid = true;

    
  } else if (field->type() == COMPOSITE_VECTOR_FIELD) {
    // vector field
    Teuchos::RCP<const CompositeVector> vec = field->GetFieldData();
    AMANZI_ASSERT(vec->HasComponent(location_));

    // only write on field's comm's rank 0
    write_ = vec->Mesh()->get_comm()->MyPID() == 0;

    // get the region
    AmanziMesh::Entity_kind entity = vec->Location(location_);
    AmanziMesh::Entity_ID_List ids;
    vec->Mesh()->get_set_entities(region_, entity, AmanziMesh::Parallel_type::OWNED, &ids);

    double value(0.);
    if (functional_ == "observation data: minimum") {
      value = 1.e20;
    } else if (functional_ == "observation data: maximum") {
      value = -1.e20;
    }

    double volume(0.);
    const Epetra_MultiVector& subvec = *vec->ViewComponent(location_, false);

    if (entity == AmanziMesh::CELL) {
      for (AmanziMesh::Entity_ID_List::const_iterator id=ids.begin();
           id!=ids.end(); ++id) {
        double vol = vec->Mesh()->cell_volume(*id);
        value = (*function_)(value, subvec[0][*id], vol);
        volume += vol;
      }
    } else if (entity == AmanziMesh::FACE) {
      for (AmanziMesh::Entity_ID_List::const_iterator id=ids.begin();
           id!=ids.end(); ++id) {
        double vol = vec->Mesh()->face_area(*id);

        // hack to orient flux to outward-normal along a boundary only
        double sign = 1;
        if (flux_normalize_) {
          if (direction_.get()) {
            // normalize to the provided vector
            AmanziGeometry::Point normal = vec->Mesh()->face_normal(*id);
            sign = (normal * (*direction_)) / AmanziGeometry::norm(normal);
            
          } else {
            // normalize to outward normal
            AmanziMesh::Entity_ID_List cells;
            vec->Mesh()->face_get_cells(*id, AmanziMesh::Parallel_type::ALL, &cells);
            AmanziMesh::Entity_ID_List faces;
            std::vector<int> dirs;
            vec->Mesh()->cell_get_faces_and_dirs(cells[0], &faces, &dirs);
            int i = std::find(faces.begin(), faces.end(), *id) - faces.begin();
            sign = dirs[i];
            
          }
        }

        value = (*function_)(value, sign*subvec[0][*id], vol);
        volume += std::abs(vol);
      }
    } else if (entity == AmanziMesh::NODE) {
      for (AmanziMesh::Entity_ID_List::const_iterator id=ids.begin();
           id!=ids.end(); ++id) {
        double vol = 1.0;
        value = (*function_)(value, subvec[0][*id], vol);
        volume += vol;
      }
    }

    // syncronize the result across processors
    if (functional_ == "observation data: point" ||
        functional_ == "observation data: integral" ||
        functional_ == "observation data: extensive integral") {
      double local[2], global[2];
      local[0] = value; local[1] = volume;
      S.GetMesh()->get_comm()->SumAll(local, global, 2);

      if (global[1] > 0) {
        if (functional_ == "observation data: point") {
          data.value = global[0] / global[1];
          data.is_valid = true;
        } else if (functional_ == "observation data: integral" ||
                   functional_ == "observation data: extensive integral") {
          data.value = global[0];
          data.is_valid = true;
        }
      } else {
        data.value = 0.;
        data.is_valid = false;
      }
    } else if (functional_ == "observation data: minimum") {
      double global;
      S.GetMesh()->get_comm()->MinAll(&value, &global, 1);
      data.value = global;
      data.is_valid = true;
    } else if (functional_ == "observation data: maximum") {
      double global;
      S.GetMesh()->get_comm()->MaxAll(&value, &global, 1);
      data.value = global;
      data.is_valid = true;
    } else {
      data.value = 0.;
      data.is_valid = false;
    }
  } else {
    write_ = false;
    data.value = 0.;
    data.is_valid = false;
  }
}


} // namespace


