/*
  Copyright 2010-202x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*

This class calculates the actual observation value.

*/

#include <cmath>
#include <string>
#include <algorithm>

#include "mpi.h"
#include "Key.hh"
#include "errors.hh"
#include "Key.hh"
#include "Mesh.hh"

// Amanzi::State
#include "Evaluator.hh"
#include "Observable.hh"
#include "State.hh"

namespace Amanzi {

namespace Impl {

double ObservableExtensiveSum(double a, double b, double vol) { return a + b; }
double ObservableIntensiveSum(double a, double b, double vol) { return a + b*vol; }
double ObservableMin(double a, double b, double vol) { return std::min(a,b); }
double ObservableMax(double a, double b, double vol) { return std::max(a,b); }

} // namespace Impl

const double Observable::nan = std::numeric_limits<double>::quiet_NaN();

Observable::Observable(Teuchos::ParameterList& plist)
  : comm_(Teuchos::null),
    old_time_(nan),
    has_eval_(false),
    has_data_(false)
{
  // process the spec
  name_ = Keys::cleanPListName(plist.name());
  variable_ = plist.get<std::string>("variable");
  region_ = plist.get<std::string>("region");
  location_ = plist.get<std::string>("location name", "cell");
  tag_ = Tag(plist.get<std::string>("tag", Tags::DEFAULT.get()));

  // Note: -1 here means either take it from the physics if possible, or if
  // this variable is not in the physics, instead will default to 1.
  num_vectors_ = plist.get<int>("number of vectors", -1);
  dof_ = plist.get<int>("degree of freedom", -1);

  if (num_vectors_ > 0 && num_vectors_ < dof_) {
    Errors::Message msg;
    msg << "Observable \"" << name_ << "\": inconsistent request of degree of freedom "
        << dof_ << " for a vector with only " << num_vectors_ << " degrees of freedom.";
    Exceptions::amanzi_throw(msg);
  }

  time_integrated_ = plist.get<bool>("time integrated", false);

  functional_ = plist.get<std::string>("functional");
  if (functional_ == "point" ||
      functional_ == "integral" ||
      functional_ == "average") {
    function_ = &Impl::ObservableIntensiveSum;
  } else if (functional_ == "extensive integral") {
    function_ = &Impl::ObservableExtensiveSum;
  } else if (functional_ == "minimum") {
    function_ = &Impl::ObservableMin;
  } else if (functional_ == "maximum") {
    function_ = &Impl::ObservableMax;
  } else {
    Errors::Message msg;
    msg << "Observable: unrecognized functional " << functional_;
    Exceptions::amanzi_throw(msg);
  }

  // hack to orient flux to outward-normal along a boundary only
  flux_normalize_ = plist.get<bool>("direction normalized flux", false);
  if (flux_normalize_ && location_ != "face") {
    Errors::Message msg;
    msg << "Observable \"" << name_ << "\": \"direction normalized flux direction\" may only be used with location \"face\"";
    Exceptions::amanzi_throw(msg);
  }
  if (flux_normalize_) {
    if (plist.isParameter("direction normalized flux direction")) {
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
        msg << "Observable \"" << name_ << "\": \"direction normalized flux direction\" cannot have dimension "
            << (int) direction.size() << ", must be 2 or 3.";
        Exceptions::amanzi_throw(msg);
      }
    } else if (plist.isParameter("direction normalized flux relative to region")) {
      flux_normalize_region_ = plist.get<std::string>("direction normalized flux relative to region");
    }
  }
}


void Observable::Setup(const Teuchos::Ptr<State>& S)
{
  // Do we participate in this communicator?
  if (comm_ == Teuchos::null) return;

  // We may be in the comm, but not have this variable.
  Key domain = Keys::getDomain(variable_);
  if (!S->HasMesh(domain)) return;

  // If we have gotten to this point, we must have data
  has_data_ = true;

  // does the observed quantity have an evaluator?  Or can we make one? Note
  // that a non-evaluator based observation must already have been created by
  // PKs by now, because PK->Setup() has already been run.
  if (!S->HasRecord(variable_, tag_)) {
    // not yet created, require evaluator
    S->RequireEvaluator(variable_, tag_);
    has_eval_ = true;
  } else {
    // does it have an evaluator or can we make one?
    if (S->HasEvaluatorList(variable_)) {
      S->RequireEvaluator(variable_, tag_);
      has_eval_ = true;
    } else {
      has_eval_ = false;
    }
  }

  // try to set requirements on the field, if they are not already set
  if (!S->HasRecord(variable_, tag_)) {
    // require the field
    auto& cvs = S->Require<CompositeVector, CompositeVectorSpace>(variable_, tag_);

    // we have to set the mesh now -- assume it is provided by the domain
    cvs.SetMesh(S->GetMesh(Keys::getDomain(variable_)));

    // was num_vectors set?  if not use default of 1
    if (num_vectors_ < 0) num_vectors_ = 1;
    if (num_vectors_ < dof_) {
      Errors::Message msg;
      msg << "Observable \"" << name_ << "\": inconsistent request of degree of freedom "
          << dof_ << " for a vector with only " << num_vectors_ << " degrees of freedom.";
      Exceptions::amanzi_throw(msg);
    }

    // require the component on location_ with num_vectors_
    cvs.AddComponent(location_, AmanziMesh::entity_kind(location_), num_vectors_);
  }

  // communicate so that all ranks know number of vectors
  int num_vectors_local = num_vectors_;
  comm_->MaxAll(&num_vectors_local, &num_vectors_, 1);
}


void Observable::FinalizeStructure(const Teuchos::Ptr<State>& S)
{
  // if we don't have a communicator, we can't participate in this
  if (comm_ == Teuchos::null) return;

  // one last check that the structure is all set up and consistent
  if (has_data_ && num_vectors_ < 0) {
    const auto& field = S->Get<CompositeVector>(variable_, tag_);

    if (!field.HasComponent(location_)) {
      Errors::Message msg;
      msg << "Observable: \"" << name_ << "\" uses variable \""
          << variable_ << "\" but this field does not have the observed component \""
          << location_ << "\"";
      Exceptions::amanzi_throw(msg);
    } else {
      num_vectors_ = field.NumVectors(location_);
    }

    if (num_vectors_ < dof_) {
      Errors::Message msg;
      msg << "Observable \"" << name_ << "\": inconsistent request of degree of freedom "
          << dof_ << " for a vector with only " << num_vectors_ << " degrees of freedom.";
      Exceptions::amanzi_throw(msg);
    }

    // communicate so that all ranks know number of vectors
    int num_vectors_local = num_vectors_;
    comm_->MaxAll(&num_vectors_local, &num_vectors_, 1);
  }

  // must communicate the number of vectors so that all in comm have the right
  // size and get_num_vectors() is valid
  int num_vectors_l(num_vectors_);
  comm_->MaxAll(&num_vectors_l, &num_vectors_, 1);
}


void Observable::Update(const Teuchos::Ptr<State>& S,
                        std::vector<double>& data, int start_loc)
{
  // if we don't have a communicator, we do not participate, so leave the value at NaN
  if (comm_ == Teuchos::null) return;

  // deal with the time integrated case for the first observation
  if (time_integrated_ && std::isnan(old_time_)) {
    for (int i=0; i!=get_num_vectors(); ++i) data[start_loc+i] = 0.;
    old_time_ = S->get_time();
    return;
  }

  // from this point forward, has_data_ may be true or false, but we know that
  // the comm is valid so we must participate in communication!  This deals
  // with the case of an observation not on comm's rank 0, but must be
  // communicated to rank 0 to write.  We do know that get_num_vectors() is valid though.
  AMANZI_ASSERT(get_num_vectors() >= 0);
  std::vector<double> value;
  if (functional_ == "minimum") {
    value.resize(get_num_vectors() + 1, 1.e20);
  } else if (functional_ == "maximum") {
    value.resize(get_num_vectors() + 1, -1.e20);
  } else {
    value.resize(get_num_vectors() + 1, 0.);
  }

  // update the variable
  if (has_eval_)
    S->GetEvaluator(variable_, tag_).Update(*S, "observation");

  bool has_record = S->HasRecord(variable_, tag_);
  if (has_record && S->GetRecord(variable_, tag_).ValidType<double>()) {
    // scalars, just return the value
    value[0] = S->GetRecord(variable_, tag_).Get<double>();
    value[1] = 1;

  } else if (has_record && S->GetRecord(variable_, tag_).ValidType<CompositeVector>()) {
    // vector field
    const auto& vec = S->GetRecord(variable_, tag_).Get<CompositeVector>();
    AMANZI_ASSERT(vec.HasComponent(location_));

    // get the region
    AmanziMesh::Entity_kind entity = AmanziMesh::entity_kind(location_);
    AmanziMesh::Entity_ID_List ids;

    // get the vector component
    vec.Mesh()->get_set_entities(region_, entity, AmanziMesh::Parallel_type::OWNED, &ids);
    const Epetra_MultiVector& subvec = *vec.ViewComponent(location_, false);

    if (entity == AmanziMesh::CELL) {
      for (auto id : ids) {
        double vol = vec.Mesh()->cell_volume(id);

        if (dof_ < 0) {
          for (int i=0; i!=get_num_vectors(); ++i) {
            value[i] = (*function_)(value[i], subvec[i][id], vol);
          }
        } else {
          value[0] = (*function_)(value[0], subvec[dof_][id], vol);
        }
        value[get_num_vectors()] += vol;
      }
    } else if (entity == AmanziMesh::FACE) {
      for (auto id : ids) {
        double vol = vec.Mesh()->face_area(id);

        // hack to orient flux to outward-normal along a boundary only
        double sign = 1;
        if (flux_normalize_) {
          if (direction_.get()) {
            // normalize to the provided vector
            AmanziGeometry::Point normal = vec.Mesh()->face_normal(id);
            sign = (normal * (*direction_)) / AmanziGeometry::norm(normal);

          } else if (!flux_normalize_region_.empty()) {
            // normalize to outward normal relative to a volumetric region
            AmanziMesh::Entity_ID_List vol_cells;
            vec.Mesh()->get_set_entities(flux_normalize_region_,
                    AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::ALL, &vol_cells);

            // which cell of the face is "inside" the volume
            AmanziMesh::Entity_ID_List cells;
            vec.Mesh()->face_get_cells(id, AmanziMesh::Parallel_type::ALL, &cells);
            AmanziMesh::Entity_ID c = -1;
            for (const auto& cc : cells) {
              if (std::find(vol_cells.begin(), vol_cells.end(), cc) != vol_cells.end()) {
                c = cc;
                break;
              }
            }
            if (c < 0) {
              Errors::Message msg;
              msg << "Observeable on face region \"" << region_
                  << "\" flux normalized relative to volumetric region \""
                  << flux_normalize_region_ << "\" but face "
                  << vec.Mesh()->face_map(true).GID(id)
                  << " does not border the volume region.";
              Exceptions::amanzi_throw(msg);
            }

            // normalize with respect to that cell's direction
            AmanziMesh::Entity_ID_List faces;
            std::vector<int> dirs;

            vec.Mesh()->cell_get_faces_and_dirs(c, &faces, &dirs);
            int i = std::find(faces.begin(), faces.end(), id) - faces.begin();

            sign = dirs[i];

          } else {
            // normalize to outward normal
            AmanziMesh::Entity_ID_List cells;
            vec.Mesh()->face_get_cells(id, AmanziMesh::Parallel_type::ALL, &cells);
            AmanziMesh::Entity_ID_List faces;
            std::vector<int> dirs;
            vec.Mesh()->cell_get_faces_and_dirs(cells[0], &faces, &dirs);
            int i = std::find(faces.begin(), faces.end(), id) - faces.begin();
            sign = dirs[i];
          }
        }

        if (dof_ < 0) {
          for (int i=0; i!=get_num_vectors(); ++i) {
            value[i] = (*function_)(value[i], sign*subvec[i][id], vol);
          }
        } else {
          value[0] = (*function_)(value[0], sign*subvec[dof_][id], vol);
        }
        value[get_num_vectors()] += std::abs(vol);
      }
    } else if (entity == AmanziMesh::NODE) {
      for (auto id : ids) {
        double vol = 1.0;

        if (dof_ < 0) {
          for (int i=0; i!=get_num_vectors(); ++i) {
            value[i] = (*function_)(value[i], subvec[i][id], vol);
          }
        } else {
          value[0] = (*function_)(value[0], subvec[dof_][id], vol);
        }
        value[get_num_vectors()] += vol;
      }
    }
  }

  // syncronize the result across all processes on the provided comm
  if (functional_ == "point" ||
      functional_ == "integral" ||
      functional_ == "average" ||
      functional_ == "extensive integral") {
    std::vector<double> global_value(value);
    comm_->SumAll(value.data(), global_value.data(), value.size());

    if (global_value[get_num_vectors()] > 0) {
      if (functional_ == "point" ||
          functional_ == "average") {
        for (int i=0; i!=get_num_vectors(); ++i) {
          value[i] = global_value[i] / global_value[get_num_vectors()];
        }
      } else if (functional_ == "integral" ||
                 functional_ == "extensive integral") {
        for (int i=0; i!=get_num_vectors(); ++i) {
          value[i] = global_value[i];
        }
      }
    } else {
      for (int i=0; i!=get_num_vectors(); ++i) {
        value[i] = nan;
      }
    }
  } else if (functional_ == "minimum") {
    std::vector<double> global_value(value);
    comm_->MinAll(value.data(), global_value.data(), value.size()-1);
    value = global_value;
  } else if (functional_ == "maximum") {
    std::vector<double> global_value(value);
    comm_->MaxAll(value.data(), global_value.data(), value.size()-1);
    value = global_value;
  } else {
    AMANZI_ASSERT(false);
  }

  // copy back from value into the data array
  for (int i=0; i!=get_num_vectors(); ++i) {
    data[start_loc+i] = value[i];
  }

  // factor of dt for time integration
  if (time_integrated_) {
    double dt = S->get_time() - old_time_;
    old_time_ = S->get_time();
    for (int i=0; i!=get_num_vectors(); ++i) {
      data[start_loc+i] *= dt;
    }
  }
}

} // namespace


