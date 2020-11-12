/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
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

#include "Key.hh"
#include "errors.hh"
#include "Mesh.hh"
#include "State.hh"
#include "Field.hh"
#include "FieldEvaluator.hh"

#include "Observable.hh"

namespace Amanzi {

namespace Impl {

double ObservableExtensiveSum(double a, double b, double vol) { return a + b; }
double ObservableIntensiveSum(double a, double b, double vol) { return a + b*vol; }
double ObservableMin(double a, double b, double vol) { return std::min(a,b); }
double ObservableMax(double a, double b, double vol) { return std::max(a,b); }

} // namespace Impl

const double Observable::nan = std::numeric_limits<double>::quiet_NaN();

Observable::Observable(Teuchos::ParameterList& plist,
                       const Teuchos::Ptr<State>& S)
  : old_time_(nan)
{
  // process the spec
  name_ = Keys::cleanPListName(plist.name());
  variable_ = plist.get<std::string>("variable");
  region_ = plist.get<std::string>("region");
  location_ = plist.get<std::string>("location name", "cell");
  num_vectors_ = plist.get<int>("number of vectors", 1);
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
    msg << "Observable " << name_ << ": \"direction normalized flux direction\" may only be used with location \"face\"";
    Exceptions::amanzi_throw(msg);
  }
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
      msg << "Observable: \"direction normalized flux direction\" cannot have dimension "
          << (int) direction.size() << ", must be 2 or 3.";
      Exceptions::amanzi_throw(msg);
    }
  }

  // ensure the field exists, and require on location with num_vectors
  S->RequireField(variable_)
    ->SetMesh(S->GetMesh(Keys::getDomain(variable_)))
    ->AddComponent(location_, AmanziMesh::entity_kind(location_), num_vectors_);

  // try to get a field evaluator -- note, this isn't especially well supported
  // by State
  if (S->HasFieldEvaluator(variable_)) {
    has_eval_ = true;
  } else if (S->FEList().isSublist(variable_)) {
    has_eval_ = true;
    S->RequireFieldEvaluator(variable_);
  } else {
    has_eval_ = false;
  }
}

void Observable::Update(const Teuchos::Ptr<State>& S,
                        std::vector<double>& data, int start_loc)
{
  // deal with the time integrated case for the first observation
  if (time_integrated_ && std::isnan(old_time_)) {
    for (int i=0; i!=num_vectors_; ++i) data[start_loc+i] = 0.;
    old_time_ = S->time();
    return;
  }

  // update the variable
  if (has_eval_) {
    S->GetFieldEvaluator(variable_)->HasFieldChanged(S, "observation");
  }

  Teuchos::RCP<const Field> field = S->GetField(variable_);
  if (field->type() == CONSTANT_SCALAR) {
    // scalars, just return the value
    data[start_loc] = *field->GetScalarData();

  } else if (field->type() == COMPOSITE_VECTOR_FIELD) {
    // vector field
    Teuchos::RCP<const CompositeVector> vec = field->GetFieldData();
    AMANZI_ASSERT(vec->HasComponent(location_));

    // get the region
    AmanziMesh::Entity_kind entity = AmanziMesh::entity_kind(location_);
    AmanziMesh::Entity_ID_List ids;
    vec->Mesh()->get_set_entities(region_, entity, AmanziMesh::Parallel_type::OWNED, &ids);

    std::vector<double> value;
    if (functional_ == "minimum") {
      value.resize(num_vectors_ + 1, 1.e20);
    } else if (functional_ == "maximum") {
      value.resize(num_vectors_ + 1, -1.e20);
    } else {
      value.resize(num_vectors_ + 1, 0.);
    }

    const Epetra_MultiVector& subvec = *vec->ViewComponent(location_, false);

    if (entity == AmanziMesh::CELL) {
      for (AmanziMesh::Entity_ID_List::const_iterator id=ids.begin();
           id!=ids.end(); ++id) {
        double vol = vec->Mesh()->cell_volume(*id);
        for (int i=0; i!=num_vectors_; ++i) {
          value[i] = (*function_)(value[i], subvec[i][*id], vol);
        }
        value[num_vectors_] += vol;
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

        for (int i=0; i!=num_vectors_; ++i) {
          value[i] = (*function_)(value[i], sign*subvec[i][*id], vol);
        }
        value[num_vectors_] += std::abs(vol);
      }
    } else if (entity == AmanziMesh::NODE) {
      for (AmanziMesh::Entity_ID_List::const_iterator id=ids.begin();
           id!=ids.end(); ++id) {
        double vol = 1.0;

        for (int i=0; i!=num_vectors_; ++i) {
          value[i] = (*function_)(value[i], subvec[i][*id], vol);
        }
        value[num_vectors_] += vol;
      }
    }

    // syncronize the result across processors
    if (functional_ == "point" ||
        functional_ == "integral" ||
        functional_ == "average" ||
        functional_ == "extensive integral") {
      std::vector<double> global_value(value);
      S->GetMesh()->get_comm()->SumAll(value.data(), global_value.data(), value.size());

      if (global_value[num_vectors_] > 0) {
        if (functional_ == "point" ||
            functional_ == "average") {
          for (int i=0; i!=num_vectors_; ++i) {
            value[i] = global_value[i] / global_value[num_vectors_];
          }
        } else if (functional_ == "integral" ||
                   functional_ == "extensive integral") {
          for (int i=0; i!=num_vectors_; ++i) {
            value[i] = global_value[i];
          }
        }
      } else {
        for (int i=0; i!=num_vectors_; ++i) {
          value[i] = nan;
        }
      }
    } else if (functional_ == "minimum") {
      std::vector<double> global_value(value);
      S->GetMesh()->get_comm()->MinAll(value.data(), global_value.data(), value.size()-1);
      value = global_value;
    } else if (functional_ == "maximum") {
      std::vector<double> global_value(value);
      S->GetMesh()->get_comm()->MaxAll(value.data(), global_value.data(), value.size()-1);
      value = global_value;
    }

    for (int i=0; i!=num_vectors_; ++i) {
      data[start_loc+i] = value[i];
    }

  } else {
    for (int i=0; i!=num_vectors_; ++i) {
      data[start_loc+i] = nan;
    }
  }

  if (time_integrated_) {
    double dt = S->time() - old_time_;
    old_time_ = S->time();
    for (int i=0; i!=num_vectors_; ++i) {
      data[start_loc+i] *= dt;
    }
  }
}


} // namespace


