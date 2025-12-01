/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*!
  Function applied to a mesh component, along with meta-data to store the values
  of this function in a ComposteVector.
*/

#include "errors.hh"
#include "VerboseObject.hh"
#include "MeshDefs.hh"
#include "CompositeVectorFunction.hh"

namespace Amanzi {
namespace Functions {

CompositeVectorFunction::CompositeVectorFunction(const Teuchos::RCP<const MeshFunction>& func,
                                                 const std::vector<std::string>& names,
                                                 bool dot_with_normal,
                                                 const std::string& spatial_dist_method,
						 bool must_cover)
  : func_(func),
    dot_with_normal_(dot_with_normal),
    spatial_dist_method_(spatial_dist_method)
{
  AMANZI_ASSERT(names.size() == func->size());

  if (dot_with_normal_) AMANZI_ASSERT(spatial_dist_method == "none");

  // Zip the two iterators, adding the resulting "tuple" to the container of pairs.
  std::vector<std::string>::const_iterator name = names.begin();
  for (MeshFunction::spec_iterator spec = func->begin(); spec != func->end(); ++spec) {
    cv_spec_list_.push_back(Teuchos::rcp(new CompositeVectorSpec(*name, *spec)));
    ++name;
  }

  if (must_cover) CheckCovers();
}

void
CompositeVectorFunction::Compute(double time,
                                 const Teuchos::Ptr<CompositeVector>& cv,
                                 const VerboseObject* vo)
{
  if (dot_with_normal_) {
    ComputeDotWithNormal_(time, cv, vo);
  } else if (spatial_dist_method_ != "none") {
    ComputeSpatiallyDistributed_(time, cv, vo);
  } else {
    Compute_(time, cv, vo);
  }
}


void
CompositeVectorFunction::Compute_(double time,
                                  const Teuchos::Ptr<CompositeVector>& cv,
                                  const VerboseObject* vo)
{
  cv->PutScalar(0.);

#ifdef ENSURE_INITIALIZED_CVFUNCS
  // ensure all components are touched
  std::map<std::string, bool> done;
  for (auto compname : *cv) {
    done[compname] = false;
  }
#endif

  // loop over the name/spec pair
  for (const Teuchos::RCP<CompositeVectorSpec>& cv_spec : cv_spec_list_) {
    const std::string& compname = cv_spec->first;
    const MeshFunction::Spec& spec = *cv_spec->second;
    Epetra_MultiVector& vec = *cv->ViewComponent(compname, false);

    AmanziMesh::Entity_kind kind = spec.first->second;
    if (vo && vo->os_OK(Teuchos::VERB_EXTREME)) {
      *vo->os() << "Writing function on component: " << compname << " and entity "
                << AmanziMesh::to_string(kind) << std::endl;
    }

    ComputeSpec_(spec, time, vec, vo);

#ifdef ENSURE_INITIALIZED_CVFUNCS
    done[compname] = true;
#endif
  }

#ifdef ENSURE_INITIALIZED_CVFUNCS
  for (auto compname : *cv) {
    if (!done[compname]) {
      Errors::Message message;
      message << "CV: component \"" << compname << "\" was not set";
      Exceptions::amanzi_throw(message);
    }
  }
#endif
}


void
CompositeVectorFunction::ComputeSpec_(const MeshFunction::Spec& spec,
                                      double time,
                                      Epetra_MultiVector& vec,
                                      const VerboseObject* vo)
{
  const AmanziMesh::Mesh& mesh = *func_->mesh();

  // create the input tuple
  int dim = mesh.getSpaceDimension();
  std::vector<double> args(1 + dim, 0.);
  args[0] = time;

  AmanziMesh::Entity_kind kind = spec.first->second;

  // loop over all regions in the spec
  for (const std::string& region : spec.first->first) {
    bool valid = mesh.isValidSetName(region, kind);
    if (vo && vo->os_OK(Teuchos::VERB_EXTREME)) {
      if (!valid) *vo->os() << "  region: " << region << " not valid!" << std::endl;
    }

    if (valid) {
      // get the indices of the domain.
      auto id_list = mesh.getSetEntities(region, kind, AmanziMesh::Parallel_kind::OWNED);
      const Epetra_BlockMap& map = vec.Map();

      if (vo && vo->os_OK(Teuchos::VERB_EXTREME)) {
        *vo->os() << "  region: " << region << " contains " << id_list.size() << " local entities"
                  << std::endl;
      }

      // loop over indices
      for (const auto& id : id_list) {
        // get the coordinate
        AmanziGeometry::Point xc = mesh.getCentroid(kind, id);
        for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];

        // evaluate the functions and stuff the result into the CV
        double* value = (*spec.second)(args);

        int g = map.FirstPointInElement(id);
        int ndofs = map.ElementSize(id);
        for (int n = 0; n < ndofs; ++n) {
          for (int i = 0; i != (*spec.second).size(); ++i) {
            vec[i][g + n] = value[i];
          }
        }
      }

    } else {
      Errors::Message message;
      message << "CV: unknown region \"" << region << "\"";
      Exceptions::amanzi_throw(message);
    }
  }
}


void
CompositeVectorFunction::ComputeDotWithNormal_(double time,
                                               const Teuchos::Ptr<CompositeVector>& cv,
                                               const VerboseObject* vo)
{
  AMANZI_ASSERT(cv->Map().NumComponents() == 1);                              // one comp
  AMANZI_ASSERT(cv->Map().HasComponent("face"));                              // is named face
  AMANZI_ASSERT(cv->Map().Location("face") == AmanziMesh::Entity_kind::FACE); // is on face
  AMANZI_ASSERT(cv->Map().NumVectors("face") == 1);                           // and is scalar

  // create a vector on faces of the appropriate dimension
  int dim = cv->Mesh()->getSpaceDimension();

  CompositeVectorSpace cvs;
  cvs.SetMesh(cv->Mesh());
  cvs.SetComponent("face", AmanziMesh::Entity_kind::FACE, dim);
  auto vel_vec = Teuchos::rcp(new CompositeVector(cvs));

  // Evaluate the velocity function
  Compute_(time, vel_vec.ptr(), vo);

  // CV's map may differ from the regular mesh map due to presense of fractures
  const auto& fmap = *cv->Map().Map("face", true);
  int nfaces_owned =
    cv->Mesh()->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);

  {
    Epetra_MultiVector& dat_f = *cv->ViewComponent("face");
    const Epetra_MultiVector& vel_f = *vel_vec->ViewComponent("face");

    // Dot the velocity with the normal
    // -- two branches: single flux per face, multiple fluxes
    int dir;
    AmanziGeometry::Point vel(dim);
    for (int f = 0; f != nfaces_owned; ++f) {
      for (int i = 0; i < dim; ++i) vel[i] = vel_f[i][f];

      int ndofs = fmap.ElementSize(f);
      int g = fmap.FirstPointInElement(f);
      if (ndofs == 1) {
        const AmanziGeometry::Point& normal = cv->Mesh()->getFaceNormal(f);
        dat_f[0][g] = vel * normal;
      } else {
        auto cells = cv->Mesh()->getFaceCells(f);

        for (int i = 0; i < ndofs; ++i) {
          const AmanziGeometry::Point& normal = cv->Mesh()->getFaceNormal(f, cells[i], &dir);
          dat_f[0][g + i] = (vel * normal) * dir;
        }
      }
    }
  }
}


void
CompositeVectorFunction::ComputeSpatiallyDistributed_(double time,
                                                      const Teuchos::Ptr<CompositeVector>& cv,
                                                      const VerboseObject* vo)
{
  cv->PutScalar(0.);

#ifdef ENSURE_INITIALIZED_CVFUNCS
  // ensure all components are touched
  std::map<std::string, bool> done;
  for (auto compname : *cv) {
    done[compname] = false;
  }
#endif

  // loop over the name/spec pair
  for (const Teuchos::RCP<CompositeVectorSpec>& cv_spec : cv_spec_list_) {
    const std::string& compname = cv_spec->first;
    const MeshFunction::Spec& spec = *cv_spec->second;
    Epetra_MultiVector& vec = *cv->ViewComponent(compname, false);

    AmanziMesh::Entity_kind kind = spec.first->second;
    if (vo && vo->os_OK(Teuchos::VERB_EXTREME)) {
      *vo->os() << "Writing function on component: " << compname << " and entity "
                << AmanziMesh::to_string(kind) << std::endl;
    }

    // key difference -- first copy construct and create a work vector, compute on that
    Epetra_MultiVector tmp_vec(vec);
    ComputeSpec_(spec, time, tmp_vec, vo);

    // now we know all regions are valid, etc, so we can compute the region
    // extent without error checking
    double total_volume = 0.;
    const AmanziMesh::Mesh& mesh = *func_->mesh();

    for (const std::string& region : spec.first->first) {
      auto id_list = mesh.getSetEntities(region, kind, AmanziMesh::Parallel_kind::OWNED);
      for (const auto& id : id_list) {
        total_volume += mesh.getExtent(kind, id);
      }
    }

    double l_total_volume(total_volume);
    mesh.getComm()->SumAll(&l_total_volume, &total_volume, 1);

    // spatial distribution method = volume:  uniform distribution across the volume
    vec.Update(1. / total_volume, tmp_vec, 1.);

#ifdef ENSURE_INITIALIZED_CVFUNCS
    done[compname] = true;
#endif
  }

#ifdef ENSURE_INITIALIZED_CVFUNCS
  for (auto compname : *cv) {
    if (!done[compname]) {
      Errors::Message message;
      message << "CV: component \"" << compname << "\" was not set";
      Exceptions::amanzi_throw(message);
    }
  }
#endif
}


void
CompositeVectorFunction::CheckCovers()
{
  // get a set of entity kinds and regions that will be used
  std::set<AmanziMesh::Entity_kind> kinds;
  std::map<AmanziMesh::Entity_kind, std::set<std::string>> regions_by_kind;
  for (const Teuchos::RCP<CompositeVectorSpec>& cv_spec : cv_spec_list_) {
    auto [regions, kind] = *cv_spec->second->first;
    kinds.insert(kind);
    for (const auto& region : regions) regions_by_kind[kind].insert(region);
  }

  // create a CV with component per entity kind, 1 dof
  CompositeVectorSpace cvs;
  cvs.SetMesh(func_->mesh());
  for (const AmanziMesh::Entity_kind& kind : kinds) {
    cvs.AddComponent(AmanziMesh::to_string(kind), kind, 1);
  }
  CompositeVector cv(cvs);
  cv.PutScalar(0.);

  // populate the cv indicating number of touches
  for (const auto& [kind, regions] : regions_by_kind) {
    std::string kind_as_string = AmanziMesh::to_string(kind);
    Epetra_MultiVector& vec = *cv.ViewComponent(kind_as_string, false);
    
    for (const std::string& region : regions) {
      const auto& region_ids = func_->mesh()->getSetEntities(region, kind, AmanziMesh::Parallel_kind::OWNED);
      for (const auto& e : region_ids) vec[0][e] += 1;
    }
  }

  // diagnostics
  for (const auto& comp : cv) {
    int count[2] = {0,0};
    AmanziMesh::Entity_GID where[2] = {-1, -1};
    AmanziMesh::Entity_kind kind = AmanziMesh::createEntityKind(comp);
    Epetra_MultiVector& vec = *cv.ViewComponent(comp, false);
    for (int i = 0; i != vec.MyLength(); ++i) {
      if (vec[0][i] == 0) {
	count[0]++;
	where[0] = func_->mesh()->getEntityGID(kind, i);
      } else if (vec[0][i] > 1) {
	count[1]++;
	where[1] = func_->mesh()->getEntityGID(kind, i);
      }
    }

    int count_g[2];
    AmanziMesh::Entity_GID where_g[2];
    cv.Comm()->SumAll(count, count_g, 2);

    if (where[0] == -1) where[0] = std::numeric_limits<AmanziMesh::Entity_GID>::max();
    if (where[1] == -1) where[1] = std::numeric_limits<AmanziMesh::Entity_GID>::max();
    cv.Comm()->MinAll(where, where_g, 2);
    if (where_g[0] == std::numeric_limits<AmanziMesh::Entity_GID>::max()) where_g[0] = -1;
    if (where_g[1] == std::numeric_limits<AmanziMesh::Entity_GID>::max()) where_g[1] = -1;

    if (count_g[0] > 0) {
      Errors::Message msg;
      msg << "CompositeVectorFunction component \"" << comp
	  << "\" is not covered by the provided functions, which are defined on "
	  << regions_by_kind[kind].size() << " regions:\n";
      for (const auto& r : regions_by_kind[kind]) {
	msg << "  \"" << r << "\"\n";
      }
      msg << "There are "
	  << count_g[0] << " of " << vec.GlobalLength()
	  << " total untouched entities, including GID " << where_g[0]
	  << ", and there are " << count_g[1] << " entities that are touched more than once, including GID "
	  << where_g[1] << ".";
      Exceptions::amanzi_throw(msg);
    }
  }    
}


} // namespace Functions
} // namespace Amanzi
