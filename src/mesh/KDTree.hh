/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov
*/

//! <MISSING_ONELINE_DOCSTRING>

#ifndef AMANZI_MESH_KD_TREE_HH_
#define AMANZI_MESH_KD_TREE_HH_

#include <memory>
#include <vector>

#include "Point.hh"

#include "Kokkos_Core.hpp"

namespace Amanzi {
namespace AmanziMesh {

// Interface class for nanoflann.
// Note that it does not own any data.
class PointCloud {
 public:
  PointCloud(){};
  ~PointCloud(){};

  // Interface for nanoflann
  // -- must return the number of points in the cloud
  inline size_t kdtree_get_point_count() const { return points_.extent(0); }

  // -- must return the i-th component of the n-th point
  inline double kdtree_get_pt(const size_t n, const size_t i) const
  {
    return points_(n)[i];
  }

  // -- must return optional bounding-box status:
  //    false to default to a standard bounding box computation loop.
  //    true if the box was already computed and returned in bb.
  //    Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3
  //    for point clouds)
  template <class BoundingBox>
  bool kdtree_get_bbox(BoundingBox& bb) const
  {
    return false;
  }

  void Init(const Kokkos::View<AmanziGeometry::Point*> points)
  {
    points_ = points;
  }

 private:
  Kokkos::View<AmanziGeometry::Point*> points_;
};


// At the moment, only one KDTree is used
typedef nanoflann::KDTreeSingleIndexAdaptor<
  nanoflann::L2_Adaptor<double, PointCloud>, PointCloud, -1>
  KDTree_L2Adaptor;

class KDTree {
 public:
  KDTree(){};
  ~KDTree(){};

  // main member function
  void Init(const Kokkos::View<AmanziGeometry::Point*> points)
  {
    int d = points(0).dim();
    cloud_.Init(points);
    tree_ = std::make_shared<KDTree_L2Adaptor>(
      d, cloud_, nanoflann::KDTreeSingleIndexAdaptorParams(10));
    tree_->buildIndex();
  }

  // find the first n points closest to the given point p
  std::vector<size_t> SearchNearest(const AmanziGeometry::Point& p,
                                    std::vector<double>& dist_sqr, int n = 1)
  {
    AMANZI_ASSERT(tree_ != NULL);

    double query[3];
    for (int i = 0; i < p.dim(); ++i) query[i] = p[i];

    std::vector<size_t> idx(n);
    dist_sqr.resize(n);

    int m = tree_->knnSearch(&query[0], n, &idx[0], &dist_sqr[0]);

    // the case of m < n
    idx.resize(m);
    dist_sqr.resize(m);

    return idx;
  }

  // find all points in the sphere of centered to the given point p
  std::vector<size_t>
  SearchInSphere(const AmanziGeometry::Point& p, std::vector<double>& dist_sqr,
                 double radius_sqr)
  {
    AMANZI_ASSERT(tree_ != NULL);

    double query[3];
    for (int i = 0; i < p.dim(); ++i) query[i] = p[i];

    std::vector<std::pair<size_t, double>> matches;
    nanoflann::SearchParams params;

    // params.sorted = false;
    int m = tree_->radiusSearch(&query[0], radius_sqr, matches, params);

    std::vector<size_t> idx(m);
    dist_sqr.resize(m);

    for (int i = 0; i < m; ++i) {
      idx[i] = matches[i].first;
      dist_sqr[i] = matches[i].second;
    }

    return idx;
  }

 private:
  PointCloud cloud_;
  std::shared_ptr<KDTree_L2Adaptor> tree_;
};

} // namespace AmanziMesh
} // namespace Amanzi

#endif
