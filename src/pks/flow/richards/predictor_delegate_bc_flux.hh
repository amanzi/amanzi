/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Delegate for modifying the predictor in the case of infiltration into dry soil.

  License: BSD
  Authors: Ethan Coon (ATS version) (ecoon@lanl.gov)
*/

#ifndef PREDICTOR_DELEGATE_BC_FLUX_
#define PREDICTOR_DELEGATE_BC_FLUX_

#include "Mesh.hh"
#include "state.hh"

#include "matrix_mfd.hh"
#include "wrm_partition.hh"

namespace Amanzi {
namespace Flow {

class PredictorDelegateBCFlux {

 public:
  PredictorDelegateBCFlux(const Teuchos::RCP<const State>& S_next,
                          const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                          const Teuchos::RCP<Operators::MatrixMFD>& matrix,
                          const Teuchos::RCP<FlowRelations::WRMPartition>& wrms,
                          std::vector<Operators::Matrix_bc>* bc_markers,
                          std::vector<double>* bc_values) :
      S_next_(S_next),
      mesh_(mesh),
      matrix_(matrix),
      wrms_(wrms),
      bc_markers_(bc_markers),
      bc_values_(bc_values)
  {}

  bool modify_predictor(double h, Teuchos::RCP<TreeVector> u);

 protected:

  class FluxBCFunctor {
   public:
    FluxBCFunctor(const Teuchos::RCP< std::vector<double> > Aff,
                  const Teuchos::RCP< std::vector<double> > lambda,
                  int face_index,
                  double cell_p,
                  double bc_flux,
                  double g_flux,
                  int dir,
                  double patm,
                  const Teuchos::RCP<FlowRelations::WRM>& wrm) :
        Aff_(Aff), lambda_(lambda), face_index_(face_index),
        cell_p_(cell_p), bc_flux_(bc_flux), g_flux_(g_flux),
        wrm_(wrm), dir_(dir), patm_(patm) {
    }

    double operator()(double face_p) {
      (*lambda_)[face_index_] = face_p;
      double Krel = wrm_->k_relative(patm_ - face_p);
      std::cout << "Fluxes: " << flux_() << ", " << g_flux_*Krel << ", " << bc_flux_ << std::endl;
      return flux_() + g_flux_*Krel - bc_flux_;
    }

   protected:

    double flux_() {
      double s = 0.;
      double Krel = wrm_->k_relative(patm_ - (*lambda_)[face_index_]);

      for (int n=0; n!=lambda_->size(); ++n)
        s += (*Aff_)[n] * Krel * (cell_p_ - (*lambda_)[n]);
      return s * dir_;
    }

   protected:
    Teuchos::RCP< std::vector<double> > Aff_;
    Teuchos::RCP< std::vector<double> > lambda_;
    int face_index_;
    double face_Mff_;
    double cell_p_;
    double bc_flux_;
    double g_flux_;
    int dir_;
    double patm_;
    Teuchos::RCP<FlowRelations::WRM> wrm_;
  };

  struct Tol_ {
    Tol_(double eps) : eps_(eps) {}
    bool operator()(const double& a, const double& b) const {
      return std::abs(a - b) <= eps_;
    }
    double eps_;
  };

 protected:
  Teuchos::RCP<FluxBCFunctor> CreateFunctor_(int f,
          const Teuchos::Ptr<const CompositeVector>& pres);
  int CalculateLambdaToms_(int f, const Teuchos::Ptr<const CompositeVector>& pres,
                             double& lambda);

 protected:
  Teuchos::RCP<const State> S_next_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<Operators::MatrixMFD> matrix_;
  Teuchos::RCP<FlowRelations::WRMPartition> wrms_;

  std::vector<Operators::Matrix_bc>* bc_markers_;
  std::vector<double>* bc_values_;

};

} // namespace
} // namespace

#endif
