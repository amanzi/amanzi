#ifndef AMANZI_FLUX_TPFA_BC_FUNC
#define AMANZI_FLUX_TPFA_BC_FUNC

#include <../pks/flow/WaterRetentionModel.hh>

namespace Amanzi {
namespace Operators {

class FlowFluxTPFA_BCFunc {
   public:
    FlowFluxTPFA_BCFunc(const double trans_f,
                  const double lambda,
                  int face_index,
                  double cell_p,
                  double bc_flux,
                  double g_flux,
                  int dir,
                  double patm,
		  const Teuchos::RCP<Flow::WaterRetentionModel>& wrm) :
        trans_f_(trans_f), lambda_(lambda), face_index_(face_index),
        cell_p_(cell_p), bc_flux_(bc_flux), g_flux_(g_flux),
        wrm_(wrm), dir_(dir), patm_(patm) {
    }

    double operator()(double face_p) {
      lambda_ = face_p;
      double Krel = wrm_->k_relative(patm_ - lambda_);
      //std::cout << "Fluxes: " << std::endl;
      double q = flux_();
      //std::cout << "  K grad p = " << q << ", K grad gz = " << g_flux_*Krel << ", bc = " << bc_flux_ << std::endl;
      return flux_() + g_flux_*Krel - bc_flux_;
    }

   protected:

    double flux_() {
      double s = 0.;
      double Krel = wrm_->k_relative(patm_ - lambda_);

      //      std::cout << "  Krel = " << Krel << std::endl;
      //      std::cout << "  lambda_bc = " << (*lambda_)[face_index_] << std::endl;
      s = dir_ * trans_f_ * Krel * (cell_p_ - lambda_);
      return s;
    }

   protected:
    double trans_f_;
    double lambda_;
    int face_index_;
    double face_Mff_;
    double cell_p_;
    double bc_flux_;
    double g_flux_;
    int dir_;
    double patm_;
  Teuchos::RCP<Flow::WaterRetentionModel> wrm_;
  };

  struct Tol_ {
    Tol_(double eps) : eps_(eps) {}
    bool operator()(const double& a, const double& b) const {
      return std::abs(a - b) <= eps_;
    }
    double eps_;
  };


}
}

#endif
