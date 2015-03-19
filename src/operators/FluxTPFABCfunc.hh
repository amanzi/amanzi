/*
  This is the operators component of the Amanzi code. 

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Daniil Svyatskiy (dasvyat@lanl.gov)
*/

#ifndef AMANZI_FLUX_TPFA_BC_FUNC
#define AMANZI_FLUX_TPFA_BC_FUNC

#include <boost/math/tools/roots.hpp>

namespace Amanzi {
namespace Operators {

template <class Nonlin_rcp_ptr>
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
		      const Nonlin_rcp_ptr& nln) :
      trans_f_(trans_f), lambda_(lambda), face_index_(face_index),
      cell_p_(cell_p), bc_flux_(bc_flux), g_flux_(g_flux),
      dir_(dir), patm_(patm) {
    nln_rcp_ptr = nln;
  }

  double operator()(double face_p) {
    lambda_ = face_p;
    double Krel = nln_rcp_ptr->k_relative(patm_ - lambda_);
    double q = flux_();
    return flux_() + g_flux_*Krel - bc_flux_;
  }

 protected:
  double flux_() {
    double s = 0.;
    double Krel = nln_rcp_ptr->k_relative(patm_ - lambda_);

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
  Nonlin_rcp_ptr nln_rcp_ptr;
};

struct Tol_ {
  Tol_(double eps) : eps_(eps) {};
  bool operator()(const double& a, const double& b) const {
    return std::abs(a - b) <= eps_;
  }
  double eps_;
};


template <class WRM>
double OperatorDiffusionFV::DeriveBoundaryFaceValue(
    int f, const CompositeVector& u, const WRM& wrm)
{
  if (u.HasComponent("face")) {
    const Epetra_MultiVector& u_face = *u.ViewComponent("face");
    return u_face[f][0];
  } else {
    const Epetra_MultiVector& trans_face = *transmissibility_->ViewComponent("face", true);
    const Epetra_MultiVector& gravity_face = *gravity_term_->ViewComponent("face", true);

    const std::vector<int>& bc_model = bcs_[0]->bc_model();
    const std::vector<double>& bc_value = bcs_[0]->bc_value();

    if (bc_model[f] == OPERATOR_BC_DIRICHLET) {
      return bc_value[f];
    } else if (bc_model[f] == OPERATOR_BC_NEUMANN) {
      AmanziMesh::Entity_ID_List cells, faces;
      std::vector<int> dirs;
      const Epetra_MultiVector& u_cell = *u.ViewComponent("cell");      
      const Epetra_MultiVector& Krel_face = *k_->ViewComponent("face");

      mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
      int c = cells[0];
      
      mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
      for (int i=0; i<faces.size(); i++) {
      	if (faces[i] == f) {
      	  double a = dirs[i] * trans_face[0][f];
      	  double b = bc_value[f]* mesh_->face_area(f);	  
      	  double face_val = u_cell[0][c] + gravity_face[0][f]/a - b/(a*Krel_face[0][f]);
	  double trans = trans_face[0][f];
	  double gflux = gravity_face[0][f];
	  double bc_flux = bc_value[f]* mesh_->face_area(f);
	  double atm_pressure_ = 101325;

	  FlowFluxTPFA_BCFunc<WRM> func(trans, face_val, i, u_cell[0][c],
					bc_flux, gflux, dirs[i], atm_pressure_, wrm);
          // -- convergence criteria
	  double eps = std::max(1.e-4 * std::abs(bc_flux), 1.e-8);
	  Tol_ tol(eps);
	  boost::uintmax_t max_it = 100;
	  boost::uintmax_t actual_it(max_it);
	  
	  double res = func(face_val);
	  double left = 0.;
	  double right = 0.;
	  double lres = 0.;
	  double rres = 0.;

	  if (res > 0.) {
	    left = face_val;
	    lres = res;
	    right = std::max(face_val, atm_pressure_);
	    rres = func(right);
	    while (rres > 0.) {
	      right += atm_pressure_;
	      rres = func(right);
	    }
	  } else {
	    right = face_val;
	    rres = res;
#if DEBUG_FLAG
	    std::cout << "RIGHT = " << right << ", " << rres << std::endl;
#endif
	    left = std::min(101325., face_val);
	    lres = func(left);
	    while (lres < 0.) {
#if DEBUG_FLAG
	      std::cout << "LEFT = " << left << ", " << lres << std::endl;
#endif
	      left -= 101325.;
	      lres = func(left);
	    }
	  }
#if DEBUG_FLAG
  std::cout << "   bracket (res): " << left << " (" << lres << "), "
            << right << " (" << rres << ")" << std::endl;
#endif

          std::pair<double,double> result =
	    boost::math::tools::toms748_solve(func, left, right, lres, rres, tol, actual_it);
	  if (actual_it >= max_it) {
	    std::cout << " Failed to converged in " << actual_it << " steps." << std::endl;
	    return 3;
	  }
	  
	  face_val = (result.first + result.second) / 2.;

#if DEBUG_FLAG
	  std::cout << "face_val = "<<face_val<<"\n";
#endif
        return face_val;
        }
      }

    } else {
      const Epetra_MultiVector& u_cell = *u.ViewComponent("cell");
      AmanziMesh::Entity_ID_List cells;
      mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
      int c = cells[0];
      return u_cell[0][c];
    }
  }
}

}  // namespace Operators
}  // namespace Amanzi

#endif
