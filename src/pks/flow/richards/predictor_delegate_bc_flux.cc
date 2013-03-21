/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Delegate for modifying the predictor in the case of infiltration into dry soil.

  License: BSD
  Authors: Ethan Coon (ATS version) (ecoon@lanl.gov)
*/

#include <boost/math/tools/roots.hpp>
#include "predictor_delegate_bc_flux.hh"

namespace Amanzi {
namespace Flow {

#define DEBUG_FLAG 1

bool PredictorDelegateBCFlux::modify_predictor(double h, Teuchos::RCP<TreeVector> u) {
  Teuchos::RCP<CompositeVector> pres = u->data();

  int nfaces = bc_values_->size();
  for (int f=0; f!=nfaces; ++f) {
    if ((*bc_markers_)[f] == Operators::MATRIX_BC_FLUX) {
      double lambda = (*pres)("face",f);
      int ierr = CalculateLambdaToms_(f, pres.ptr(), lambda);
      ASSERT(!ierr);
      if (!ierr) (*pres)("face",f) = lambda;
    }
  }
  return true;
}


Teuchos::RCP<PredictorDelegateBCFlux::FluxBCFunctor>
PredictorDelegateBCFlux::CreateFunctor_(int f,
        const Teuchos::Ptr<const CompositeVector>& pres) {
  // inner cell
  AmanziMesh::Entity_ID_List cells;
  mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
  ASSERT(cells.size() == 1);
  int c = cells[0];

  // that cell's faces
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;
  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);

  // index within that cell's faces
  int n = std::find(faces.begin(), faces.end(), f) - faces.begin();
  ASSERT(n != faces.size());

  Teuchos::RCP< std::vector<double> > Aff = Teuchos::rcp(new std::vector<double>());
  Aff->resize(faces.size());
  Teuchos::RCP< std::vector<double> > lambda = Teuchos::rcp(new std::vector<double>());
  lambda->resize(faces.size());

  // fill the arrays
  for (int i=0; i!=faces.size(); ++i) {
    (*Aff)[i] = matrix_->Aff_cells()[c](n,i);
    (*lambda)[i] = (*pres)("face",faces[i]);
  }

  // create and return
  return Teuchos::rcp(new FluxBCFunctor(Aff, lambda, n, (*pres)("cell",c),
          (*bc_values_)[f], dirs[n], 101325.0, wrms_->second[(*wrms_->first)[c]]));
}


int PredictorDelegateBCFlux::CalculateLambdaNewton_(int f,
        const Teuchos::Ptr<const CompositeVector>& pres, double& lambda) {
  Teuchos::RCP<FluxBCFunctor> func = CreateFunctor_(f,pres);

  // -- scaling for the norms
  double tol = std::max(1.e-4 * std::abs((*bc_values_)[f]), 1.e-8);
  double max_steps = 100;
  double stepnum = 0;

  // get the initial residual
  double res = (*func)(lambda);
  double norm = std::abs(res);
#if DEBUG_FLAG
  std::cout << " Flux correcting: " << std::endl;
  std::cout << "   Iter: " << stepnum;
  std::cout << " lambda (res) = " << lambda << " ("
            << res << ")" << std::endl;
#endif

  // check convergence
  bool converged = norm < tol;

  // iterate
  while (!converged) {
    stepnum++;

    double deriv = func->derivative(lambda);
    if (std::abs(deriv) < 1.e-30) {
      std::cout << "  Zero-derivative in calculating flux" << std::endl;
      return 1;
    }


    double correction = res / deriv;
    double new_lambda = lambda - correction;
    std::cout << "    lambda=" << lambda << ", correction=" << correction << ", new_lambda=" << new_lambda << std::endl;

    res = (*func)(new_lambda);
#if DEBUG_FLAG
      std::cout << "  Iter: " << stepnum;
      std::cout << " corrected lambda (res) = " << new_lambda << " ("
                << res << ")" << std::endl;
#endif

    double new_norm = std::abs(res);

    // backtrack
    double damp = 1.;
    while (new_norm > norm) {
      damp *= 0.5;
      new_lambda = lambda - damp * correction;
      res = (*func)(new_lambda);

#if DEBUG_FLAG
      std::cout << "    Damping: " << stepnum;
      std::cout << " corrected lambda (res) = " << new_lambda << " ("
                << res << ")" << std::endl;
#endif

      new_norm = std::abs(res);
    }
    norm = new_norm;
    lambda = new_lambda;

    converged = res < tol;
    if (stepnum > max_steps && !converged) {
      std::cout << "  Nonconverged after " << max_steps << " steps with norm=" << norm << std::endl;
      return 2;
    }
  }

  return 0;
}


int PredictorDelegateBCFlux::CalculateLambdaToms_(int f,
        const Teuchos::Ptr<const CompositeVector>& pres, double& lambda) {
  Teuchos::RCP<FluxBCFunctor> func = CreateFunctor_(f,pres);

  // -- convergence criteria
  double eps = std::max(1.e-4 * std::abs((*bc_values_)[f]), 1.e-8);
  Tol_ tol(eps);
  int max_it = 100;
  uintmax_t actual_it(max_it);

  double res = (*func)(lambda);
  double left = 0.;
  double right = 0.;
  double lres = 0.;
  double rres = 0.;

  if (res > 0.) {
    left = lambda;
    lres = res;
    right = std::max(lambda,101325.);
    rres = (*func)(right);
    while (rres > 0.) {
      right += 101325.;
      rres = (*func)(right);
    }

  } else {
    right = lambda;
    rres = res;

    left = std::min(101325., lambda);
    lres = (*func)(left);
    while (lres < 0.) {
      left -= 101325.;
      lres = (*func)(left);
    }
  }

#if DEBUG_FLAG
  std::cout << " Flux correcting: " << std::endl;
  std::cout << "   bracket (res): " << left << " (" << lres << "), "
            << right << " (" << rres << ")" << std::endl;
#endif

  std::pair<double,double> result =
      boost::math::tools::toms748_solve(*func, left, right, lres, rres, tol, actual_it);
  if (actual_it >= max_it) {
    std::cout << " Failed to converged in " << actual_it << " steps." << std::endl;
    return 3;
  }

  lambda = (result.first + result.second) / 2.;
#if DEBUG_FLAG
  std::cout << "  Converged to " << lambda << " in " << actual_it << " steps." << std::endl;
#endif
  return 0;
}



} // namespace
} // namespace
