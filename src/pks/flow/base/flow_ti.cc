/*
This is the flow component of the Amanzi code. 
License: BSD
Authors: Neil Carlson (version 1) 
         Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
         Ethan Coon (ATS version) (ecoon@lanl.gov)
*/

#include "flow.hh"

namespace Amanzi {
namespace Flow {

// computes a norm on u-du and returns the result
double Flow::enorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> du) {
  double enorm_val = 0.0;

  std::vector<std::string> names = u->names();
  for (int i=0; i!=names.size(); ++i) {
    Teuchos::RCP<const Epetra_MultiVector> comp = u->data()->ViewComponent(names[i],false);
    Teuchos::RCP<const Epetra_MultiVector> comp_dot = du->data()->ViewComponent(names[i],false);

    for (unsigned int lcv_vec=0; lcv!=comp->NumVector(); ++lcv_vec) {
      for (unsigned int i=0; i!=comp->MyLength(); ++i) {
        double tmp = abs((*(*comp_dot)(lcv_vec))[i]) /
          (atol_ + rtol_*abs((*(*comp)(lcv_vec))[i]));
        enorm_val = std::max<double>(enorm_val, tmp);
      }
    }
  }

#ifdef HAVE_MPI
  double buf = enorm_val;
  MPI_Allreduce(&buf, &enorm_val, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif

  return enorm_val;
};

} // namespace
} // namespace
