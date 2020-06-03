/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/
//!  Base class for preconditioners.

/*!

Provides approximate inverses of matrices.

Note that preconditioners here differ from linear operators not in the
approximate nature of their inverse, but in their interface.  Preconditioners
work with raw vectors and matrices, and may need assembled matrices.  A `Linear
Solver`_ works with the action of matrices only, and never need assembled
matrices.  As such they are templated with an arbitrary Matrix and Vector type,
whereas Preconditioners are not.


.. _preconditioner-typed-spec:
.. admonition:: preconditioner-typed-spec

    * `"preconditioner type`" ``[string]`` **identity** Iterative method to be used.
    * `"_preconditioner_type_ parameters`" ``[_preconditioner_type_-spec]``
      Parameters associated with the requested preconditioner.

Example:

.. code-block:: xml

     <ParameterList name="my preconditioner">
       <Parameter name="type" type="string" value="trilinos ml"/>
        <ParameterList name="trilinos ml parameters"> ?????? check me!
            ... 
        </ParameterList>
     </ParameterList>
      
 */

#ifndef AMANZI_PRECONDITIONER_HH_
#define AMANZI_PRECONDITIONER_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_RowMatrix.h"

#include "exceptions.hh"

namespace Amanzi {
namespace AmanziPreconditioners {

class Preconditioner {
 public:
  virtual ~Preconditioner() = default;

  // Initializes the solver with provided parameters.
  // This need not be called by preconditioners created using the factory.
  virtual void Init(const std::string& name,
                    const Teuchos::ParameterList& list) = 0;

  // Rebuild the preconditioner using the given matrix A.
  virtual void Update(const Teuchos::RCP<Epetra_RowMatrix>& A) = 0;

  // Destroy the preconditioner and auxiliary data structures.
  virtual void Destroy() = 0;

  // Apply the preconditioner.
  virtual int ApplyInverse(const Epetra_MultiVector& v,
                           Epetra_MultiVector& hv) = 0;

  virtual int returned_code() = 0;
};

}  // namespace AmanziPreconditioners
}  // namespace Amanzi


#endif


