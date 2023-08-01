/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

//! A base class for assembled preconditioners.

/*!

This sublist contains entries for various
preconditioners required by a simulation. At the moment, we support Trilinos multilevel
preconditioner, Hypre BoomerAMG preconditioner, ILU preconditioner, Hypre's ILU
preconditioner, and identity preconditioner. 

* `"preconditioning method`" [string] defines preconditioner algorithm.

* `"xxx parameters`" [list] provides parameters for the preconditioner specified 
  by parameter `"preconditioning method`".
 
.. code-block:: xml

  <ParameterList>  <!-- parent list -->
  <ParameterList name="preconditioners">
    <ParameterList name="_TRILINOS ML">
      <Parameter name="preconditioning method" type="string" value="ml"/>
      <ParameterList name="ml parameters">
        ... 
      </ParameterList>
    </ParameterList>

    <ParameterList name="_HYPRE AMG">
      <Parameter name="preconditioning method" type="string" value="boomer amg"/>
      <ParameterList name="boomer amg parameters">
        ...
      </ParameterList>
    </ParameterList>

    <ParameterList name="_BLOCK ILU">
      <Parameter name="preconditioning method" type="string" value="block ilu"/>
      <ParameterList name="block ilu parameters">
        ...
      </ParameterList>
    </ParameterList>

    <ParameterList name="_DIAGONAL">
      <Parameter name="preconditioning method" type="string" value="diagonal"/>
    </ParameterList>
  </ParameterList>

*/

#pragma once

#include "Inverse.hh"

class Epetra_CrsMatrix;
class Epetra_MultiVector;
class Epetra_Map;

namespace Amanzi {
namespace AmanziSolvers {

using Preconditioner = Inverse<Epetra_CrsMatrix, Epetra_CrsMatrix, Epetra_Vector, Epetra_Map>;

} // namespace AmanziSolvers
} // namespace Amanzi
