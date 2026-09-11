/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon, Bo Gao
*/
/*!

This evaluator is typically used for providing tensor type data that are functions of space 
and time. Data is provided, discretely (e.g. with one data point per cell/face/node), at a 
series of time slices. The time slices are interpolated linearly in time to provide the value.

Within the file, data is expected to meet the following (HDF5) layout::

   /time : a 1D array of length NTIMES, providing the time in seconds.
   /variable_name.ENTITY.DOF  (group)

      /0 : a 1D array of length NENTITIES, providing the values for each entity
           at time /time[0]
      /1 : ...
      /NTIMES-1 : 1D array at time /time[NTIMES-1]


`"evaluator type`" == `"independent variable tensor from file`"

.. _evaluator-independent-variable-tensor-from-file-spec:
.. admonition:: evaluator-independent-variable-tensor-from-file-spec

   * `"tensor type`" ``[string]`` One of:
     - `"scalar`" A single scalar, isotropic (1 DoF)
     - `"horizontal and vertical`" Diagonal with identical x-y entries and
       different z entries. (2 DoF)
     - `"diagonal`" Diagonal (3 DoF for 3D, 2 for 2D)
     - `"full symmetric`" Full (symmetric) tensor, (6 DoF for 3D, 3 for 2D)
     - `"full`" Full (nonsymmetric) tensor, (9 DoF for 3D, 6 for 2D)
   * `"filename`" ``[string]`` Path to the file.
   * `"variable name`" ``[string]`` Name of the dataset to read from the file.
   * `"domain name`" ``[string]`` **domain** Name of the domain on which the
     field is defined.
   * `"component name`" ``[string]`` **cell** Name of the component in the
     field to populate.
   * `"mesh entity`" ``[string]`` **cell** Name of the entity on which the
     component is defined.
   * `"constant in time`" ``[bool]`` **true** Is the value constant throughout all time?
   * `"checkpoint file`" ``[bool]`` **false** If this is true, then it is
     expected that `"filename`" is a checkpoint-file-like object, where
     /variable_name.ENTITY.DOF is itself a vector, and not a group.  Note this
     forces `"constant in time`" to be true.
   * `"time function`" ``[function-spec]`` **optional** If provided, time is
     first manipulated by this function before interpolation.  This is useful
     for things like cyclic data, which can use a modulo time function to
     repeat the same data.

.. code-block:: xml

   <ParameterList name="permeability">
     <Parameter name="evaluator type" type="string" value="independent variable tensor from file"/>
     <Parameter name="tensor type" type="string" value="full"/>
     <Parameter name="constant in time" type="bool" value="true"/>
     <Parameter name="filename" type="string" value="_DATA_FILE.h5"/>
     <Parameter name="domain name" type="string" value="domain"/>
     <Parameter name="variable name" type="string" value="permeability"/>
     <Parameter name="component name" type="string" value="cell"/>
     <Parameter name="mesh entity" type="string" value="cell"/>
   </ParameterList>

The field *permeability* is defined as a cell-based constant-in-time variable with 9 DoF.
The file *_DATA_FILE.h5* should include a time dataset and one data group for each tensor 
DoF (9 data groups for a full 3D tensor). Each DoF group contains a single 1D dataset 
storing the values for all mesh entities (e.g., cells).

*/

//


#pragma once

#include "EvaluatorIndependent.hh"
#include "Evaluator_Factory.hh"
#include "TensorVector.hh"
#include "EvaluatorFromFile_Helpers.hh"

namespace Amanzi {

class Function;

class EvaluatorIndependentTensorFromFile
  : public EvaluatorIndependent<TensorVector, TensorVector_Factory> {
 public:
  // ---------------------------------------------------------------------------
  // Constructors
  // ---------------------------------------------------------------------------
  explicit EvaluatorIndependentTensorFromFile(Teuchos::ParameterList& plist);

  EvaluatorIndependentTensorFromFile(const EvaluatorIndependentTensorFromFile& other) = default;

  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual Evaluator& operator=(const Evaluator& other) override;

  EvaluatorIndependentTensorFromFile& operator=(const EvaluatorIndependentTensorFromFile& other);

  virtual void EnsureCompatibility(State& S) override;

 protected:
  // ---------------------------------------------------------------------------
  // Update the value in the state.
  // ---------------------------------------------------------------------------
  virtual void Update_(State& S) override;

  template<typename EvaluatorType>
  friend void EvaluatorFromFile_Helpers::LoadFile(int, EvaluatorType&);

  template<typename EvaluatorType>
  friend void EvaluatorFromFile_Helpers::UpdateTimeInterpolation(double, EvaluatorType&, CompositeVector&);

 protected:
  double t_before_, t_after_;
  Teuchos::RCP<CompositeVector> val_before_, val_after_;
  std::string filename_;
  std::vector<double> times_;
  int current_interval_;

  std::string meshname_;
  std::string compname_;
  std::string varname_;
  AmanziMesh::Entity_kind locname_;
  int num_funcs_;
  int ndofs_;
  int dimension_;
  std::string tensor_type_;
  int rank_;
  double rescaling_;

  bool checkpoint_file_;
  Teuchos::RCP<Function> time_func_;

 private:
  static Utils::RegisteredFactory<Evaluator, EvaluatorIndependentTensorFromFile> fac_;
};

} // namespace Amanzi
