/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! An evaluator for a WhetStone::Tensor field, specified by discrete data in
//! an HDF5 file.
/*
  State

*/

/*!

This evaluator is typically used for providing tensor-valued data (e.g. a
permeability tensor) that are functions of space and time.  Data is
provided, discretely (e.g. with one data point per cell), at a series of
time slices.  The time slices are interpolated linearly in time to provide
the value.

The tensor's shape is controlled by `"tensor type`", which determines how
many degrees of freedom (DoFs) are read per entity, and how those DoFs are
mapped into the tensor's entries:

* `"scalar`" -- isotropic, a single value on the diagonal (1 DoF).
* `"horizontal and vertical`" -- diagonal, with identical x-y entries and a
  distinct z entry (2 DoF, 3D only).
* `"diagonal`" -- diagonal, one value per dimension (2 DoF in 2D, 3 DoF in 3D).
* `"full symmetric`" -- symmetric tensor: diagonal entries followed by the
  off-diagonal entries (3 DoF in 2D, 6 DoF in 3D).
* `"full`" -- general, nonsymmetric tensor, all entries in row-major order
  (4 DoF in 2D, 9 DoF in 3D).

Within the file, data is expected to meet the following (HDF5) layout::

   /time : a 1D array of length NTIMES, providing the time in seconds.
   /variable_name.ENTITY.DOF  (group)

      /0 : a 1D array of length NENTITIES, providing the values for each entity
           at time /time[0]
      /1 : ...
      /NTIMES-1 : 1D array at time /time[NTIMES-1]

where DOF ranges over the number of degrees of freedom implied by `"tensor
type`" above, in the order given above (e.g. for `"full symmetric`" in 3D:
xx, yy, zz, xy, xz, yz).

This evaluator is used by providing the option:

`"evaluator type`" == `"independent variable tensor from file`"

.. _independent-variable-tensor-from-file-evaluator-spec:
.. admonition:: independent-variable-tensor-from-file-evaluator-spec

   * `"tensor type`" ``[string]`` One of `"scalar`", `"horizontal and
     vertical`", `"diagonal`", `"full symmetric`", or `"full`" (see above).
   * `"filename`" ``[string]`` Path to the file.
   * `"variable name`" ``[string]`` Name of the dataset to read from the file.
   * `"domain name`" ``[string]`` **domain** Name of the domain on which the
     field is defined.
   * `"component name`" ``[string]`` **cell** Name of the component in the
     field to populate.
   * `"mesh entity`" ``[string]`` **cell** Name of the entity on which the
     component is defined.
   * `"constant in time`" ``[bool]`` **false** Is the value constant throughout all time?
   * `"page size`" ``[int]`` **2** Number of consecutive time slices to keep
     resident in memory at once.  Larger values amortize file reads across
     more timesteps at the cost of more memory; must be at least 2.
   * `"file format`" ``[string]`` **optional** `"HDF5`" or `"NetCDF`" --
     inferred from the filename's extension (".h5" or ".nc") if not provided.
   * `"checkpoint file`" ``[bool]`` **false** If this is true, then it is
     expected that `"filename`" is a checkpoint-file-like object, where
     /variable_name.ENTITY.DOF is itself a single 1D array, and not a group
     of per-time-index arrays.  Note this forces `"constant in time`" to be
     true.
   * `"time function`" ``[function-spec]`` **optional** If provided, time is
     first manipulated by this function before interpolation.  This is useful
     for things like cyclic data, which can use a modulo time function to
     repeat the same data.

.. code-block:: xml

  <ParameterList name="field_evaluators">  <!-- parent list -->
  <ParameterList name="permeability">
    <Parameter name="field evaluator type" type="string" value="independent variable tensor from file"/>
    <Parameter name="tensor type" type="string" value="diagonal"/>
    <Parameter name="filename" type="string" value="_DATA_FILE.h5"/>
    <Parameter name="domain name" type="string" value="domain"/>
    <Parameter name="variable name" type="string" value="permeability"/>
    <Parameter name="component name" type="string" value="cell"/>
    <Parameter name="mesh entity" type="string" value="cell"/>
  </ParameterList>
  </ParameterList>

The field *permeability* is defined as a cell-based, diagonal tensor
(dimension DoFs per cell), interpolated in time from the time slices stored
in the file.

*/

#pragma once

#include "EvaluatorIndependent.hh"
#include "Evaluator_Factory.hh"
#include "TensorVector.hh"
#include "EvaluatorFromFile_Helpers.hh"
#include "IOFileFormat.hh"

namespace Amanzi {

class Function;

class EvaluatorIndependentTensorFromFile
  : public EvaluatorIndependent<TensorVector, TensorVector_Factory> {
 public:
  // ---------------------------------------------------------------------------
  // Constructors
  // ---------------------------------------------------------------------------
  explicit EvaluatorIndependentTensorFromFile(const Teuchos::RCP<Teuchos::ParameterList>& plist);

  EvaluatorIndependentTensorFromFile(const EvaluatorIndependentTensorFromFile& other) = default;

  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual Evaluator& operator=(const Evaluator& other) override;

  EvaluatorIndependentTensorFromFile& operator=(const EvaluatorIndependentTensorFromFile& other);

  virtual void EnsureCompatibility(State& S) override;

  static const std::string eval_type;
  virtual std::string getType() const override { return eval_type; }

 protected:
  // ---------------------------------------------------------------------------
  // Update the value in the state.
  // ---------------------------------------------------------------------------
  virtual void Update_(State& S) override;

 protected:
  std::string filename_;

  std::string meshname_;
  std::string compname_;
  std::string varname_;
  std::string locname_;

  std::string tensor_type_;
  int dimension_;
  int rank_;
  int num_funcs_;
  int page_size_;
  Utils::IOFileFormat_kind file_format_;
  double rescaling_;
  bool checkpoint_file_;

  Teuchos::RCP<Function> time_func_;
  Teuchos::RCP<EvaluatorFromFile_Helpers::TimeSeriesInterpolatorCV> interpolator_;

 private:
  static Utils::RegisteredFactory<Evaluator, EvaluatorIndependentTensorFromFile> reg_;
};

} // namespace Amanzi
