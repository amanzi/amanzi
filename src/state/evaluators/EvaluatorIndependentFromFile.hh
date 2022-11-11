/*
  State

  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon
*/
//! An evaluator with no dependencies specified by discrete data in a file.

/*!

This evaluator is typically used for providing data that are functions of space
and time.  Data is provided, discretely (e.g. with one data point per
cell/face/node), at a series of time slices.  The time slices are interpolated
linearly in time to provide the value.

Within the file, data is expected to meet the following (HDF5) layout:

..code::

   /time : a 1D array of length NTIMES, providing the time in seconds.
   /variable_name.ENTITY.DOF  (group)

      /0 : a 1D array of length NENTITIES, providing the values for each entity
           at time /time[0]
      /1 : ...
      /NTIMES-1 : 1D array at time /time[NTIMES-1]

This evaluator is used by providing the option:

`"evaluator type`" == `"independent variable`"

.. _independent-variable-from-file-evaluator-spec:
.. admonition:: independent-variable-from-file-evaluator-spec

   * `"filename`" ``[string]`` Path to the file.
   * `"variable name`" ``[string]`` Name of the dataset to read from the file.
   * `"domain name`" ``[string]`` **domain** Name of the domain on which the
      field is defined.
   * `"component name`" ``[string]`` **cell** Name of the component in the
     field to populate.
   * `"mesh entity`" ``[string]`` **cell** Name of the entity on which the
     component is defined.
   * `"number of dofs`" ``[int]`` **1** Number of degrees of freedom to read.
   * `"time function`" ``[function-spec]`` **optional** If provided, time is
     first manipulated by this function before interpolation.  This is useful
     for things like cyclic data, which can use a modulo time function to
     repeat the same data.

 */


  TODO: This needs a test and documentation! --etc
*/

#ifndef AMANZI_STATE_INDEPENDENT_FIELD_EVALUATOR_FROMFILE_
#define AMANZI_STATE_INDEPENDENT_FIELD_EVALUATOR_FROMFILE_

#include "EvaluatorIndependent.hh"
#include "Evaluator_Factory.hh"

namespace Amanzi {

class Function;

class EvaluatorIndependentFromFile
    : public EvaluatorIndependent<CompositeVector, CompositeVectorSpace> {

public:
  // ---------------------------------------------------------------------------
  // Constructors
  // ---------------------------------------------------------------------------
  explicit EvaluatorIndependentFromFile(Teuchos::ParameterList &plist);

  EvaluatorIndependentFromFile(const EvaluatorIndependentFromFile &other) =
      default;

  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual Evaluator &operator=(const Evaluator &other) override;

  EvaluatorIndependentFromFile &
  operator=(const EvaluatorIndependentFromFile &other);

  virtual void EnsureCompatibility(State &S) override;

protected:
  // ---------------------------------------------------------------------------
  // Update the value in the state.
  // ---------------------------------------------------------------------------
  virtual void Update_(State &S) override;

  // implementation
  void LoadFile_(int i);
  void Interpolate_(double time, CompositeVector &v);

protected:
  double t_before_, t_after_;
  Teuchos::RCP<CompositeVector> val_before_, val_after_;
  std::string filename_;
  std::vector<double> times_;
  int current_interval_;

  std::string meshname_;
  std::string compname_;
  std::string varname_;
  std::string locname_;
  int ndofs_;

  Teuchos::RCP<Function> time_func_;

private:
  static Utils::RegisteredFactory<Evaluator, EvaluatorIndependentFromFile> fac_;
};

} // namespace Amanzi

#endif
