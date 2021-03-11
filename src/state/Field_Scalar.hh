/*
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
*/
// Interface for a Field containing a CompositeVector.

/*

  Field_Scalar is metadata for a scalar, but also includes functionality to
initialize, restart, visualize, etc the data.  Parameters available here are
provided for use in InitialConditions_ specs.

*/

/*!

``[initial-conditions-scalar-spec]``

* `"value`" ``[double]`` Set the value.


Example:

.. code-block:: xml

  <ParameterList name="atmospheric_pressure">
    <Parameter name="value" type="double" value="101325."/>
  </ParameterList>


*/




#ifndef AMANZI_STATE_FIELD_SCALAR_HH_
#define AMANZI_STATE_FIELD_SCALAR_HH_

#include <string>
#include "Teuchos_RCP.hpp"

#include "Field.hh"

namespace Amanzi {

class Function;

class Field_Scalar : public Field {

public:
  // constructors
  Field_Scalar(std::string fieldname, std::string owner);
  Field_Scalar(std::string fieldname, std::string owner, Teuchos::RCP<double>& data);

  // copy constructor and assignment
  explicit Field_Scalar(const Field_Scalar& other);
  virtual Teuchos::RCP<Field> Clone() const;
  virtual Teuchos::RCP<Field> Clone(std::string fieldname) const;
  virtual Teuchos::RCP<Field> Clone(std::string fieldname, std::string owner) const;

  // destructor
  ~Field_Scalar() {}

  // Creation of the actual data (data is created lazily, allowing empty fields).
  virtual void CreateData();

  // data access and mutators
  virtual Teuchos::RCP<const double> GetScalarData() const { return data_; }
  virtual Teuchos::RCP<double> GetScalarData();

  virtual void SetData(const Teuchos::RCP<double>& data);
  virtual void SetData(const double& data);

  // initialization
  virtual void Initialize(Teuchos::ParameterList& plist);

  // vis writing
  virtual void WriteVis(Visualization& vis);
  virtual void WriteCheckpoint(Checkpoint& ckp);

  // Compute from a function
  virtual void Compute(double time);

 protected:

  Teuchos::RCP<double> data_;
  Teuchos::RCP<Function> func_;

private:
  // operator= disabled
  Field_Scalar& operator=(const Field_Scalar&);

}; // class Field

} // namespace Amanzi

#endif
