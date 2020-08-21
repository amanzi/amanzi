/* -------------------------------------------------------------------------
  Amanzi & ATS

  License: see $ATS_DIR/COPYRIGHT
  Author: Ethan Coon

  Interface for a Field.  A field contains meta-data about a data structure
  and that data structure, encapsulating things like vis, checkpointing,
  ownership, initialization, etc.
------------------------------------------------------------------------- */

#ifndef AMANZI_STATE_FIELD_HH_
#define AMANZI_STATE_FIELD_HH_

#include <string>
#include <vector>

// TPLs
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"

// State
#include "Checkpoint.hh"
#include "CompositeVector.hh"
#include "StateDefs.hh"
#include "Visualization.hh"

namespace Amanzi {

class Field {

 private: 
  typedef std::map<Key, Teuchos::RCP<Field> > FieldMap;
  typedef std::string Units;

 public:
  // virtual destructor
  virtual ~Field() {/*std::cout<<"delete "<<fieldname_<<"\n";*/};

  // copy constructor and assignment
  virtual Teuchos::RCP<Field> Clone() const = 0;
  virtual Teuchos::RCP<Field> Clone(std::string fieldname) const = 0;
  virtual Teuchos::RCP<Field> Clone(std::string fieldname, std::string owner) const = 0;

  // access
  std::string fieldname() const { return fieldname_; }
  std::string owner() const { return owner_; }
  std::string vis_key() const { return vis_key_; }
  FieldType type() const { return type_; }
  bool io_checkpoint() const { return io_checkpoint_; }
  bool io_vis() const { return io_vis_; }
  bool initialized() const { return initialized_; }

  // mutators
  void set_owner(std::string owner) { owner_ = owner; }
  void set_vis_key(std::string key) { vis_key_ = key; }
  void set_io_checkpoint(bool io_checkpoint=true) { io_checkpoint_ = io_checkpoint; }
  void set_io_vis(bool io_vis=true) { io_vis_ = io_vis; }
  void set_initialized(bool initialized=true) { initialized_ = initialized; }

  void RequireCopy(Key tag);
  void RequireCopy(Key tag, Key new_owner);
  Teuchos::RCP<const Field> GetCopy(Key tag) const;
  Teuchos::RCP<Field> GetCopy(Key tag, Key pk_name);
  void SwitchCopies(Key tag1, Key tag2);
  bool HasCopy(Key tag) const;
  bool DeleteCopy(Key tag) const;
  void SetCopy(Key tag, const Teuchos::RCP<Field>& field);


  // Creation of the actual data (data is created lazily, allowing empty fields).
  virtual void CreateData() = 0;

  // data access and mutators
  virtual Teuchos::RCP<const double> GetScalarData() const {
    assert_type_or_die_(CONSTANT_SCALAR);
    not_implemented_error_();
    return Teuchos::null;
  }
  virtual Teuchos::RCP<double> GetScalarData() {
    assert_type_or_die_(CONSTANT_SCALAR);
    not_implemented_error_();
    return Teuchos::null;
  }
  virtual Teuchos::RCP<const Epetra_Vector> GetConstantVectorData() const {
    assert_type_or_die_(CONSTANT_VECTOR);
    not_implemented_error_();
    return Teuchos::null;
  }
  virtual Teuchos::RCP<Epetra_Vector> GetConstantVectorData() {
    assert_type_or_die_(CONSTANT_VECTOR);
    not_implemented_error_();
    return Teuchos::null;
  }
  virtual Teuchos::RCP<const CompositeVector> GetFieldData() const {
    assert_type_or_die_(COMPOSITE_VECTOR_FIELD);
    not_implemented_error_();
    return Teuchos::null;
  }
  virtual Teuchos::RCP<CompositeVector> GetFieldData() {
    assert_type_or_die_(COMPOSITE_VECTOR_FIELD);
    not_implemented_error_();
    return Teuchos::null;
  }

  // set data by pointer -- does not copy
  virtual void SetData(const Teuchos::RCP<CompositeVector>& data) {
    assert_type_or_die_(COMPOSITE_VECTOR_FIELD);
    not_implemented_error_();
  }
  virtual void SetData(const Teuchos::RCP<Epetra_Vector>& data) {
    assert_type_or_die_(CONSTANT_VECTOR);
    not_implemented_error_();
  }
  virtual void SetData(const Teuchos::RCP<double>& data) {
    assert_type_or_die_(CONSTANT_SCALAR);
    not_implemented_error_();
  }

  // set data by reference -- copies
  virtual void SetData(const CompositeVector& data) {
    assert_type_or_die_(COMPOSITE_VECTOR_FIELD);
    not_implemented_error_();
  }
  virtual void SetData(const Epetra_Vector& data) {
    assert_type_or_die_(CONSTANT_VECTOR);
    not_implemented_error_();
  }
  virtual void SetData(const double& data) {
    assert_type_or_die_(CONSTANT_SCALAR);
    not_implemented_error_();
  }

  // Initialize from a parameter list.
  virtual void Initialize(Teuchos::ParameterList& plist) = 0;

  // Initialize from a checkpoint file
  virtual bool ReadCheckpoint(HDF5_MPI& file) { return true; }

  // write data to the visualization file
  virtual void WriteVis(Visualization& vis) = 0;
  virtual void WriteCheckpoint(Checkpoint& chk) = 0;

  // Compute from a function, if implemented by a subclass
  virtual void Compute(double time) {};

  virtual long int GetLocalElementCount() { return 0L; }

  // Iterate over field_copies.
  typedef FieldMap::const_iterator copy_iterator;
  copy_iterator copy_begin() const { return field_copy_.begin(); }
  copy_iterator copy_end() const { return field_copy_.end(); }
  FieldMap::size_type copy_count() { return field_copy_.size(); }
  void add_copy(Key tag, Teuchos::RCP<Field> field_copy) {field_copy_[tag] = field_copy;}


 protected:
  // constructors protected, should only be called by derived classes
  Field(std::string fieldname, std::string owner="state");
  Field(const Field& other);


  // consistency checking
  void assert_type_or_die_(FieldType type) const;
  void not_implemented_error_() const;

  Teuchos::RCP<const Field> GetCopy_(Key tag) const;
  Teuchos::RCP<Field> GetCopy_(Key tag);

  FieldType type_;
  std::string fieldname_;
  std::string owner_;
  std::string vis_key_;

  bool io_checkpoint_;
  bool io_vis_;
  bool initialized_;

  FieldMap field_copy_;

  Units units_;
  

private:
  // operator= disabled
  Field& operator=(const Field&);
};

} // namespace Amanzi

#endif
