/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Ethan Coon

   Interface for the State.  State is a simple data-manager, allowing PKs to
   require, read, and write various fields, including:
    -- Acts as a factory for fields through the various require methods.
    -- Provides some data protection by providing both const and non-const
       data pointers to PKs.
    -- Provides some initialization capability -- this is where all
       independent variables can be initialized (as independent variables
       are owned by state, not by any PK).
   ------------------------------------------------------------------------- */

#ifndef STATE_STATE_HH_
#define STATE_STATE_HH_

#include <string>
#include <vector>
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_Map.h"
#include "Epetra_Export.h"
#include "Mesh.hh"
#include "Vis.hpp"
#include "CompositeVector.hh"

#include "Field.hh"
#include "Field_Scalar.hh"
#include "Field_ConstantVector.hh"
#include "Field_CV.hh"

namespace Amanzi {

class State {

public:
  State(Teuchos::RCP<Amanzi::AmanziMesh::Mesh>& mesh);

  State(Teuchos::ParameterList& state_plist,
        Teuchos::RCP<Amanzi::AmanziMesh::Mesh>& mesh);

  // Copy constructor, copies memory not pointers.
  explicit State(const State& other);

  // Assignment operator, copies memory not pointers.  Note this
  // implementation requires the State being copied has the same structure (in
  // terms of fields, order of fields, etc) as *this.  This really means that
  // it should be a previously-copy-constructed version of the State.  One and
  // only one State should be instantiated and populated -- all other States
  // should be copy-constructed from that initial State.
  State& operator=(const State& other);

  // initialize values over blocks
  void Initialize();
  bool CheckAllInitialized();

  // Add things to the state.
  //
  // Note that multiple PKs may require a field, but only one may own it.
  void Require(std::string fieldname, FieldType type, std::string owner="state");

  void RequireScalar(std::string fieldname, std::string owner,
                     Teuchos::RCP<double>& data_ptr);
  void RequireScalar(std::string fieldname, std::string owner, double data);

  void RequireConstantVector(std::string fieldname, std::string owner,
          Teuchos::RCP<Epetra_Vector>& data_ptr);
  void RequireConstantVector(std::string fieldname, std::string owner,
          const Epetra_Vector& data);
  void RequireConstantVector(std::string fieldname, std::string owner, int dimension);

  void RequireField(std::string fieldname, std::string owner,
                     Teuchos::RCP<CompositeVector>& data_ptr);
  void RequireField(std::string fieldname, std::string owner,
                     const CompositeVector& data);
  void RequireField(std::string fieldname, std::string owner, AmanziMesh::Entity_kind location,
                    int num_dofs=1, bool ghosted=true);
  void RequireField(std::string fieldname, std::string owner,
                     std::vector<std::string> names,
                     std::vector<AmanziMesh::Entity_kind> locations,
                     int num_dofs=1, bool ghosted=true);
  void RequireField(std::string fieldname, std::string owner,
                     std::vector<std::string> names,
                     std::vector<AmanziMesh::Entity_kind> locations,
                     std::vector<int> num_dofs, bool ghosted=true);

  // -- access methods -- Const methods should be used by PKs who don't own
  // the field, i.e.  flow accessing a temperature field if an energy PK owns
  // the temperature field.  This ensures a PK cannot mistakenly alter data it
  // doesn't own.  Non-const methods get used by the owning PK.
  Teuchos::RCP<const double> GetScalarData(std::string fieldname) const;
  Teuchos::RCP<double> GetScalarData(std::string fieldname, std::string pk_name);
  Teuchos::RCP<const Epetra_Vector> GetConstantVectorData(std::string fieldname) const;
  Teuchos::RCP<Epetra_Vector> GetConstantVectorData(std::string fieldname,
          std::string pk_name);
  Teuchos::RCP<const CompositeVector> GetFieldData(std::string fieldname) const;
  Teuchos::RCP<CompositeVector> GetFieldData(std::string fieldname, std::string pk_name);

  // Access to the full field record, not just the data.
  Teuchos::RCP<Field> GetRecord(std::string fieldname, std::string pk_name);
  Teuchos::RCP<const Field> GetRecord(std::string fieldname) const;

  Teuchos::RCP<AmanziMesh::Mesh> mesh() { return mesh_; }
  Teuchos::RCP<const AmanziMesh::Mesh> mesh() const { return mesh_; }

  double time () const { return time_; }
  int cycle () const { return cycle_; }

  // modify methods
  // -- modify by pointer, no copy
  void SetData(std::string fieldname, std::string pk_name,
                Teuchos::RCP<double>& data);
  void SetData(std::string fieldname, std::string pk_name,
                Teuchos::RCP<Epetra_Vector>& data);
  void SetData(std::string fieldname, std::string pk_name,
                Teuchos::RCP<CompositeVector>& data);

  // -- modify by reference, with copy
  void SetData(std::string fieldname, std::string pk_name, const double& data);
  void SetData(std::string fieldname, std::string pk_name, const Epetra_Vector& data);
  void SetData(std::string fieldname, std::string pk_name, const CompositeVector& data);

  void set_time ( double new_time ) { time_ = new_time; }
  void advance_time(double dT) { time_ += dT; }
  void set_cycle ( int cycle ) { cycle_ = cycle; }
  void advance_cycle ( int dcycle=1 ) { cycle_ += dcycle; }

  // vis and restart functions
// void WriteVis(Amanzi::Vis& vis);

private:

  void InitializeFromParameterList_();
  Teuchos::RCP<Field> GetRecord_(std::string fieldname);
  Teuchos::RCP<const Field> GetRecord_(std::string fieldname) const;
  Teuchos::RCP<Field> CheckMayCreateOrOwn_or_die_(std::string fieldname, FieldType type);
  void PushBackNewField_(std::string fieldname, FieldType type, std::string owner);

  // field container and fieldname map from name -> container location
  std::vector< Teuchos::RCP<Field> > fields_;
  std::map<std::string, std::vector< Teuchos::RCP<Field> >::size_type> field_name_map_;

  double time_;
  int cycle_;

  // mesh
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh_;

  // parameter list
  Teuchos::ParameterList state_plist_;
};

inline
Teuchos::RCP<Field> State::GetRecord_(std::string fieldname) {
  return fields_[field_name_map_[fieldname]];
};

inline
Teuchos::RCP<const Field> State::GetRecord_(std::string fieldname) const {
  return fields_[field_name_map_.find(fieldname)->second];
};

} // namespace amanzi
#endif
