/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Ethan Coon

   Interface for CompositeVector, an implementation of a slightly improved
   Epetra_MultiVector which spans multiple simplices and knows how to
   communicate itself.
   ------------------------------------------------------------------------- */

#ifndef COMPOSITEVECTOR_HH_
#define COMPOSITEVECTOR_HH_

#include <vector>
#include "Teuchos_RCP.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_CombineMode.h"
#include "Epetra_Import.h"

#include "Mesh.hh"
#include "Vector.hh"
namespace Amanzi {

  typedef enum { CONSTRUCT_WITHOUT_DATA } ConstructMode;

  class CompositeVector : public Vector {

  public:
    // Constructors
    CompositeVector(Teuchos::RCP<AmanziMesh::Mesh>& mesh, std::vector<std::string> names,
                    std::vector<AmanziMesh::Entity_kind> locations,
                    std::vector<int> num_dofs, bool ghosted=true);

    CompositeVector(Teuchos::RCP<AmanziMesh::Mesh>& mesh, std::vector<std::string> names,
                    std::vector<AmanziMesh::Entity_kind> locations,
                    int num_dofs=1, bool ghosted=true);

    CompositeVector(Teuchos::RCP<AmanziMesh::Mesh>& mesh, std::string name,
                    AmanziMesh::Entity_kind location, int num_dofs=1, bool ghosted=true);

    // copy constructor / assignment
    CompositeVector(const CompositeVector& other);
    CompositeVector(const CompositeVector& other, ConstructMode mode);
    CompositeVector& operator=(const CompositeVector& other);

    // access
    std::vector<std::string> names() const { return names_; }
    std::vector< std::vector<std::string> > subfield_names() const { return subfield_names_; }
    std::vector<int> num_dofs() const { return num_dofs_; }
    int num_components() const { return num_components_; }
    Teuchos::RCP<AmanziMesh::Mesh> mesh() const { return mesh_; }

    // mutators
    void set_subfield_names(std::vector< std::vector<std::string> > subfield_names) {
      subfield_names_ = subfield_names; }

    // view data
    // -- Access a view of a single component's data.
    // Ghosted views are simply the vector itself, while non-ghosted views are
    // lazily generated.
    Teuchos::RCP<const Epetra_MultiVector> ViewComponent(std::string name="",
            bool ghosted=false) const;
    Teuchos::RCP<Epetra_MultiVector> ViewComponent(std::string="", bool ghosted=false);

    // -- Access a view of the owned composite data.
    // All meta-data is copied, but pointers to the data are replaced by
    // pointers to owned data.
    Teuchos::RCP<const CompositeVector> ViewOwned() const;
    Teuchos::RCP<CompositeVector> ViewOwned();

    // communicate
    // -- Scatter master values to ghosted values.
    // Modes shown in Epetra_CombineMode.h, but the default is Insert, which
    // overwrites the current ghost value with the (unique) new master value.
    void ScatterMasterToGhosted(Epetra_CombineMode mode=Insert);
    void ScatterMasterToGhosted(std::string name, Epetra_CombineMode mode=Insert);

    // -- Combine ghosted values back to master values.
    // Modes shown in Epetra_CombineMode.h, but the default is InsertAdd,
    // where off-process values are first summed, then replace the current
    // value.
    void GatherGhostedToMaster(Epetra_CombineMode mode=InsertAdd);
    void GatherGhostedToMaster(std::string name, Epetra_CombineMode mode=InsertAdd);

    // assorted vector operations for use by time integrators?
    // These may make life much easier for time integration of the full
    // CompositeVector.  For instance, one could implement the BDF integrators
    // doing all operations by first calling ViewOwned() to get a non-ghosted
    // CompositeVector, and then using these operations.

    // -- Insert value into data.
    int PutScalar(double scalar);
    int PutScalar(std::vector<double> scalar);
    int PutScalar(std::string name, double scalar);
    int PutScalar(std::string name, std::vector<double> scalar);
    int PutScalar(int set_id, double scalar);
    int PutScalar(int set_id, std::vector<double> scalar);
    int PutScalar(std::string name, int set_id, double scalar);
    int PutScalar(std::string name, int set_id, std::vector<double> scalar);

    // -- this <- value*this
    int Scale(double value);
    int Scale(std::string name, double value);

    // -- this <- this + scalarA
    int Shift(double scalarA);
    int Shift(std::string name, double scalarA);

    // -- result <- other \dot this
    int Dot(const CompositeVector& other, double* result) const;

    // -- this <- scalarA*A + scalarThis*this
    CompositeVector& Update(double scalarA, const CompositeVector& A, double scalarThis);

    // -- this <- scalarA*A + scalarB*B + scalarThis*this
    CompositeVector& Update(double scalarA, const CompositeVector& A,
                            double scalarB, const CompositeVector& B, double scalarThis);

    // -- this <- scalarAB * A@B + scalarThis*this  (@ is the elementwise product
    int Multiply(double scalarAB, const CompositeVector& A, const CompositeVector& B,
                 double scalarThis);

    // -- norms
    int NormInf(double* norm) const;
    int Norm1(double* norm) const;
    int Norm2(double* norm) const;

    // Extras
    void Print(ostream &os) const;

  private:
    void VerifyAndCreateData_();
    void CreateOwnedView_(unsigned int index) const;

    Teuchos::RCP<AmanziMesh::Mesh> mesh_;

    unsigned int num_components_;
    mutable std::map< std::string, unsigned int > indexmap_;
    bool ghosted_;

    std::vector< std::string > names_;
    std::vector< std::vector<std::string> > subfield_names_;
    std::vector< Teuchos::RCP<Epetra_MultiVector> > data_;
    mutable std::vector< Teuchos::RCP<Epetra_MultiVector> > owned_data_; // generated lazily
    std::vector< AmanziMesh::Entity_kind > locations_;
    std::vector< int > num_dofs_;
    std::vector< int > cardinalities_;
    std::vector< Teuchos::RCP<Epetra_Import> > importers_;

    mutable Teuchos::RCP<CompositeVector> owned_composite_; // generated lazily
    mutable Teuchos::RCP<Epetra_Vector> owned_vector_; // generated lazily
  };

} // namespace

#endif
