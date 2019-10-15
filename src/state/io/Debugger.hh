/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//!

#ifndef AMANZI_DEBUGGER_HH_
#define AMANZI_DEBUGGER_HH_

#include <string>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Ptr.hpp"

#include "Mesh.hh"
#include "VerboseObject.hh"

namespace Amanzi {

template <typename Scalar>
class CompositeVector_;
using CompositeVector = CompositeVector_<double>;

class Debugger {
 public:
  // Constructor
  Debugger(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh, std::string name,
           Teuchos::ParameterList& plist,
           Teuchos::EVerbosityLevel verb_level = Teuchos::VERB_HIGH);

  // Write cell + face info
  void WriteCellInfo(bool include_faces = false);

  // Write a vector individually.
  void WriteVector(const std::string& name,
                   const Teuchos::Ptr<const CompositeVector>& vec,
                   bool include_faces = false);

  // Write boundary condition data.
  // void WriteBoundaryConditions(const std::vector<int> &flag,
  //                              const std::vector<double> &data);

  // Write list of vectors.
  void
  WriteVectors(const std::vector<std::string>& names,
               const std::vector<Teuchos::Ptr<const CompositeVector>>& vecs,
               bool include_faces = false);

  // call MPI_Comm_Barrier to sync between writing steps
  void Barrier();

  // write a line of ----
  void WriteDivider();

  void SetPrecision(int prec) { precision_ = prec; }
  void SetWidth(int width) { width_ = width; }

  // reverse order -- instead of passing in vector, do writing externally
  Teuchos::RCP<VerboseObject> GetVerboseObject(AmanziMesh::Entity_ID, int rank);

 protected:
  std::string Format_(double dat);
  std::string FormatHeader_(std::string name, int c);

 protected:
  Teuchos::EVerbosityLevel verb_level_;
  Teuchos::RCP<VerboseObject> vo_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  std::vector<AmanziMesh::Entity_ID> dc_;
  std::vector<AmanziMesh::Entity_ID> dc_gid_;
  std::vector<Teuchos::RCP<VerboseObject>> dcvo_;

  int width_;
  int header_width_;
  int precision_;
  int cellnum_width_;
  int decimal_width_;
};

} // namespace Amanzi

#endif
