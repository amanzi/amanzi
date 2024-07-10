/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! <MISSING_ONELINE_DOCSTRING>
#ifndef AMANZI_OUTPUT_SILO_HH_
#define AMANZI_OUTPUT_SILO_HH_

#include <string>
#include <vector>

#include "Teuchos_ParameterList.hpp"

extern "C"
{
#include "silo.h"
};

#include "errors.hh"
#include "dbc.hh"
#include "Mesh.hh"

#include "Output.hh"

namespace Amanzi {

class OutputSilo : public Output {
 public:
  OutputSilo(Teuchos::ParameterList& plist, const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);

  // destructor must release file resource on non-finalized
  ~OutputSilo();

  // open and close files
  virtual void createTimestep(double time, int cycle) override;
  virtual void finalizeTimestep() override;
  virtual std::string getFilename(int cycle) const override { return formatter_(cycle) + ".silo"; }

  // how nice it would be to template virtual methods...
  virtual void write(const Teuchos::ParameterList& attrs, const int& val) const override {}
  virtual void write(const Teuchos::ParameterList& attrs, const double& val) const override {}
  virtual void write(const Teuchos::ParameterList& attrs, const std::string& val) const override {}

  virtual void write(const Teuchos::ParameterList& attrs, const Vector_type& vec) const override
  {
    writeVector_(attrs, vec);
  }
  virtual void write(const Teuchos::ParameterList& attrs, const IntVector_type& vec) const override
  {
    writeVector_(attrs, vec);
  }

 protected:
  void init_(Teuchos::ParameterList& plist);
  void closeFile_() const;
  void writeMesh_();
  std::string fixName_(const std::string& instring) const;
  int getSiloLocation_(AmanziMesh::Entity_kind) const;

  template <typename vector_type>
  void writeVector_(const Teuchos::ParameterList& attrs, const vector_type& vec) const;

 protected:
  int count_;
  FilenameFormatter formatter_;

  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  mutable DBfile* fid_;
};


namespace Impl {

// a compile-time map from type to Silo datatype enum
template <typename T>
struct Silo_DatatypeMap;
template <>
struct Silo_DatatypeMap<int> {
  static const int type = DB_INT;
};
template <>
struct Silo_DatatypeMap<double> {
  static const int type = DB_DOUBLE;
};


} // namespace Impl


template <typename vector_type>
void
OutputSilo::writeVector_(const Teuchos::ParameterList& attrs, const vector_type& vec) const
{
  int comm_size = mesh_->getComm()->getSize();
  int comm_rank = mesh_->getComm()->getRank();
  std::string varname = fixName_(attrs.name());
  std::string fname = getFilename(count_);
  AmanziMesh::Entity_kind entity_kind = attrs.get<AmanziMesh::Entity_kind>("location");

  for (int rank = 0; rank != comm_size; ++rank) {
    if (rank == comm_rank) {
      fid_ = DBOpen(fname.c_str(), DB_HDF5, DB_APPEND);

      // directory
      std::string dirname = "/domain_";
      dirname = dirname + std::to_string(rank);
      int ierr = DBSetDir(fid_, dirname.c_str());
      AMANZI_ASSERT(!ierr);

      auto vecv = vec.get1dView();
      ierr = DBPutUcdvar1(fid_,
                          varname.c_str(),
                          "mesh",
                          (void*)vecv.get(),
                          vec.getLocalLength(),
                          NULL,
                          0,
                          Impl::Silo_DatatypeMap<typename vector_type::scalar_type>::type,
                          getSiloLocation_(entity_kind),
                          NULL);
      AMANZI_ASSERT(!ierr);
      closeFile_();
    }
    mesh_->getComm()->barrier();
  }

  if (comm_rank == 0) {
    fid_ = DBOpen(fname.c_str(), DB_HDF5, DB_APPEND);
    DBSetDir(fid_, "/");
    std::vector<std::string> varnames_str(comm_size);
    std::vector<char*> varnames(comm_size, nullptr);
    std::vector<int> vartypes(comm_size);
    for (int i = 0; i != comm_size; ++i) {
      varnames_str[i] = std::string("/domain_") + std::to_string(i) + "/" + varname;
      varnames[i] = const_cast<char*>(varnames_str[i].c_str());
      vartypes[i] = DB_UCDVAR;
    }

    int ierrl = DBPutMultivar(fid_,
                              fixName_(attrs.name()).c_str(),
                              comm_size,
                              (char**)varnames.data(),
                              (int*)vartypes.data(),
                              nullptr);
    AMANZI_ASSERT(!ierrl);
    closeFile_();
  }
}


} // namespace Amanzi

#endif
