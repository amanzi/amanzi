/*
  Output

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Silo implementation of an Output object.
*/

#ifndef AMANZI_OUTPUT_SILO_HH_
#define AMANZI_OUTPUT_SILO_HH_

#include <string>
#include <vector>

#include "Teuchos_ParameterList.hpp"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"

extern "C" {
#include "silo.h"
};

#include "errors.hh"
#include "dbc.hh"
#include "Mesh.hh"

#include "Output.hh"

namespace Amanzi {

class OutputSilo : public Output {

 public:

  OutputSilo(Teuchos::ParameterList& plist,
             const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
             bool is_vis,
             bool is_dynamic);

  // destructor must release file resource on non-finalized
  ~OutputSilo();
  
  // open and close files
  virtual void InitializeCycle(double time, int cycle, const std::string& tag);
  virtual void FinalizeCycle();

  // write data to file
  virtual void WriteVector(const Epetra_Vector& vec, const std::string& name,
                           const AmanziMesh::Entity_kind& kind) const;
  virtual void WriteMultiVector(const Epetra_MultiVector& vec, const std::vector<std::string>& names,
                                const AmanziMesh::Entity_kind& kind) const;

  // can we template this?
  virtual void WriteAttribute(const double& val, const std::string& name) const;
  virtual void WriteAttribute(const int& val, const std::string& name) const;
  virtual void WriteAttribute(const std::string& val, const std::string& name) const;

  // read data from file
  virtual void ReadVector(Epetra_Vector& vec, const std::string& name) const {
    ReadThrowsError_();
  }
      
  virtual void ReadMultiVector(Epetra_MultiVector& vec,
          const std::vector<std::string>& name) const {
    ReadThrowsError_();
  }

  virtual void ReadAttribute(double& val, const std::string& name) const {
    ReadThrowsError_();
  }
    
  virtual void ReadAttribute(int& val, const std::string& name) const {
    ReadThrowsError_();
  }
    
  virtual void ReadAttribute(std::string& val, const std::string& name) const {
    ReadThrowsError_();
  }

 protected:
  void Init_(Teuchos::ParameterList& plist);
  void ReadThrowsError_() const;
  void CloseFile_()  const;
  void WriteMesh_();
  std::string FixName_(const std::string& instring) const;
  
 protected:

  std::string filenamebase_;
  int count_;
  int sigfigs_;
  
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  mutable DBfile* fid_;
};
  
} // namespace Amanzi

#endif
