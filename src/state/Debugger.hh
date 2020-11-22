/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/
//! A mesh and vector structure aware utility for printing info.

/*!

This is a utility that makes it easier for the user to control output written
to the screen.  It allows the user to provide element IDs, and then provides
functionality for a PK to write mesh geometry information and vector values of
those elements to screen based upon verbosity levels.

Note, most information is only written if the owning object's verbosity level
from the `"Verbose Object`" spec is set to `"high`" or higher.

.. debugger-spec:
.. admonition:: debugger-spec

    * `"debug cells`" ``[Array(int)]`` For each global ID of a cell provided
      here, controls writing of vectors inside of the using PK.

    * `"debug faces`" ``[Array(int)]`` For each global ID of a face provided
      here, writes all adjoining cell information as if each cell was included
      in `"debug cells`".

*/
      
#ifndef AMANZI_DEBUGGER_HH_
#define AMANZI_DEBUGGER_HH_

#include <string>

#include "Teuchos_Ptr.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_MultiVector.h"

#include "VerboseObject.hh"
#include "Mesh.hh"


namespace Amanzi {

class CompositeVector;

class Debugger {

 public:
  // Constructor
  Debugger(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
           std::string name,
           Teuchos::ParameterList& plist,
           Teuchos::EVerbosityLevel verb_level=Teuchos::VERB_HIGH
           );

  // Write cell + face info
  void WriteCellInfo(bool include_faces=false);

  // Write a vector individually.
  void WriteVector(const std::string& name,
                   const Teuchos::Ptr<const CompositeVector>& vec,
                   bool include_faces=false);

  void WriteCellVector(const std::string& name,
                       const Epetra_MultiVector& vec);
  
  // Write boundary condition data.
  void WriteBoundaryConditions(const std::vector<int>& flag,
			       const std::vector<double>& data);

  // Write list of vectors.
  void WriteVectors(const std::vector<std::string>& names,
                    const std::vector<Teuchos::Ptr<const CompositeVector> >& vecs,
                    bool include_faces=false);

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
  std::vector<Teuchos::RCP<VerboseObject> > dcvo_;

  int width_;
  int header_width_;
  int precision_;
  int cellnum_width_;
  int decimal_width_;
};

} // namespace Amanzi


#endif
