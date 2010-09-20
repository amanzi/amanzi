#ifndef STKMESH_HPP
#define STKMESH_HPP

#include <vector>

#include "Teuchos_ParameterList.hpp"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"

// Start of STK stuff
#include <stk_util/parallel/Parallel.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/fem/FieldTraits.hpp>



class STKMesh1D {

public:
  typedef stk::mesh::Field<double,stk::mesh::Cartesian> VectorFieldType ;
  typedef stk::mesh::Field<double>                      ScalarFieldType ;  

  STKMesh1D(const Teuchos::RCP<const Epetra_Comm>& epetra_comm,
	    const Teuchos::RCP<Teuchos::ParameterList>& params);
  
  ~STKMesh1D();

  Teuchos::RCP<const Teuchos::ParameterList> 
  getValidDiscretizationParameters() const;
  
  stk::mesh::MetaData* metaData;
  stk::mesh::BulkData* bulkData;
  std::vector<stk::mesh::Part*> partVec;
  std::vector<stk::mesh::Part*> nsPartVec;
  VectorFieldType* coordinates_field;
};

#endif
