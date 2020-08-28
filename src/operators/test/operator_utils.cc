/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  Tests nonlinear, coupled diffusion problem.
*/


#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Epetra_Vector.h"
#include "EpetraExt_RowMatrixOut.h"
#include "UnitTest++.h"

// Amanzi
#include "Mesh_MSTK.hh"
#include "TreeVector.hh"
#include "SuperMap.hh"

// Operators
#include "OperatorDefs.hh"
#include "OperatorUtils.hh"

using namespace Teuchos;
using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;
using namespace Amanzi::Operators;

struct Maps {
  Maps() {
    comm = Amanzi::getDefaultComm();
    
    // create a mesh
    mesh = Teuchos::rcp(new Mesh_MSTK(0.,0.,1.,1.,10,10, comm));

    // create a vector
    cvs = Teuchos::rcp(new CompositeVectorSpace());
    cvs->SetMesh(mesh)->SetGhosted(true)
        ->AddComponent("cell", AmanziMesh::CELL, 1)
        ->AddComponent("face", AmanziMesh::FACE, 1);

    Teuchos::RCP<TreeVectorSpace> tvs0 = Teuchos::rcp(new TreeVectorSpace());
    tvs0->SetData(cvs);

    tvs = Teuchos::rcp(new TreeVectorSpace());
    tvs->PushBack(tvs0);
    tvs->PushBack(tvs0);

    // create a supermap, vec
    map = createSuperMap(*tvs);
    
  }

  Comm_ptr_type comm;
  Teuchos::RCP<Mesh> mesh;
  Teuchos::RCP<CompositeVectorSpace> cvs;
  Teuchos::RCP<TreeVectorSpace> tvs;
  Teuchos::RCP<SuperMap> map;
};


TEST(SUPERMAP_COPY_INVERTIBLE) {
  Maps maps;
  Teuchos::RCP<TreeVector> tv = Teuchos::rcp(new TreeVector(*maps.tvs));

  // initialize randomly
  tv->SubVector(0)->Data()->Random();
  tv->SubVector(1)->Data()->Random();

  Epetra_Vector vec(*maps.map->Map());

  // copy forward, backward
  Teuchos::RCP<TreeVector> tv2 = Teuchos::rcp(new TreeVector(*maps.tvs));
  int ierr = copyToSuperVector(*maps.map, *tv, vec);
  CHECK(!ierr);

  ierr = copyFromSuperVector(*maps.map, vec, *tv2);
  CHECK(!ierr);

  // check the same
  tv2->Update(-1., *tv, 1.);
  double norm;
  tv2->Norm2(&norm);
  CHECK_CLOSE(0., norm, 1.e-16);
}


TEST(SUPERMAP_COPY_INTS) {

  Maps maps;
  Teuchos::RCP<TreeVector> tv = Teuchos::rcp(new TreeVector(*maps.tvs));
  tv->SubVector(0)->Data()->ViewComponent("cell", false)->PutScalar(3.);
  tv->SubVector(1)->Data()->ViewComponent("cell", false)->PutScalar(4.);
  tv->SubVector(0)->Data()->ViewComponent("face", false)->PutScalar(5.);
  tv->SubVector(1)->Data()->ViewComponent("face", false)->PutScalar(6.);

  Epetra_Vector vec(*maps.map->Map());

  // copy forward
  Teuchos::RCP<TreeVector> tv2 = Teuchos::rcp(new TreeVector(*maps.tvs));
  int ierr = copyToSuperVector(*maps.map, *tv, vec);
  CHECK(!ierr);

  // check values
  int ncells = maps.mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces = maps.mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);

  // check sizes
  CHECK_EQUAL(2*ncells + 2*nfaces, vec.MyLength());

  for (int i = 0; i!=ncells; ++i) {
    CHECK_EQUAL(3., vec[i*2]);
    CHECK_EQUAL(4., vec[i*2+1]);
  }

  for (int i = 0; i!=nfaces; ++i) {
    CHECK_EQUAL(5., vec[2*ncells + i*2]);
    CHECK_EQUAL(6., vec[2*ncells + i*2+1]);
  }
  
}
