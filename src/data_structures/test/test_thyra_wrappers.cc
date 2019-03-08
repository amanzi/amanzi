#include <ios>
#include "UnitTest++.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraits.hpp"

#include "Epetra_MpiComm.h"
#include "Thyra_VectorBase.hpp"
#include "Thyra_VectorStdOps.hpp"

#include "MeshFactory.hh"
#include "CompositeVector.hh"
#include "TreeVector.hh"
#include "amanzi_thyra_wrappers.hh"

#include "dbc.hh"
#include "errors.hh"
#include "exceptions.hh"

using namespace Amanzi;

struct test_data {
  Comm_ptr_type comm;
  Teuchos::RCP<AmanziMesh::Mesh> mesh;

  Teuchos::RCP<TreeVector> tv1;
  Teuchos::RCP<TreeVector> tv2;
  Teuchos::RCP<CompositeVector> cv;

  test_Data() {
    comm = new Epetra_MpiComm(MPI_COMM_SELF);
    AmanziMesh::MeshFactory mesh_fact(comm);
    mesh = mesh_fact(0.0, 0.0, 0.0, 2.0, 1.0, 1.0, 2, 1, 1);

    std::vector<AmanziMesh::Entity_kind> locations(2);
    locations[0] = AmanziMesh::CELL;
    locations[1] = AmanziMesh::FACE;

    std::vector<std::string> names(2);
    names[0] = "cell";
    names[1] = "face";

    std::vector<int> num_dofs(2);
    num_dofs[0] = 2;
    num_dofs[1] = 1;

    cv = Teuchos::rcp(new CompositeVector(mesh, names, locations, num_dofs, true));
    cv->CreateData();
    cv->PutScalar(1.0);

    tv1 = Teuchos::rcp(new TreeVector("TV with data"));
    tv1->SetData(cv);

    tv2 = Teuchos::rcp(new TreeVector("TV with subvecs"));
    tv2->PushBack(tv1);
    tv2->PushBack(tv1);

  }

  ~test_Data() {  }
};


SUITE(AmanziThyraWrappersTests) {
  // data structures for testing
  TEST_FIXTURE(test_data, CreateThyraFromCV) {
    // assign a value to the CV and get norms
    double norm_cv = 0;
    int ierr = cv->Norm1(&norm_cv);
    CHECK(!ierr);

    // create the thyra object
    Teuchos::RCP< Thyra::VectorBase<double> > thyra_vec =
      ThyraWrappers::CreateThyraVector(cv);

    // check the thyra norm is close to the CV norm.
    //    Teuchos::ScalarTraits<double>::magnitudeType norm_th = Thyra::norm_1<double>(*thyra_vec);
    double norm_th = Thyra::norm_1<double>(*thyra_vec);
    CHECK_CLOSE(norm_th, norm_cv, 0.000001);

    // check the const version
    Teuchos::RCP<const CompositeVector> cv2 = Teuchos::rcp(new CompositeVector(*cv));
    Teuchos::RCP< const Thyra::VectorBase<double> > thyra_vec2 =
      ThyraWrappers::CreateThyraVector(cv2);
    double norm_cv2 = 0;
    ierr = cv2->Norm1(&norm_cv2);
    CHECK(!ierr);
    double norm_th2 = Thyra::norm_1<double>(*thyra_vec2);
    CHECK_CLOSE(norm_th2, norm_cv2, 0.000001);
  }

  TEST_FIXTURE(test_data, CreateThyraFromTVwData) {
    // assign a value to the CV and get norms
    double norm_tv = 0;
    int ierr = tv1->Norm1(&norm_tv);
    CHECK(!ierr);

    // create the thyra object
    Teuchos::RCP< Thyra::VectorBase<double> > thyra_vec =
      ThyraWrappers::CreateThyraVector(tv1);

    // check the thyra norm is close to the TV norm.
    double norm_th = Thyra::norm_1<double>(*thyra_vec);
    CHECK_CLOSE(norm_th, norm_tv, 0.000001);
  }

  TEST_FIXTURE(test_data, CreateThyraFromTVwSubVecs) {
    // assign a value to the CV and get norms
    double norm_tv = 0;
    int ierr = tv2->Norm1(&norm_tv);
    CHECK(!ierr);

    // create the thyra object
    Teuchos::RCP< Thyra::VectorBase<double> > thyra_vec =
      ThyraWrappers::CreateThyraVector(tv2);

    // check the thyra norm is close to the TV norm.
    double norm_th = Thyra::norm_1<double>(*thyra_vec);
    CHECK_CLOSE(norm_th, norm_tv, 0.000001);
  }

  TEST_FIXTURE(test_data, CreateCVFromThyra) {
    // create the thyra object
    Teuchos::RCP< Thyra::VectorBase<double> > thyra_vec =
      ThyraWrappers::CreateThyraVector(cv);

    // check the thyra norm is close to the CV norm.
    //    Teuchos::ScalarTraits<double>::magnitudeType norm_th = Thyra::norm_1<double>(*thyra_vec);
    double norm_th = Thyra::norm_1<double>(*thyra_vec);

    // convert back to a new CV.
    Teuchos::RCP<CompositeVector> cv2 =
      ThyraWrappers::CreateCompositeVector(thyra_vec, cv);

    // check norms
    double norm_cv = 0;
    int ierr = cv2->Norm1(&norm_cv);
    CHECK(!ierr);
    CHECK_CLOSE(norm_th, norm_cv, 0.000001);
  }

  TEST_FIXTURE(test_data, CreateTVFromThyra) {
    // create the thyra object
    Teuchos::RCP< Thyra::VectorBase<double> > thyra_vec =
      ThyraWrappers::CreateThyraVector(tv2);

    // check the thyra norm is close to the CV norm.
    //    Teuchos::ScalarTraits<double>::magnitudeType norm_th = Thyra::norm_1<double>(*thyra_vec);
    double norm_th = Thyra::norm_1<double>(*thyra_vec);

    // convert back to a new TV.
    Teuchos::RCP<TreeVector> tv3 =
      ThyraWrappers::CreateTreeVector("test", thyra_vec, tv2);
    double norm_tv = 0;
    int ierr = tv3->Norm1(&norm_tv);
    CHECK(!ierr);
    CHECK_CLOSE(norm_th, norm_tv, 0.000001);
  }

}
