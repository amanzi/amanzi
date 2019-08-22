#include <iostream>
#include <string>

#define HAVE_TPETRA_PRECONDITIONERS

#include "UnitTest++.h"

#include "exceptions.hh"
#include "Teuchos_RCP.hpp"

#include "AmanziComm.hh"
#include "VectorWrapper.hh"
#include "LinearOperatorPCG.hh"
#include "PreconditionerFactory.hh"
#include "PreconditionerDiagonal.hh"
#include "PreconditionerIdentity.hh"

SUITE(SOLVERS) {

const int N = 125;
using namespace Amanzi;
using namespace Amanzi::AmanziPreconditioners;

class Matrix {
  using CrsMatrix_type = Tpetra::CrsMatrix<double,int,int>;
  using Vector_type = VectorWrapper<Amanzi::Vector_type>;
  using cVector_type = VectorWrapper<const Amanzi::Vector_type>;
  
 public:
  Matrix() {}
  Matrix(const Map_ptr_type& map) : map_(map) {};
  ~Matrix() {};
  Matrix(const Matrix& other) : map_(other.map_) {}

  void Init(const std::string& name) {
    Teuchos::ParameterList plist;
    plist.set<std::string>("preconditioner type", name);
    std::string params(name);
    Teuchos::ParameterList& tmp = plist.sublist(params.append(" parameters"));

    if (name == "diagonal") {
      preconditioner_ = Teuchos::rcp(new PreconditionerDiagonal<Matrix_type,Vector_type>());
    } else if (name == "identity") {
      preconditioner_ = Teuchos::rcp(new PreconditionerIdentity<Matrix_type,Vector_type>());
    } else if (name == "ml") {
      PreconditionerFactory<Matrix_type,Vector_type> factory;
      tmp.set<int>("coarse: max size", 5);
      tmp.set<int>("cycle applications", 2);
      tmp.set<int>("ML output", 0);
      preconditioner_ = factory.Create(plist);
    } else {
      tmp.set<int>("max coarse size", 5);
      tmp.set<int>("cycle applications", 1);
      tmp.set<int>("verbosity", 0);
      PreconditionerFactory<Matrix_type,Vector_type> factory;
      preconditioner_ = factory.Create(plist);
    }
    A_ = Teuchos::rcp(new CrsMatrix_type(map_, map_, 3));
    for (int i = 0; i < N; i++) {
      int indices[3];
      double values[3] = {double(-i), double(2 * i + 1), double(-i - 1)};
      for (int k = 0; k < 3; k++) indices[k] = i + k - 1; 
      A_->insertLocalValues(i, 3, values, indices);
    }
    A_->fillComplete(map_, map_);
    preconditioner_->Update(A_);
  };    

  virtual int Apply(const Vector_type& v, Vector_type& mv) const { 
    std::cout << "   in apply: ";
    double norm(0.);
    v.Norm2(&norm);
    std::cout << "v = " << norm;

    A_->apply(*v.get(), *mv.get());

    norm = 0.;
    mv.Norm2(&norm);
    std::cout << " Av = " << norm << std::endl;
    return 0;
  }
  virtual int ApplyInverse(const Vector_type& v, Vector_type& hv) const {
    std::cout << "   in apply: ";
    double norm(0.);
    v.Norm2(&norm);
    std::cout << "v = " << norm;

    auto ret = preconditioner_->ApplyInverse(v, hv);
    
    norm = 0.;
    hv.Norm2(&norm);
    std::cout << " A^-1 v = " << norm << std::endl;
    return ret;
  }

  virtual const Map_ptr_type& DomainMap() const { return map_; }
  virtual const Map_ptr_type& RangeMap() const { return map_; }

 private:
  Map_ptr_type map_;
  Teuchos::RCP<CrsMatrix_type> A_;
  Teuchos::RCP<Preconditioner<Matrix_type,Vector_type> > preconditioner_;
};


TEST(DIAGONAL_PRECONDITIONER) {
  std::cout << "\nComparison of preconditioners for N=125" << std::endl;

  auto comm = Amanzi::getDefaultComm();
  auto map = Teuchos::rcp(new Map_type(N, 0, comm));

  // create the pcg operator
  Teuchos::RCP<Matrix> m = Teuchos::rcp(new Matrix(map));
  AmanziSolvers::LinearOperatorPCG<Matrix, VectorWrapper<Vector_type>, Map_ptr_type> pcg(m, m);
  pcg.Init();
  pcg.set_tolerance(1e-12);
  pcg.set_max_itrs(200);

  Vector_type u(map), v(map);
  auto u_ptr = Teuchos::rcpFromRef<Vector_type>(u);
  auto v_ptr = Teuchos::rcpFromRef<Vector_type>(v);
  VectorWrapper<Vector_type> vec_u(u_ptr);
  VectorWrapper<Vector_type> vec_v(v_ptr);
  
  {
    auto uv = u.getLocalViewHost();
    for (int i = 0; i < N; i++) uv(i,0) = 1.0 / (i + 2.0);
  }
  u.sync_device();

  // solving with preconditioner
  std::string prec_names[4];
  prec_names[0] = "identity";
  prec_names[1] = "diagonal";
  for (int n = 0; n < 4; n++) {
    m->Init(prec_names[n]);

    vec_v.PutScalar(0.0);
    printf("Preconditioner: %s\n", prec_names[n].c_str());

    double norms_vec = 0;
    vec_u.Norm2(&norms_vec);
    printf("  initial u norm: %g\n", norms_vec);

    pcg.ApplyInverse(vec_u, vec_v);

    v.sync_host();
    {
      auto vv = v.getLocalViewHost();
      CHECK_CLOSE(11.03249773994628, vv(0,0), 1e-6);
      CHECK_CLOSE(10.53249773994628, vv(1,0), 1e-6);
    }
  }

};

}




