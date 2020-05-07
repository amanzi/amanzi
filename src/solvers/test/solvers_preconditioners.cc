/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

#include <iostream>
#include <string>

#ifdef _OPENMP
#  include "omp.h"
#endif

#define HAVE_EPETRA_PRECONDITIONERS
#define HAVE_TRILINOS_PRECONDITIONERS
#define HAVE_HYPRE_PRECONDITIONERS

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "UnitTest++.h"

#include "exceptions.hh"
#include "AmanziTypes.hh"
#include "matrix.hh"
#include "Preconditioner.hh"
#include "PreconditionerFactory.hh"
#include "LinearOperatorPCG.hh"

SUITE(PRECONDITIONERS)
{

  TEST(PC_ONLY)
  {
    std::cout << "One pass of preconditioner:" << std::endl
              << "------------------------------" << std::endl;

    auto main_list = Teuchos::rcp(new Teuchos::ParameterList());
    Teuchos::ParameterXMLFileReader xmlreader("test/solvers_preconditioners.xml");
    *main_list = xmlreader.getParameters();
    
    Matrix m(10);
    std::vector<std::string> prec_names = {"identity", "diagonal", "ifpack2: ILUT"};

    Vector_type x(m.getDomainMap());
    Vector_type y(m.getRangeMap());

    auto pc_list = Teuchos::sublist(main_list, "preconditioners");
    
    for (const auto& name : prec_names) {
      m.Init(name, pc_list);
      y.putScalar(1.);
      x.putScalar(0.);

      m.applyInverse(y,x);

      if (name == "identity") {
        x.update(-1., y, 1.);
        CHECK_CLOSE(0., x.normInf(), 1.e-10);
      } else if (name == "diagonal") {
        x.sync_host();
        auto x_v = x.getLocalViewHost();
        for (int i=0; i!=10; ++i) {
          CHECK_CLOSE(1.0 / (2*i+1), x_v(i,0), 1.e-10);
        }

      } else {
        CHECK(x.normInf() > 0.);
      }

      std::cout << name << ":" << std::endl;
      Teuchos::FancyOStream out(Teuchos::rcpFromRef(std::cout));
      x.describe(out, Teuchos::VERB_EXTREME);
    }
  }


  TEST(PC_PCG)
  {
    int N = 125;
    std::cout << "PCG + preconditioner:" << std::endl
              << "------------------------------" << std::endl;

    auto main_list = Teuchos::rcp(new Teuchos::ParameterList());
    Teuchos::ParameterXMLFileReader xmlreader("test/solvers_preconditioners.xml");
    *main_list = xmlreader.getParameters();
    
    auto m = Teuchos::rcp(new Matrix(N));
    std::vector<std::string> prec_names = {"identity", "diagonal", "ifpack2: ILUT"};

    // create the pcg operator
    AmanziSolvers::LinearOperatorPCG<Matrix, Vector_type, Map_type> pcg(m,m);
    pcg.Init();
    pcg.set_tolerance(1.e-12);
    pcg.set_max_itrs(200);

    
    Vector_type y(m->getRangeMap());
    {
      auto yv = y.getLocalViewHost();
      for (int i = 0; i < N; i++) yv(i,0) = 1.0 / (i + 2.0);
    }
    y.sync_device();
    
    Vector_type x(m->getDomainMap());

    auto pc_list = Teuchos::sublist(main_list, "preconditioners");
    
    for (const auto& name : prec_names) {
      m->Init(name, pc_list);
      x.putScalar(0.);

      pcg.applyInverse(y,x);

      {
        x.sync_host();
        auto x_v = x.getLocalViewHost();
        CHECK_CLOSE(11.03249773994628, x_v(0,0), 1e-6);
        CHECK_CLOSE(10.53249773994628, x_v(1,0), 1e-6);
      }

      if (Keys::startsWith(name, "ifpack2")) {
        CHECK_EQUAL(1, pcg.num_itrs());
      }
      
      std::cout << name << " : " << pcg.num_itrs() << std::endl;
      // Teuchos::FancyOStream out(Teuchos::rcpFromRef(std::cout));
      // x.describe(out, Teuchos::VERB_EXTREME);
    }
  }
}


