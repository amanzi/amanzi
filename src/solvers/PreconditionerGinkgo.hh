/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.
*/

#ifndef AMANZI_PRECONDITIONER_GINKGO_HH_
#define AMANZI_PRECONDITIONER_GINKGO_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#ifdef HAVE_IFPACK2_GINKGO
#include "ginkgo/ginkgo.hpp"
#endif

#include "exceptions.hh"
#include "Preconditioner.hh"

namespace Amanzi {
namespace AmanziSolvers {

using Ifpack2_PC_type = Ifpack2::Preconditioner<Matrix_type::scalar_type,
                                               Matrix_type::local_ordinal_type,
                                               Matrix_type::global_ordinal_type,
                                               Matrix_type::node_type>;


class PreconditionerGinkgo : public Preconditioner {
 public:

#ifdef HAVE_IFPACK2_GINKGO
  // Ginkgo Matrix type 
  using SC = Matrix_type::scalar_type; 
  using LO = Matrix_type::local_ordinal_type; 
  using mtx = gko::matrix::Csr<SC,LO>;
  using vec = gko::matrix::Dense<SC>;

  // Preconditioners
  // Block Jacobi 
  using bj = gko::preconditioner::Jacobi<SC,LO>;
  using ic = gko::preconditioner::Ic<gko::solver::LowerTrs<>,LO>;
  using ilu = gko::preconditioner::Ilu<gko::solver::LowerTrs<>,gko::solver::UpperTrs<>,false,LO>;
  using isai = gko::preconditioner::Isai<gko::preconditioner::isai_type::general,SC,LO>;

  PreconditionerGinkgo():
    Preconditioner(), initialized_(false), exec_(gko::ReferenceExecutor::create())
    {}

  //std::shared_ptr<bj> 
  void 
  copyTpetraToGinkgo(){

    const auto nrows = h_->getNodeNumRows(); 
    // Compute the cols per rows directly on the device 
    Kokkos::View<LO*,Kokkos:: DefaultExecutionSpace> colsperrow("ColsPerRow",nrows+1); 
    auto rowPtrs = h_->getCrsGraph()->getLocalGraph().row_map; 
    Kokkos::parallel_for(nrows,
      KOKKOS_LAMBDA(const int i){
        colsperrow(i) = rowPtrs[i+1]-rowPtrs[i]; 
    }); 
    Kokkos::parallel_for(nrows,
      KOKKOS_LAMBDA(const int i){
        if(i != 0)
          colsperrow(i) = colsperrow(i)+colsperrow(i-1); 
    }); 
    Kokkos::parallel_for(nrows+1,
      KOKKOS_LAMBDA(const int i){
        if(i != nrows)
          colsperrow(nrows-i) = colsperrow(nrows-1-i);
        else 
          colsperrow(0) = 0; 
    }); 
    auto rowindices = h_->getRowMap()->getMyGlobalIndices();

    Kokkos::View<SC*,Kokkos:: DefaultExecutionSpace> values = h_->getLocalValuesView(); 
    Kokkos::View<LO*,Kokkos:: DefaultExecutionSpace> colindices = h_->getCrsGraph()->getLocalGraph().entries; 
    
    std::cout<<"nrows: "<<nrows<<std::endl;
    std::cout<<"values: "<<values.size()<<std::endl;
    std::cout<<"colindices: "<<colindices.size()<<std::endl;
    std::cout<<"colsperrow: "<<colsperrow.size()<<std::endl;

    matrix_ = mtx::create(
        exec_, gko::dim<2>(nrows,nrows),
        gko::Array<Matrix_type::scalar_type>::view(exec_, values.size(), values.data()),
        gko::Array<Matrix_type::local_ordinal_type>::view(exec_, colindices.size(), colindices.data()),
        gko::Array<Matrix_type::local_ordinal_type>::view(exec_, colsperrow.size(), colsperrow.data())); 
        //std::make_shared<typename Mtx::load_balance>(2));

    if(method_ == "jacobi"){
      prec_ = bj_factory_->generate(exec_,lend(matrix_));     
    }else if(method_ == "ic"){
      prec_ = ic_factory_->generate(exec_,lend(matrix_));     
    }else if(method_ == "ilu"){
      prec_ = ilu_factory_->generate(exec_,lend(matrix_));     
    }else if(method_ == "isai"){
      prec_ = isai_factory_->generate(exec_,lend(matrix_));     
    }else{
      Errors::Message msg("\"Ginkgo\" method unknown");
      Exceptions::amanzi_throw(msg);
    }
    
    #if 0 
    auto v = matrix_->get_values(); 
    std::cout<<"Matrix: "; 
    for(int i = 0 ; i <  values.size(); ++i){
      std::cout<<v[i]<<" ";
    }
    std::cout<<std::endl;
    #endif 

  }

  virtual void set_inverse_parameters(Teuchos::ParameterList& plist) override {
    plist_ = plist;
    std::string vo_name = this->name()+" ("+plist_.get<std::string>("method")+")";
    vo_ = Teuchos::rcp(new VerboseObject(vo_name, plist_));
    method_ = plist_.get<std::string>("method");
    initialized_ = true;
  }

  virtual void initializeInverse() override {
    std::cout<<"initializeInverse"<<std::endl;
    if(method_ == "jacobi"){
      bool skip_sorting = false; 
      std::cout<<"jacobi method"<<std::endl;
      //factory_ = bj::build().with_max_block_size(3).on(exec_); 
      bj_factory_ = bj::build() 
        .with_max_block_size(32u) // Default 32u
        .with_max_block_stride(0u) // Default 0u
        .with_skip_sorting(false) // Default False
        .with_block_pointers(nullptr) // Default nullptr
        .with_accuracy(1e-1) // Default 1e-1
        .on(exec_); 
    }else if(method_ == "ic"){
      std::cout<<"ICholevsky method"<<std::endl;
      ic_factory_ = ic::build() 
        .with_l_solver_factory(nullptr) // Default nullptr
        .with_factorization_factory(nullptr) // Default nullptr
        .on(exec_); 
    }else if(method_ == "ilu"){
      std::cout<<"Ilu method"<<std::endl;
      ilu_factory_ = ilu::build()
        .with_l_solver_factory(nullptr) // Default nullptr
        .with_u_solver_factory(nullptr) // Default nullptr
        .with_factorization_factory(nullptr) // Default nullptr
        .on(exec_); 
    }else if(method_ == "isai"){
      std::cout<<"Isai method"<<std::endl;
      isai_factory_ = isai::build()
        .with_skip_sorting(false) // Default false 
        .with_sparsity_power(1) // Default 1
        .with_excess_limit(0u) // Default 0u
        .with_excess_solver_factory(nullptr) // Default nullptr
        .on(exec_); 
    }else{
      Errors::Message msg("\"Ginkgo\" method unknown");
      Exceptions::amanzi_throw(msg);
    }

    std::cout<<"initializeInverse DONE"<<std::endl;
  }

  virtual void computeInverse() override {
    std::cout<<"ComputeInverse"<<std::endl;
    // Change matrix values 
    #if 0 
    Kokkos::View<SC*,Kokkos:: DefaultExecutionSpace> values = h_->getLocalValuesView();
    std::cout<<"New values: "<<values.size()<<std::endl; 

    auto v = matrix_->get_values(); 

    std::cout<<"Matrix: "; 
    for(int i = 0 ; i <  values.size(); ++i){
      std::cout<<v[i]<<" ";
    }
    std::cout<<std::endl;

    std::cout<<"Values: "; 
    for(int i = 0 ; i <  values.size(); ++i){
      std::cout<<values[i]<<" ";
    }
    std::cout<<std::endl;

    std::cout<<"Copy values in matrix"<<std::endl;
    for(int i = 0 ; values.size(); ++i){
      v[i] = values(i);
    }
    prec_ = bj_factory_->generate(exec_,lend(matrix_)); 
    #endif 

    copyTpetraToGinkgo(); 
    
    std::cout<<"ComputeInverse DONE"<<std::endl;
  }

  virtual int returned_code() const override final { return returned_code_; }

  virtual int applyInverse(const Vector_type& v, Vector_type& hv) const override {
    assert(v.getData().size() == hv.getData().size()); 
    // Convert vector to linOp: if contiguous just need pointer 
    auto lo_v = vec::create(exec_,gko::dim<2>(v.getData().size(),1));//,v.getData().getRawPtr()); 
    auto lo_hv = vec::create(exec_,gko::dim<2>(hv.getData().size(),1));//,hv.getData().getRawPtr()); 
    // Copy data 
    for(int i = 0 ; i < v.getData().size(); ++i){ lo_v->at(i,0) = v.getData()[i]; }

    prec_->apply(lend(lo_v),lend(lo_hv)); 
    
    // Copy result back 
    for(int i = 0 ; i < v.getData().size(); ++i){ hv.getDataNonConst()[i] = lo_hv->at(i,0); }
    return 0;
  }

  virtual std::string returned_code_string() const override
  { return "NULL"; }

 protected:
  Teuchos::ParameterList plist_;

  Teuchos::RCP<VerboseObject> vo_;

  std::shared_ptr<mtx> matrix_; 

  const std::shared_ptr<gko::ReferenceExecutor>  exec_;

  //std::unique_ptr<typename gko::AbstractFactory> factory_; 

  std::unique_ptr<typename bj::Factory> bj_factory_;
  std::unique_ptr<typename ic::Factory> ic_factory_;
  std::unique_ptr<typename ilu::Factory> ilu_factory_;
  std::unique_ptr<typename isai::Factory> isai_factory_;

  std::shared_ptr<gko::LinOp> prec_; 

  std::string method_;
  bool initialized_;
  mutable int returned_code_;
#endif
};


} // namespace AmanziPreconditioners
} // namespace Amanzi

#endif
