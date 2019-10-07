#include "CombinativeTreeOperator.hh"


namespace Amanzi {
namespace Operators {


CombinativeTreeOperator::CombinativeTreeOperator(const Teuchos::RCP<TreeOperator>& tree_op, 
                                                 Teuchos::RCP<Teuchos::ParameterList> plist, bool use_asc):
  tree_op_(tree_op)
{
  block_id_ = plist->get<Teuchos::Array<int> >("correction blocks").toVector();
  global_solve_ = plist->get<bool>("global solve", true);
  use_asc_ = use_asc;
}

void CombinativeTreeOperator::InitPreconditionerGlobal(const std::string& prec_name, 
        const Teuchos::ParameterList& plist) {
  tree_op_->InitPreconditioner(prec_name, plist);
}

void CombinativeTreeOperator::SymbolicAssembleMatrix() {
  tree_op_->SymbolicAssembleMatrix();
  for (std::vector<int>::const_iterator it = block_id_.begin(); it != block_id_.end(); it++) {
    auto block = &*tree_op_->GetBlock(*it,*it);
    const_cast<Operator*>(block)->SymbolicAssembleMatrix();
  }
}

void CombinativeTreeOperator::AssembleMatrix() {
  tree_op_->AssembleMatrix();
  for (std::vector<int>::const_iterator it = block_id_.begin(); it != block_id_.end(); it++) {
    auto block = &*tree_op_->GetBlock(*it,*it);
    const_cast<Operator*>(block)->AssembleMatrix();
  }
}

void CombinativeTreeOperator::InitPreconditionerBlock(const std::vector<std::string>& prec_names, 
        const Teuchos::ParameterList& plist) {
  for (std::vector<int>::const_iterator it = block_id_.begin(); it != block_id_.end(); it++) {
    auto block = &*tree_op_->GetBlock(*it,*it);
    const_cast<Operator*>(block)->InitPreconditioner(prec_names[*it], plist);
  }
}

int CombinativeTreeOperator::Apply(const TreeVector& X, TreeVector& Y) const {
  return tree_op_->Apply(X, Y);
}

// Y = inv(A) * X
int CombinativeTreeOperator::ApplyInverse(const TreeVector& X, TreeVector& Y) const
{
  if (!use_asc_) {

    Teuchos::RCP<TreeVector> tmp_rhs = Teuchos::rcp(new TreeVector(X));

    if (global_solve_) {
      // obtain u_{k+1/2} = x_k + ILU(A)^-1 * r_k
      tree_op_->ApplyInverse(X, Y);

      // compute the residual r_{k+1/2} = A * u_{k+1/2}
      tmp_rhs->PutScalar(0.0);
      tree_op_->Apply(Y, *tmp_rhs);
      tmp_rhs->Update(1.0, X, -1.0);
    }

    // compute the correction(i) = AMG(Aii)^-1 * r_{k+1/2}(i)
    Teuchos::RCP<TreeVector> tmp_pu = Teuchos::rcp(new TreeVector(Y));
    tmp_pu->PutScalar(0.0);
    for (std::vector<int>::const_iterator it = block_id_.begin(); it != block_id_.end(); it++) {
      tree_op_->GetBlock(*it,*it)->ApplyInverse(*tmp_rhs->SubVector(*it)->Data(), *tmp_pu->SubVector(*it)->Data());
    }

    // update x_k = u_{k+1/2} + correction
    for (std::vector<int>::const_iterator it = block_id_.begin(); it != block_id_.end(); it++) {
      Y.SubVector(*it)->Data()->Update(1.0, *tmp_pu->SubVector(*it)->Data(), 1.0);
    }

    /*
    Y.PutScalar(0.0);
    Teuchos::RCP<TreeVector> tmp_rhs = Teuchos::rcp(new TreeVector(X));

    // compute the correction(i) = AMG(Aii)^-1 * rhs(i)
    Teuchos::RCP<TreeVector> tmp_pu = Teuchos::rcp(new TreeVector(Y));
    for (std::vector<int>::const_iterator it = block_id_.begin(); it != block_id_.end(); it++) {
      tree_op_->GetBlock(*it,*it)->ApplyInverse(*tmp_rhs->SubVector(*it)->Data(), *tmp_pu->SubVector(*it)->Data());
    }

    if (global_solve_) {
      // compute the residual r_{k+1/2} = rhs - A*correction(i)
      tree_op_->Apply(*tmp_pu, *tmp_rhs);
      tmp_rhs->Update(1.0, X, -1.0);

      // obtain u_{k+1/2} = u_k + ILU(A)^-1 * r_{k+1/2}
      tree_op_->ApplyInverse(*tmp_rhs, Y);
    }

    // update x_k = u_{k+1/2} + correction
    for (std::vector<int>::const_iterator it = block_id_.begin(); it != block_id_.end(); it++) {
      Y.SubVector(*it)->Data()->Update(1.0, *tmp_pu->SubVector(*it)->Data(), 1.0);
    }
    */
  }

  /***************************************************
  *        Testing Schur Complement method           *
  ****************************************************/
  else {
    Y.PutScalar(0.0);
    // compute delta_s = inv(A_ss) * rhs_s
    //X.SubVector(1)->Data()->Print(std::cout);
    tree_op_->GetBlock(1,1)->ApplyInverse(*X.SubVector(1)->Data(), *Y.SubVector(1)->Data());
    //Y.SubVector(1)->Data()->Print(std::cout);

    // compute rhs_p = rhs_p - A_ps * delta_s
    Teuchos::RCP<TreeVector> tmp_rhs = Teuchos::rcp(new TreeVector(*Y.SubVector(1)));
    tmp_rhs->PutScalar(0.0);
    tree_op_->GetBlock(0,1)->Apply(*Y.SubVector(1)->Data(), *tmp_rhs->Data());
    tmp_rhs->Update(1.0, *X.SubVector(0), -1.0);

    // compute delta_p = inv(S) * rhs_p where S = A_pp - A_ps * inv(diag(A_ss)) * A_sp
    tree_op_->GetBlock(0,0)->ApplyInverse(*tmp_rhs->Data(), *Y.SubVector(0)->Data());
    //Y.SubVector(0)->Data()->Print(std::cout);
  }

  return 0; 
}

}
}
