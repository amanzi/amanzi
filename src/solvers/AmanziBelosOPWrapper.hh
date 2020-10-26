/*
  Solvers

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Alicia Klinvex (amklinv@sandia.gov)

  This tells Belos how to work with the Amanzi operators.
*/

#ifndef AMANZI_BELOS_OP_WRAPPER_HH_
#define AMANZI_BELOS_OP_WRAPPER_HH_

#include<BelosOperator.hpp>
#include<AmanziBelosMVWrapper.hh>

namespace Amanzi
{
// This is not an optimal implementation, but we can improve it later
template<class Op, class Vector>
class AmanziBelosOp : public Belos::Operator<double>
{
private:
  typedef double                                                    ScalarType;
  typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;
public:
  // Constructor
  AmanziBelosOp(Teuchos::RCP<const Op> op, bool applyInverse = false)
  { op_ = op; applyInverse_ = applyInverse; }

  // Apply method
  void Apply(const Belos::MultiVec<ScalarType> &x, Belos::MultiVec<ScalarType> &y,
      Belos::ETrans trans=Belos::NOTRANS) const
  {
    const CompositeMultiVector<Vector> *cmx = dynamic_cast<const CompositeMultiVector<Vector>*>(&x);
    CompositeMultiVector<Vector> *cmy = dynamic_cast<CompositeMultiVector<Vector>*>(&y);

    int nvecs = cmx->GetNumberVecs();
    for(int i=0; i<nvecs; i++)
    {
      Teuchos::RCP<Vector> singleX = cmx->getVector(i);
      Teuchos::RCP<Vector> singleY = cmy->getVector(i);
      if (applyInverse_)
        op_->applyInverse(*singleX,*singleY);
      else
        op_->apply(*singleX,*singleY);
    }
  }

  bool hasApplyTranspose() const { return false; }

private:
  Teuchos::RCP<const Op> op_;
  bool applyInverse_;
};

}

#endif
