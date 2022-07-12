/*
    This is the partition fn base class

    This will provide an interface to access slow, fast and full solutions and time derivatives

    To be used with both IMIM and EXIM schemes


*/

#ifndef AMANZI_PARTITION_FNBASE_HH_
#define AMANZI_PARTITION_FNBASE_HH_

template<class Vector>
class PartitionFnBase
{
public:
  //Modify the FULL solution
  // TODO:
  virtual void ModifySolutionFull(const double t, Teuchos::RCP<Vector> u) {};

  // computes FULL functional f = f(t,u_1, u_2, ..., u_n)
  // TODO:
  virtual void FunctionalTimeDerivativeFull(const double t, const Teuchos::RCP<Vector> u, Teuchos::RCP<Vector> u_eval) = 0;


  //Modify the Slow of the solution
  // TODO:
  virtual void ModifySolutionSlow(const double t, Teuchos::RCP<Vector> u) {};

  // Computes the Slow RHS ys' = f^{s}(t, y)
  // TODO:
  virtual void FunctionalTimeDerivativeSlow(const double t, const Teuchos::RCP<Vector> u, Teuchos::RCP<Vector> u_eval) = 0;

  //Modify the Slow of the solution
    //TODO:
    virtual void ModifySolutionFast(const double t, Teuchos::RCP<Vector> u) {};

    // Computes the Slow RHS yf' = f^{f}(t, y)
    // TODO:
    virtual void FunctionalTimeDerivativeFast(const double t, const Teuchos::RCP<Vector> u, Teuchos::RCP<Vector> u_eval) = 0;

    
  // computes a norm on u-du and returns the result
  virtual double ErrorNorm(Teuchos::RCP<const Vector> u, Teuchos::RCP<const Vector> du) = 0;

};

#endif



