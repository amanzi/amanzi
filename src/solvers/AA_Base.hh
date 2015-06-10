/*

*/

#ifndef AMANZI_AA_BASE_HH_
#define AMANZI_AA_BASE_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "VerboseObject.hh"
#include "dbc.hh"

namespace Amanzi {
namespace AmanziSolvers {

#define AA_TRUE 1
#define AA_FALSE 0
#define AA_EOL -1

template<class Vector, class VectorSpace>
class AA_Base {
 public:
  //AA_Base(int mvec, double vtol, const VectorSpace& map);
  AA_Base(int mvec, double vtol, double beta, const VectorSpace& map);
  ~AA_Base();

  void Init(Teuchos::ParameterList& plist) {
    vo_ = Teuchos::rcp(new VerboseObject("AA_Base", plist));
  }

  void Relax();
  void Restart();
  void Correction(const Vector&, Vector&,
                  const Teuchos::Ptr<const Vector>& oldv=Teuchos::null);

 private:
  int subspace_;  // boolean: a nonempty subspace
  int pending_;   // contains pending vectors -- boolean
  int mvec_;      // maximum number of subspace vectors
  double vtol_;   // vector drop tolerance
  double beta_;   // relaxation parameter

  Teuchos::RCP<Vector> *u_;  // previous iterates
  Teuchos::RCP<Vector> *dF_;  // preconditioner residuals
  Teuchos::RCP<Vector> *F_test;  // preconditioner residuals

  double **h_;  // matrix of w vector inner products 

  // Linked-list organization of the vector storage.
  int first_f_;  // index of first_f subspace vector
  int last_f_;   // index of last_f subspace vector
  int free_f_;   // index of the initial vector in free storage linked list
  int new_f_;
  int num_vec_;
  int *next_f_;  // next_f index link field
  int *prev_f_;  // previous index link field in doubly-linked subspace f

  Teuchos::RCP<VerboseObject> vo_;
};


/* ******************************************************************
 * Allocate memory
 ***************************************************************** */
template<class Vector, class VectorSpace>
AA_Base<Vector, VectorSpace>::AA_Base(int mvec, double vtol, double beta, const VectorSpace& map)
{
  mvec_ = std::max(mvec, 1);  // we cannot have mvec_ < 1
  vtol_ = vtol;
  beta_ = beta;

  u_ = new Teuchos::RCP<Vector> [mvec_ + 1];
  dF_ = new Teuchos::RCP<Vector> [mvec_ + 1];
  F_test = new Teuchos::RCP<Vector> [mvec_ + 1];

  for (int i = 0; i < mvec_ + 1; i++) {
    u_[i] = Teuchos::rcp(new Vector(map));
    dF_[i] = Teuchos::rcp(new Vector(map));
    F_test[i] = Teuchos::rcp(new Vector(map));
  }

  h_ = new double* [mvec_];
  for (int j = 0; j < mvec_; j++) {
     h_[j] = new double[mvec_];
  }


  // next_f_ = new int [mvec_ + 1];
  // prev_f_ = new int [mvec_ + 1];
};






/* ******************************************************************
 * Destroy memory
 ***************************************************************** */
template<class Vector, class VectorSpace>
AA_Base<Vector, VectorSpace>::~AA_Base()
{
  delete [] u_;
  delete [] dF_;
  for (int j = 0; j <mvec_; ++j) {
    delete [] h_[j];
  }
  delete [] h_;
  // delete [] next_f_;
  // delete [] prev_f_;
};


// /* ******************************************************************
//  * TBW
//  ***************************************************************** */
template<class Vector, class VectorSpace>
void AA_Base<Vector, VectorSpace>::Relax()
{
  // if (pending_) {
  //   // Drop the initial slot where the pending_ vectors are stored.
  //   assert(first_f_ >= 0);

  //   int new_f = first_f_;
  //   first_f_ = next_f_[first_f_];
  //   if (first_f_ == AA_EOL) {
  //     last_f_ = AA_EOL;
  //   } else {
  //     prev_f_[first_f_] = AA_EOL;
  //   }

  //   // Update the free storage list.
  //   next_f_[new_f] = free_f_;
  //   free_f_ = new_f;
  //   pending_ = AA_FALSE;
  // }
}


// /* ******************************************************************
//  * TBW
//  ***************************************************************** */
template<class Vector, class VectorSpace>
void AA_Base<Vector, VectorSpace>::Restart()
{
  // No vectors are stored.
  first_f_  = AA_EOL;
  last_f_   = AA_EOL;
  subspace_ = AA_FALSE;
  pending_  = AA_FALSE;

  // // Initialize the free storage linked list.
  // free_f_ = 0;
  // for (int k = 0; k < mvec_; k++) {
  //   next_f_[k] = AA_EOL;
  // }
  // next_f_[mvec_] = AA_EOL;

  new_f_ = 0;
  num_vec_ = 0;

}


/* ******************************************************************
 * TBW
 ***************************************************************** */
template<class Vector, class VectorSpace>
void AA_Base<Vector, VectorSpace>::Correction(const Vector& f, Vector &dir,
        const Teuchos::Ptr<const Vector>& u_old)
{
  int i, j, k, nvec, new_f;
  double s, hkk, hkj, cj;
  double  *hk, *hj, *c;

  Teuchos::RCP<Vector> vp, wp;
  Teuchos::RCP<Vector> ff = Teuchos::rcp(new Vector(f));
  Teuchos::RCP<Vector> fun = Teuchos::rcp(new Vector(f));

  std::cout.precision(12);

  std::cout<<"mvec "<<mvec_<<"\n";
  std::cout<<"new "<<new_f_<<" first "<<first_f_<<" last "<<last_f_<<"\n";
  std::cout<<"\n";

  // std::cout<<"F\n"; ff->Print(std::cout);
  // std::cout<<"u\n"; u_old->Print(std::cout);

  assert(new_f_ >= 0);
  // Save the accelerated correction for the next_f_ call.

  //ff->Scale(-1.);

  *dF_[new_f_] = *ff;
  *F_test[new_f_] = *ff;
  if (u_old != Teuchos::null) *u_[new_f_] = *u_old; 

  if (last_f_ >= 0){
    k = last_f_;
    dF_[last_f_] -> Update(1., *dF_[new_f_], -1.); 

    while (k != first_f_) {
      k--;
      if (k < 0) k = mvec_;
      //std::cout<<"k "<<k<<"\n";
      dF_[k]->Update(1., *dF_[last_f_], 1. );
      if (k == first_f_) break;
    }      
  }

  double *b  = new double[num_vec_];
  // double *d  = new double[num_vec_];
  // double **l  = new double*[num_vec_];
  double *th = new double[num_vec_];
  double *th_test = new double[num_vec_];
  double *y = new double[num_vec_];
  int *row = new int[num_vec_];
  double alp;

  for (int i=0; i<num_vec_; i++){
    th[i]=0.;
    row[i] = i;
    // l[i] = new double[num_vec_];
    // for (int j=0;j<num_vec_;j++) l[i][j] = 0.;
  }

  // std::cout<<"NEW\n";
  // F_test[new_f_]->Print(std::cout);

  // /// Construct system of normal equations
  for (i=0; i<num_vec_; i++){
    int I = (first_f_ + i)%(mvec_ + 1);
    //    F_test[I]->Print(std::cout);
    for (j=i; j<num_vec_; j++){
      int J = (first_f_ + j)%(mvec_ + 1);
      dF_[I]->Dot( *dF_[J], &h_[i][j]);
      h_[j][i] = h_[i][j];
    }
    //    d[i] = h_[i][i];
    for (j=0; j<num_vec_; j++) std::cout<<h_[i][j]<<" ";
    dF_[I]->Dot(*dF_[new_f_], &b[i]);
    std::cout<<"  =  "<<b[i]<<"\n";
  }

  bool singular = false;
  double eps = 1e-16;

  // /// Choletsky decomposition
  for (int m=0; m<num_vec_; m++){
    double max_val = 0;
    int i_max_val;
    for (int i=m;i<num_vec_;i++){
      if (h_[row[i]][m] > max_val){
        i_max_val = i;
        max_val = h_[row[i]][m];
      }
    }
    if (max_val< eps) {
      h_[row[m]][m] = 0.;
      singular = true;
      break;
    }
    // int k = row[m];
    // row[m] = i_max_val;
    // row[i_max_val] = k;

    std::cout<<"h_[row[m]][m] "<<h_[row[m]][m]<<"\n";
    h_[row[m]][m] = sqrt(h_[row[m]][m]);
    for (int j=m+1; j<num_vec_; j++){
      h_[row[m]][j] = h_[row[m]][j]/ h_[row[m]][m];
    }
    for (int i=m+1; i<num_vec_; i++){
      for (int j=i; j<num_vec_; j++){
        h_[row[i]][j] -= h_[row[m]][j]*h_[row[m]][i];
        //        std::cout<<h_[row[i]][j]<<" ";
      }       
      //std::cout<<"\n";
    }                  

    // l[m][row[m]] = sqrt(d[row[m]]);

    // for (int i=m+1; i<num_vec_;i++){
    //   l[m][row[i]] = h_[row[m]][row[i]];       
    //   for (int j=0;j<m;j++) {
    //     l[m][row[i]] -= l[j][row[m]]*l[j][row[i]];
    //   }
    //   l[m][row[i]] /=  l[m][row[m]];
    //   d[row[i]] = d[row[i]] - l[m][row[m]]*l[m][row[i]];
    // }    
  }
                                                
  std::cout<<"LLLLL\n";
  for (int i=0;i<num_vec_;i++) std::cout<<row[i]<<" ";std::cout<<"\n";
  
  for (int i=0;i<num_vec_;i++){
    for (int j=0; j<num_vec_; j++) std::cout<<h_[i][j]<<" ";
    std::cout<<"\n";
  }

  //if (num_vec_ > 2) exit(0);

  
   

  //if (!singular){
  for (int i=0; i<num_vec_; i++){
    y[row[i]] = b[row[i]];
    for (int j=0; j<i; j++) y[row[i]] -= h_[row[j]][i]*y[row[j]];
    if (h_[row[i]][i] > eps)
      y[row[i]] /= h_[row[i]][i];
    else
      y[row[i]] = 0.;
  }

  // std::cout<<"y val\n";
  // //  for (int i=0;i<num_vec_;i++) th[i] = 0.001;
  // for (int i=0;i<num_vec_;i++)  std::cout<<y[i]<<" "; std::cout<<"\n\n";

  for (int i=num_vec_ - 1; i>=0; i--){
    th[row[i]] = y[row[i]];
    for (int j=num_vec_ - 1; j>i; j--) th[row[i]] -= h_[row[i]][j]*th[row[j]];
    if (h_[row[i]][i] > eps)
      th[row[i]] /= h_[row[i]][i];
    else
      th[row[i]] = 0.;
  }
    //}
  // else{
  //   exit(0);
  //}

  std::cout<<"theta\n";
  //  for (int i=0;i<num_vec_;i++) th[i] = 0.001;
  for (int i=0;i<num_vec_;i++)  std::cout<<th[i]<<" "; std::cout<<"\n\n";
  //if (singular) exit(0);
  /// TEST
  double tt, tt_new;
  if (num_vec_ > 0){
    double alp = 1.;
    for (int i=0; i<num_vec_; i++) alp -= th[i];

    *fun = *dF_[new_f_];
    fun->Scale(alp);

    for (int i=0; i<num_vec_; i++){
      int I = (first_f_ + i)%(mvec_ + 1);
      std::cout<<"I="<<I<<"\n";
      fun->Update(th[i], *F_test[I], 1.);
    }
    fun->Norm2(&tt);
    std::cout<<"fun orig"<<tt<<"\n";
    std::srand(time(NULL));   
    for (int k=0;k<20;k++){

      for (int i=0; i<num_vec_; i++) th_test[i] = th[i] - 1 + 2*(double)std::rand()/(RAND_MAX);
      alp=1.;
      for (int i=0; i<num_vec_; i++) alp -= th_test[i];
      tt_new = 0.;

      *fun = *dF_[new_f_];
      fun->Scale(alp);
      
      for (int i=0; i<num_vec_; i++){
        int I = (first_f_ + i)%(mvec_ + 1);
        fun->Update(th_test[i], *F_test[I], 1.);
      }
      fun->Norm2( &tt_new);
      std::cout<<"fun new "<<tt_new<<"\n"; 
      if (tt_new - tt <0)    exit(0);
    }    
  }

  // std::cout<<"theta\n";                               
  // for (int i=0;i<num_vec_;i++) std::cout<<th[i]<<" "; std::cout<<"\n\n";

  alp = 0.;
  for (int i=0;i<num_vec_;i++) alp += th[i];
  alp = 1 - alp;

  // alp = 1.;// - 0.001*num_vec_;
  // for (int i=0;i<num_vec_;i++) th[i] = 0.;
  std::cout<<"alpha "<<alp<<"\n";
  

  // dir =  *u_[new_f_];
  // dir.Update(-alp, *u_[new_f_], 1.);
  // dir.Update(alp, *F_test[new_f_], 1.); 

  // for (int i=0; i<num_vec_; i++) {
  //   int I = (first_f_ + i)%(mvec_ + 1);
  //   std::cout<<"dir I "<<I<<"\n";
  //   dir.Update(-th[i], *u_[I], 1.);
  //   dir.Update(th[i], *F_test[I], 1.);
  // }

  beta_ = 0.95;

  dir =  *u_[new_f_];
  dir.Update(-alp, *u_[new_f_], 1.);
  dir.Update(beta_, *dF_[new_f_], 1.);

  for (int i=0; i<num_vec_; i++) {
    int I = (first_f_ + i)%(mvec_ + 1);
    dir.Update(-th[i]*beta_, *dF_[I], 1.);
    dir.Update(-th[i], *u_[I], 1.);
  }

  // dir.Scale(beta_);

  // if (num_vec_ > 0){
  //   for (int i=0; i<num_vec_; i++) {
  //     int I = (first_f_ + i)%(mvec_ + 1);
  //     dir.Update(-th[i], *u_[I], 1.);
  //   }
  //   dir.Update(alp, *u_[new_f_], 1.);
  // }


  //  std::cout<<"dir\n"; dir.Print(std::cout);

  delete[] y;
  delete[] b;
  delete[] th;



  if (num_vec_ < mvec_){
    first_f_ = 0;
    last_f_ = new_f_;
    if (new_f_++ > mvec_) new_f_ = 0;
    num_vec_ ++;
  }
  else {
    if (first_f_ >= mvec_) first_f_ = 0;
    else first_f_++;

    if (last_f_ >= mvec_) last_f_ = 0;
    else last_f_++;

    if (new_f_ >= mvec_) new_f_ = 0;
    else new_f_++;
  }




  std::cout<<"new_f "<<new_f_<<" first "<<first_f_<<" last "<<last_f_<<"\n";
  std::cout<<"\n";

};


}  // namespace AmanziSolvers
}  // namespace Amanzi

#endif
