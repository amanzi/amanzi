/*
  Solvers

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Daniil Svyatskiy
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

  void QRdelete();
  void TestQR(int nv);
  void Relax();
  void Restart();
  void Correction(const Vector&, Vector&, const Teuchos::Ptr<const Vector>& oldv=Teuchos::null);

 private:
  int subspace_;  // boolean: a nonempty subspace
  int pending_;   // contains pending vectors -- boolean
  int mvec_;      // maximum number of subspace vectors
  double vtol_;   // vector drop tolerance
  double beta_;   // relaxation parameter

  Teuchos::RCP<Vector> *u_;  // previous iterates
  Teuchos::RCP<Vector> *dF_;  // preconditioner residuals
  Teuchos::RCP<Vector> *dG_;  // preconditioner residuals
  Teuchos::RCP<Vector> *Q_;  // Q matrix in QR decomposition
  Teuchos::RCP<Vector> *F_test;  // preconditioner residuals

  //  double **h_;  // matrix of w vector inner products 
  double *R_;


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
  dG_ = new Teuchos::RCP<Vector> [mvec_ + 1];
  F_test = new Teuchos::RCP<Vector> [mvec_ + 1];
  Q_ = new Teuchos::RCP<Vector> [mvec_ + 1];

  for (int i = 0; i < mvec_ + 1; i++) {
    u_[i] = Teuchos::rcp(new Vector(map));
    dF_[i] = Teuchos::rcp(new Vector(map));
    dG_[i] = Teuchos::rcp(new Vector(map));
    F_test[i] = Teuchos::rcp(new Vector(map));
    Q_[i] = Teuchos::rcp(new Vector(map));
  }

  // h_ = new double* [mvec_];
  // for (int j = 0; j < mvec_; j++) {
  //    h_[j] = new double[mvec_];
  // }

  R_ = new double[mvec_*(mvec_+1)/2];

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
  // for (int j = 0; j <mvec_; ++j) {
  //   delete [] h_[j];
  // }
  // delete [] h_;
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


template<class Vector, class VectorSpace>
void AA_Base<Vector, VectorSpace>::QRdelete() {

  //std::cout<<"*************** QR DELETE ******************\n";

  Teuchos::RCP<Vector> temp_vec = Teuchos::rcp(new Vector(*Q_[0]));
  int m = num_vec_ - 1;

  // for (int i=0;i< m;i++) {
  //   //int I = (first_f_ + i)%(mvec_ + 1);
  //   std::cout<<"Q "<<i<<"\n";
  //   Q_[i]->Print(std::cout);
  // }

  // std::cout<<"R matrix "<<"\n";
  // for (int i=0;i< m; i++) {
  //   for (int j=i;j< m; j++) {
  //     int loc_id=j*(j+1)/2 + i;
  //     std::cout<<R_[loc_id]<<" ";
  //   }
  //   std::cout<<"\n";
  // }


  for (int i=0; i < m - 1; i++) {
    int loc_id = (i+2)*(i+3)/2 - 1; 

    double temp =  sqrt( R_[loc_id]*R_[loc_id] + R_[loc_id-1]*R_[loc_id-1]);
    double c = R_[loc_id-1]/temp;
    double s =   R_[loc_id]/temp;

    R_[loc_id-1] = temp;
    R_[loc_id] = 0.;
    if (i < m - 2) {
      for (int j=i+2; j<m; j++) {
        loc_id = j*(j+1)/2 + i;
        temp = c*R_[loc_id] + s*R_[loc_id + 1];
        R_[loc_id + 1] = -s*R_[loc_id] + c*R_[loc_id + 1];
        R_[loc_id] = temp;          
      }
    }

    *temp_vec = *Q_[i];
    temp_vec->Update(s, *Q_[i+1], c);
    Q_[i+1]->Update(-s, *Q_[i], c);
   *Q_[i] = *temp_vec;
  }


  for (int i=0;i<m-1; i++) {
    for (int j=0;j<=i ; j++) {
      int new_id =i*(i+1)/2 + j;
      int old_id =(i+1)*(i+2)/2 + j;
      R_[new_id] = R_[old_id];
    }
  }

  first_f_=(first_f_ + 1)%(mvec_ + 1);
  num_vec_--;


  // for (int i=0;i<m-1;i++) {
  //   //int I = (first_f_ + i)%(mvec_ + 1);
  //   std::cout<<"Q "<<i<<"\n";
  //   Q_[i]->Print(std::cout);
  // }

  // std::cout<<"R matrix "<<"\n";
  // for (int i=0;i<m-1; i++) {
  //   for (int j=i;j<m-1; j++) {
  //     int loc_id=j*(j+1)/2 + i;
  //     std::cout<<R_[loc_id]<<" ";
  //   }
  //   std::cout<<"\n";
  // }
  //exit(0);

}


template<class Vector, class VectorSpace>
void AA_Base<Vector, VectorSpace>::TestQR(int nv) {

  std::cout<<"*************** Test QR ******************\n";
  Teuchos::RCP<Vector> ff = Teuchos::rcp(new Vector(*Q_[0]));

  double norm2 = 0;
  int loc_id = 0.;
  ff->PutScalar(0.);
  std::cout << "num_vec_ "<<num_vec_<<"\n";
  for (int i=0; i<nv; i++) {
    ff->PutScalar(0.);
    for (int j=0; j<=i; j++) {
      ff->Update(R_[loc_id], *Q_[j], 1.);
      loc_id++;
      // Q_[j]->Print(std::cout);
    }
    int I = (first_f_ + i)%(mvec_ + 1);
    ff->Update(-1., *dF_[I], 1.);

    //dF_[I]->Print(std::cout);

    ff->Norm2(&norm2);
    std::cout<<"norm2 "<<norm2<<"\n";
  }

  if (norm2 > 0.0001) exit(0);
}


/* ******************************************************************
 * TBW
 ***************************************************************** */
template<class Vector, class VectorSpace>
void AA_Base<Vector, VectorSpace>::Correction(const Vector& f, Vector &dir,
        const Teuchos::Ptr<const Vector>& u_old)
{
  Teuchos::RCP<Vector> vp, wp;
  Teuchos::RCP<Vector> ff = Teuchos::rcp(new Vector(f));
  Teuchos::RCP<Vector> fun = Teuchos::rcp(new Vector(f));
  Teuchos::RCP<Vector> tmp = Teuchos::rcp(new Vector(f));

  // double norm_u, norm_f;
  // ff->Norm2(&norm_f);
  // u_old->Norm2(&norm_u);

  assert(new_f_ >= 0);
  // Save the accelerated correction for the next_f_ call.

  //ff->Scale(-1.);

  *dF_[new_f_] = *ff;                             // dF_new = ff
  //*F_test[new_f_] = *ff;
  if (u_old != Teuchos::null) {
    *u_[new_f_] = *u_old;                         // u_new = u_old

    *dG_[new_f_] = *u_old;
    dG_[new_f_]->Update(-1., *ff, 1);            // dG_new = u - ff
  }

  if (last_f_ >= 0) {
    dF_[last_f_]->Update(1., *dF_[new_f_], -1.);   //dF_last =  dF_new - dF_last
    dG_[last_f_]->Update(1., *dG_[new_f_], -1.);   //dG_last = dG_new  - dG_last
  }

  if (num_vec_ == 1) {
    if (last_f_ == 0) {
      double norm2;
      *Q_[last_f_] = *dF_[last_f_];
      Q_[last_f_]->Norm2(&norm2);
      Q_[last_f_]->Scale(1./norm2);
      R_[0] = norm2;
    }
  } else if (num_vec_ > 1) {
    if (num_vec_ == mvec_) {
    // Delete old Vector
      QRdelete();
      //TestQR(num_vec_ - 1);
    }

    double norm2 = 0.;
    int last_col_i;
    while ((norm2 < 1e-12)&&(num_vec_ > 1)) {
      last_col_i = num_vec_*(num_vec_ - 1)/2;    
      *tmp = *dF_[last_f_];
      for (int i=0; i<num_vec_ - 1; i++) {
        double val;
        tmp->Dot(*Q_[i], &val);
        R_[last_col_i + i] = val;
        tmp->Update(-val, *Q_[i], 1.);
      }    
      tmp->Norm2(&norm2); 

      if (norm2 < 1e-12) {
        QRdelete();       
        //TestQR(num_vec_ - 1);
        //exit(0);
      }
    }

    *Q_[num_vec_ - 1] = *tmp;
    Q_[num_vec_ - 1]->Scale(1./norm2);
    R_[last_col_i + num_vec_ - 1] = norm2;     

  }

  // for (int i=0;i<num_vec_;i++) {
  //   //    int I = (first_f_ + i)%(mvec_ + 1);
  //   std::cout<<"Q "<<i<<"\n";
  //   Q_[i]->Print(std::cout);
  // }

  // std::cout<<"R matrix "<<"\n";
  // for (int i=0;i<num_vec_; i++) {
  //   for (int j=i;j<num_vec_; j++) {
  //     int loc_id=j*(j+1)/2 + i;
  //     std::cout<<R_[loc_id]<<" ";
  //   }
  //   std::cout<<"\n";
  // }

  //TestQR(num_vec_);

      
  double *b  = new double[num_vec_];
  double *th = new double[num_vec_];

  // Solve system R*th = Q^t * F_new

  for (int i=0; i<num_vec_; i++) {
    th[i] = 0.;
    dF_[new_f_]->Dot(*Q_[i], &b[i]);
  }

  for (int i = num_vec_ - 1; i>=0; i--) {
    th[i] = b[i];
    int diag_id = (i+1)*(i+2)/2 - 1;
    for (int j=i+1; j<num_vec_; j++) {
      int loc_id = j*(j+1)/2 + i;
      th[i] -= R_[loc_id]*th[j];
    }
    th[i] = th[i]/R_[diag_id];
  }
  
  // double alp;
  // alp = 1.;// - 0.001*num_vec_;

  // std::cout<<"theta\n";
  // for (int i=0;i<num_vec_;i++) std::cout<<th[i]<<" ";
  // std::cout<<"\n";
  // std::cout<<"alpha "<<alp<<"\n";
  

  // dir =  *u_[new_f_];
  // dir.Update(-alp, *u_[new_f_], 1.);
  // dir.Update(alp, *F_test[new_f_], 1.); 

  //  dir.PutScalar(0.);

  dir = *dF_[new_f_];
  for (int i=0; i<num_vec_; i++) {
    int I = (first_f_ + i)%(mvec_ + 1);
    dir.Update(th[i], *dG_[I], 1.);
  }
  // std::cout<<"dir\n";
  // dir.Print(std::cout);

  if ((beta_ >0)&&(beta_ != 1.)) {
    dir.Update(beta_ - 1, *dF_[new_f_], 1);
    for (int i=0; i<num_vec_; i++) {
      int I = (first_f_ + i)%(mvec_ + 1);
      dir.Update((1-beta_)*th[i], *dF_[I], 1.);
    }
  }


  delete[] b;
  delete[] th;


  if (num_vec_ < mvec_) {
    if (first_f_ < 0) {
      first_f_ = 0;
    }
    last_f_ = new_f_;
    new_f_++;
    if (new_f_ > mvec_) new_f_ = 0;
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

  // std::cout<<"new_f "<<new_f_<<" first "<<first_f_<<" last "<<last_f_<<"\n";
  // std::cout<<"\n";
};

}  // namespace AmanziSolvers
}  // namespace Amanzi

#endif
