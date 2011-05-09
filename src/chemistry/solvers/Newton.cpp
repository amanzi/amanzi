#include "Newton.hpp"

Newton::Newton(const int n) {
  
  size(n);
  x_.resize(n);
  r_.resize(n);
  J_ = new Block(n);
  indices_.resize(n);
  vv_.resize(n);

} // end Newton constructor

void Newton::solve() {
  cout << "Solved!\n";
} // end solve()


void Newton::LUDecomposition()
{
  const double small_number = 1.e-20;
  int i,imax,j,k;
  double big,dum,sum,temp;

	d_=1.0;
	for (i=0;i<size();i++) {
		big=0.0;
		for (j=0;j<size();j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
    if (big == 0.0) std::cout << "Singular matrix in routine ludcmp";
		vv[i]=1.0/big;
	}
	for (j=0;j<size();j++) {
		for (i=0;i<j;i++) {
			sum=a[i][j];
			for (k=0;k<i;k++) sum -= a[i][k]*a[k][j];
			(J_->getValue())[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<size();i++) {
			sum=a[i][j];
			for (k=0;k<j;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ((dum=vv_[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=0;k<size();k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			d_ = -d_;
			vv_[imax]=vv[j];
		}
		indx_[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=small_number;
		if (j != size()-1) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<size();i++) a[i][j] *= dum;
		}
  }
  delete [] vv;
} // end ludcmp()


#undef TINY

void LUBackSolve(double **a, int n, int *indx, std::vector<double> &b)
{
  int i,ii=0,ip,j;
  double sum;

	for (i=0;i<n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii != 0)
			for (j=ii-1;j<i;j++) sum -= a[i][j]*b[j];
		else if (sum != 0.0)
			ii=i+1;
		b[i]=sum;
	}
	for (i=n-1;i>=0;i--) {
		sum=b[i];
		for (j=i+1;j<n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
} // end lubksb()

Newton::~Newton() {
  if (J) delete J;
  J = NULL;
} // end Newton destructor
