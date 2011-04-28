#ifndef   __POINT_HXX
#define   __POINT_HXX

#include  <iostream>
#include  <math.h>


namespace Amanzi
{

  namespace AmanziGeometry
  {


    class Point 
    {
    private:
      int    d;
      double *xyz;

    public:
      Point()                                                    { d = 0; xyz = NULL; }
      Point( const int N )                                       { init(N); }
      Point( const double& x, const double& y )                  { init(2); xyz[0] = x; xyz[1] = y;  }
      Point( const double& x, const double& y, const double& z ) { init(3); xyz[0] = x; xyz[1] = y; xyz[2] = z; }

      ~Point() { if ( xyz != NULL ) delete [] xyz; }

      /* initialization */
      void init( const int N ) { d = N;  xyz = new double[d]; for( int i=0; i<d; i++) xyz[i] = 0.0; }

      void set( const double& val ) { for( int i=0; i<d; i++ ) xyz[i] = val;  }
      void set( const Point& p )    { for( int i=0; i<d; i++ ) xyz[i] = p[i]; }
      void set( const double& x, const double& y ) { xyz[0] = x; xyz[1] = y; }
      void set( const double& x, const double& y, const double& z ) { xyz[0] = x; xyz[1] = y; xyz[2] = z; }

      int is_valid() { return (d==2 || d==3) ? 1 : 0; } 

      /* access */
      double& operator[] (const int i)       { return xyz[i]; }
      const double& operator[] (const int i) const { return xyz[i]; }

      const double x() const { return xyz[0]; }
      const double y() const { return xyz[1]; }
      const double z() const { return (this->d==3) ? xyz[2] : NULL; }

      const int dim() const { return d; }

      /* operators */
      Point& operator=( const double& val ) { for( int i=0; i<d; i++ ) xyz[i] = val;  return *this; }
      Point& operator=( const Point& p)     { for( int i=0; i<d; i++ ) xyz[i] = p[i]; return *this; }         

      Point& operator+=( const Point& p )  { for( int i=0; i<d; i++ ) xyz[i] += p[i]; return *this; }
      Point& operator-=( const Point& p )  { for( int i=0; i<d; i++ ) xyz[i] -= p[i]; return *this; }
      Point& operator*=( const double& d ) { for( int i=0; i<d; i++ ) xyz[i] *= d;    return *this; }
      Point& operator/=( const double& d ) { for( int i=0; i<d; i++ ) xyz[i] /= d;    return *this; }

      friend Point  operator*( const double& r, const Point& p )  { return (p.d==2) ? Point( r*p[0], r*p[1] ) : Point( r*p[0], r*p[1], r*p[2] ); }
      friend Point  operator*( const Point& p,  const double& r ) { return r*p; }
      friend double operator*( const Point& p,  const Point& q )  { double s = 0.; for( int i=0; i<p.d; i++ ) s += p[i]*q[i]; return s; }

      friend Point  operator/( const Point& p,  const double& r ) { return p * (1./r); }

      friend Point  operator+( const Point& p, const Point& q ) { return (p.d==2) ? Point( p[0]+q[0], p[1]+q[1] ) : Point( p[0]+q[0], p[1]+q[1], p[2]+q[2] ); }
      friend Point  operator-( const Point& p, const Point& q ) { return (p.d==2) ? Point( p[0]-q[0], p[1]-q[1] ) : Point( p[0]-q[0], p[1]-q[1], p[2]-q[2] ); } 
      friend Point  operator-( const Point& p )                 { return (p.d==2) ? Point( -p[0], -p[1] ) : Point( -p[0], -p[1], -p[2] ); }

      friend Point  operator^( const Point& p, const Point& q ) {
	Point pq(p.d);
	if( p.d==2 ) {
	  pq[0] = p[0] * q[1] - q[0] * p[1];
	} 
	else if( p.d==3 ) {
	  pq[0] = p[1] * q[2] - p[2] * q[1];
	  pq[1] = p[2] * q[0] - p[0] * q[2];
	  pq[2] = p[0] * q[1] - p[1] * q[0];
	}
	return pq;
      }

      /* miscaleneous */
      friend double L22(  const Point& p ) { return p*p; }
      friend double norm( const Point& p ) { return sqrt( p*p ); }

      friend ostream& operator<<( ostream& os, Point& p) {
	os << p.x() << " " << p.y();
	if( p.d == 3 ) os << " " << p.z();
	return os;
      }

    }; // end class Point

  } // end namespace AmanziGeometry

} // end namespace Amanzi

#endif

