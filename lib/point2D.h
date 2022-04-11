////////////////////////////////////////////////////////////////////////
//
// point2D.h
//
// contains everything that the polyclip algorithm needs for
// handling 2D points
//
////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cassert>

using namespace std;

////////////////////////////////////////////////////////////////////////
//
// a simple 2D point class with essential functionality
//
////////////////////////////////////////////////////////////////////////

class point2D
{
 public:
  // x and y coordinates of the point
  double x,y;

  // default-constructor to generate the 2D point (0,0)
  point2D() { x=y=0; }
   // copy-constructor
  point2D(const point2D& b) { x=b.x; y=b.y; }
  // initialize 2D point with explicit coordinates
  point2D(double X, double Y) { x=X; y=Y; }
  
  // destructor
  ~point2D() {}

  // assigns the values of a 2D point to another
  inline point2D operator=(const point2D& b) { x=b.x; y=b.y; return *this; }

  // sum of two 2D points
  inline point2D operator+(const point2D& b) const { return point2D(x+b.x, y+b.y); }
  // difference of two 2D points
  inline point2D operator-(const point2D& b) const { return point2D(x-b.x, y-b.y); }
  // adds 2D point b to the operand
  inline void operator+=(const point2D& b) { x+=b.x; y+=b.y; }
  // subtracts 2D point b from the operand
  inline void operator-=(const point2D& b) { x-=b.x; y-=b.y; }

  // multiply a 2D point by a scalar
  inline point2D operator*(double r) const { return point2D(r*x,r*y); }
  // multiply the operand by a scalar
  inline void operator*=(double r) { x*=r; y*=r; }
  // divide a 2D point by a scalar
  inline point2D operator/(double r) const { assert(r!=0.0); return point2D(x/r,y/r); }
  // divide the operand by a scalar
  inline void operator/=(double r) { assert(r!=0.0); x/=r; y/=r; }

  // calculate the dot-product
  inline double operator*(const point2D& b) const { return (x*b.x+y*b.y); }
  // calculate the cross-product
  inline double operator%(const point2D& b) const { return (x*b.y-y*b.x); }
};

// multiplication operator that allows the scalar value to preceed the 2D point
inline point2D operator*(double r, const point2D& v) { return v*r; }

// write a 2D point to a stream
inline ostream& operator<<(ostream& s, const point2D& v) { return (s << v.x << " " << v.y); }
// read a 2D point from a stream
inline istream& operator>>(istream& s, point2D& v) { return (s >> v.x >> v.y); }

// compute twice the signed area of the triange [P,Q,R]
inline double A(const point2D& P, const point2D& Q, const point2D& R) {
	return (Q.x-P.x) * (R.y-P.y) - (Q.y-P.y) * (R.x-P.x);
}

