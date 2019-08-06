// Mat & Vec definitions                 1.2.8 (unor 1997)


#ifndef _MatVec_h_

#define _MatVec_h_

/*
#ifdef MatVec
   #include <assert.h>
#else
   #define assert(expr) // expr 
#endif
*/

#include <iostream>
//using namespace std;
using std::cout;
using std::endl;
using std::cin;
using std::cerr;

class Mat;      // matrix class

class Vec;      // vector class

class SVD;      // singular value decomposition class


typedef double Float;
typedef double Double;

// ----------------------------------------------------------------------------
#define MAT_NONE		0
#define MAT_ZERO		1
#define MAT_IDENTITY		2


class Mat {
public:
   friend class Vec;
   friend class SVD;
   enum init {none, zero, identity};
   Mat() { Mat_malloc(0, 0); };
   Mat(int r, int c, int mi=MAT_NONE);//Mat::init mi=Mat::none);
   Mat(const Mat&);
   virtual ~Mat();
   Mat& operator=(const Mat&);

   Float* operator[](int r) { return row_ptr[r]; }
   const Float* const operator[](int r) const { return row_ptr[r]; }

   int rows() const { return row; }
   int cols() const { return col; }

   Mat& reset(int r=0, int c=0, int mi=MAT_NONE);

   Mat& add(const Mat&);
   Mat& sub(const Mat&);
   Mat& operator+=(const Mat& A) { return add(A); }
   Mat& operator-=(const Mat& A) { return sub(A); }

   Mat& mult(Float);
   Mat& operator*=(Float s) { return mult(s); }

   Mat& trans(const Mat&);

   Mat& mult_AB(const Mat&, const Mat&);
   Mat& mult_trA_B(const Mat&, const Mat&);
   Mat& mult_trA_A(const Mat& A) { return mult_trA_B(A, A); }
   Mat& mult_A_trB(const Mat&, const Mat&);
   Mat& mult_A_d_B(const Mat&, const Vec&, const Mat&);
   Mat& mult_trA_d_B(const Mat&, const Vec&, const Mat&);
   Mat& mult_A_d_trB(const Mat&, const Vec&, const Mat&);

   Mat& set_diag(const Vec&);
   Mat& set_row(int, const Vec&);
   Mat& set_col(int, const Vec&);

   // a submatrix is just a reference to the part of an existing matrix


   Mat& submatrix(const Mat&, int rl, int ru, int cl, int cu);

   // basic IO methods


   friend std::istream& operator>>(std::istream&, Mat&);
   friend std::ostream& operator<<(std::ostream&, const Mat&);

   // methods with temporary Mat


   Mat operator+(const Mat& A) const { Mat T = *this; return T += A; }
   Mat operator-(const Mat& A) const { Mat T = *this; return T -= A; }
   Mat operator*(Float s) const { Mat T = *this; return T *= s; }
   friend Mat operator*(Float s, const Mat& A) { Mat T = A; return T *= s; }
   Mat operator*(const Mat& A) const { Mat T; return T.mult_AB(*this, A); }
   friend Mat trans(const Mat& A) { Mat T; return T.trans(A); }
   friend Mat inv(const Mat& A);   // inv(A) * A == identity


private:
   int row, col;
   Float* mem;
   Float** row_ptr;
   void Mat_malloc(int r, int c);
   void Mat_free();
   void Mat_copy(const Mat&);
   void Mat_init(int);

};   // Mat


// ----------------------------------------------------------------------------


class Vec {
public:
   friend class Mat;
   friend class SVD;
   Vec() : dim_v(0), mem_v(0) {}
   Vec(int r) : dim_v(r), mem_v(new Float[r] - 1) {}
   Vec(const Vec&);
   virtual ~Vec();
   Vec& operator=(const Vec&);

   Float& operator[](int n) { return  mem_v[n]; }
   const Float& operator[](int n) const { return mem_v[n]; }

   Vec& reset(int=0);

   int dim() const { return dim_v; }

   Float L2_norm() const;
   Float L1_norm() const;
   Float Lmax_norm() const;

   Vec& add(const Vec&);
   Vec& sub(const Vec&);
   Vec& operator+=(const Vec& A) { return add(A); }
   Vec& operator-=(const Vec& A) { return sub(A); }

   Vec& mult(Float);
   Vec& operator*=(Float s) { return mult(s); }

   Vec& add(const Vec&, const Vec&);
   Vec& sub(const Vec&, const Vec&);
   Vec& add(const Vec&, const Vec&, Float);
   Vec& sub(const Vec&, const Vec&, Float);

   Vec& mult_Ab(const Mat& A, const Vec& b);
   Vec& mult_trA_b(const Mat& A, const Vec& b);
   Vec& solve(SVD& svd, const Vec& b);

   friend std::istream& operator>>(std::istream&, Vec&);
   friend std::ostream& operator<<(std::ostream&, const Vec&);

   // methods with temporary Vec


   Vec operator+(const Vec& v) const { Vec tmp = *this; tmp += v; return tmp; }
   Vec operator-(const Vec& v) const { Vec tmp = *this; tmp -= v; return tmp; }
   Vec operator*(Float s) const { Vec tmp = *this; tmp *= s; return tmp; }
   friend Vec operator*(Float s, const Vec& v) { Vec t = v; t *= s; return t; }
   friend Vec operator*(const Mat& A, const Vec& b)
      { Vec T; return T.mult_Ab(A, b); }

private:
   int dim_v;
   Float* mem_v;

};   // Vec



// ----------------------------------------------------------------------------


/*

 Class SVD is a simple wrapper for the Singular Value Decomposition algorithm
 and is based on the SVD source in C/C++ from the MatClass by C. R. Birchenhall

 Ales Cepek <cepek@fsv.cvut.cz>

*/

#include "matvec.h"

#include <stdlib.h>


class SVD {

public:
   friend class Mat;
   friend class Vec;
   SVD() : U(), m(0), n(0), W(), V(), decomposed(0), W_tol(0), inv_W() {}
//   SVD(const Mat& A) : U(A), m(A.rows()), n(A.cols(),MAT_NONE), W(n), V(n, n),
//                       decomposed(0), W_tol(0), inv_W(n) {}
   SVD(const Mat& A);

   Float tol() { return W_tol; }
   Float tol(Float t) { W_tol = t; set_inv_W(); return W_tol; }

   SVD& decompose() { svd(); return *this; }
   SVD& reset(const Mat& A);
   SVD& clear();

   Float cov_xx(int, int);   // covariance (xi, xj) of adjusted unknowns

   Float cov_bb(int, int);   // covariance (bi, bj) of adjusted observations

   Float cov_bx(int, int);   // covariance (bi, xj)


   const Mat& SVD_U() { svd(); return U; }
   const Vec& SVD_W() { svd(); return W; }
   const Mat& SVD_V() { svd(); return V; }
   const Vec& SVD_invW() {svd();return inv_W;}

protected:
   SVD& reset(const Mat& A, const Vec& w);
   virtual void error(char* s) {cerr << s << endl; exit(1); }

private:
   Mat U;
   int m, n;
   Vec W;
   Mat V;

   int decomposed;
   void svd();

   Float W_tol;
   Vec inv_W;
   void set_inv_W();
};

#endif
