#include <iomanip>
#include <math.h>
#include <assert.h>
#include "matvec.h"


// ---------  Mat methods  ---------------------------------------------


Mat::Mat(int r, int c, int mi) {
   Mat_malloc(r, c);
   if (mi == MAT_NONE)
      return;
   Mat_init(mi);
}

Mat::Mat(const Mat& A) {
   Mat_malloc(A.row, A.col);
   Mat_copy(A);
}

Mat::~Mat() {
   Mat_free();
}

Mat& Mat::operator=(const Mat& A) {
   if (row != A.row || col != A.col) {
      Mat_free();
      Mat_malloc(A.row, A.col);
   }
   Mat_copy(A);
   return *this;
}

Mat& Mat::reset(int r, int c, int m) {
   assert(r >= 0 && c >= 0);
   Mat_free();
   Mat_malloc(r, c);
   Mat_init(m);
   return *this;   
}

Mat& Mat::add(const Mat& A) {
   assert(row == A.row && col == A.col);
   for (int i = 1; i <= A.row; i++)
      for (int j = 1; j <= A.col; j++)
	 (*this)[i][j] += A[i][j];
   return *this;
}

Mat& Mat::sub(const Mat& A) {
   assert(row == A.row && col == A.col);
   for (int i = 1; i <= A.row; i++)
      for (int j = 1; j <= A.col; j++)
	 (*this)[i][j] -= A[i][j];
   return *this;
}

Mat& Mat::mult(Float s) {
   for (int i = 1; i <= row; i++)
      for (int j = 1; j <= col; j++)
	 (*this)[i][j] *= s;
   return *this;
}

Mat& Mat::mult_AB(const Mat& A, const Mat& B) {
   assert(A.col == B.row);
   if (row != A.row || col != B.col) {
      Mat_free();
      Mat_malloc(A.row, B.col);
   }
   for (int k = 1; k <= row; k++)
      for (int l = 1; l <= col; l++)
	 (*this)[k][l] = 0;
   Float sum;
   for (int i = 1; i <= row; i++)
      for (int j = 1; j <= col; j++) {
	 sum = 0;
	 for (int s = 1; s <= A.col; s++)
	    sum += A[i][s] * B[s][j];
	 (*this)[i][j] = sum;
      }
   return *this;
}

Mat& Mat::mult_trA_B(const Mat& A, const Mat& B) {
   assert(A.row == B.row);
   if (row != A.col || col != B.col) {
      Mat_free();
      Mat_malloc(A.col, B.col);
   }
   for (int k = 1; k <= row; k++)
      for (int l = 1; l <= col; l++)
	 (*this)[k][l] = 0;
   Float sum;
   for (int i = 1; i <= row; i++)
      for (int j = 1; j <= col; j++) {
	 sum = 0;
	 for (int s = 1; s <= A.row; s++)
	    sum += A[s][i] * B[s][j];
	 (*this)[i][j] = sum;
      }
   return *this;
}

Mat& Mat::mult_A_trB(const Mat& A, const Mat& B) {
   assert(A.col == B.col);
   if (row != A.row || col != B.row) {
      Mat_free();
      Mat_malloc(A.row, B.row);
   }
   for (int k = 1; k <= row; k++)
      for (int l = 1; l <= col; l++)
	 (*this)[k][l] = 0;
   Float sum;
   for (int i = 1; i <= row; i++)
      for (int j = 1; j <= col; j++) {
	 sum = 0;
	 for (int s = 1; s <= A.col; s++)
	    sum += A[i][s] * B[j][s];
	 (*this)[i][j] = sum;
      }
   return *this;
}

Mat& Mat::mult_A_d_B(const Mat& A, const Vec& d, const Mat& B) {
   assert(A.col == d.dim_v && d.dim_v == B.row);
   if (row != A.row || col != B.col) {
      Mat_free();
      Mat_malloc(A.row, B.col);
   }
   for (int k = 1; k <= row; k++)
      for (int l = 1; l <= col; l++)
	 (*this)[k][l] = 0;
   Float sum;
   for (int i = 1; i <= row; i++)
      for (int j = 1; j <= col; j++) {
	 sum = 0;
	 for (int s = 1; s <= A.col; s++)
	    sum += A[i][s] * d[s] * B[s][j];
	 (*this)[i][j] = sum;
      }
   return *this;
}

Mat& Mat::mult_trA_d_B(const Mat& A, const Vec& d, const Mat& B) {
   assert(A.row == d.dim_v && d.dim_v == B.row);
   if (row != A.col || col != B.col) {
      Mat_free();
      Mat_malloc(A.col, B.col);
   }
   for (int k = 1; k <= row; k++)
      for (int l = 1; l <= col; l++)
	 (*this)[k][l] = 0;
   Float sum;
   for (int i = 1; i <= row; i++)
      for (int j = 1; j <= col; j++) {
	 sum = 0;
	 for (int s = 1; s <= A.row; s++)
	    sum += A[s][i] * d[s] * B[s][j];
	 (*this)[i][j] = sum;
      }
   return *this;
}

Mat& Mat::mult_A_d_trB(const Mat& A, const Vec& d, const Mat& B) {
   assert(A.col == d.dim_v && d.dim_v == B.col);
   if (row != A.row || col != B.row) {
      Mat_free();
      Mat_malloc(A.row, B.row);
   }
   for (int k = 1; k <= row; k++)
      for (int l = 1; l <= col; l++)
	 (*this)[k][l] = 0;
   Float sum;
   for (int i = 1; i <= row; i++)
      for (int j = 1; j <= col; j++) {
	 sum = 0;
	 for (int s = 1; s <= A.col; s++)
	    sum += A[i][s] * d[s] * B[j][s];
	 (*this)[i][j] = sum;
      }
   return *this;
}

Mat& Mat::trans(const Mat& A) {
   if (row != A.col || col != A.row) {
      Mat_free();
      Mat_malloc(A.col, A.row);
   }
   for (int i = 1; i <= A.row; i++)
      for (int j = 1; j <= A.col; j++)
	 (*this)[j][i] = A[i][j];
   return *this;
}

Mat& Mat::set_diag(const Vec& a) {
   assert(row >= a.dim() && col >= a.dim() && (row == a.dim() || col == a.dim()));
   for (int i=1; i <= a.dim(); i++)
      (*this)[i][i] = a[i];
   return *this;
}

Mat& Mat::set_row(int n, const Vec& a) {
   assert(n >= 0 && n <= row && col == a.dim());
   for (int i=1; i <= a.dim(); i++)
      (*this)[n][i] = a[i];
   return *this;
}

Mat& Mat::set_col(int n, const Vec& a) {
   assert(n >= 0 && n <= row && row == a.dim());
   for (int i=1; i <= a.dim(); i++)
      (*this)[i][n] = a[i];
   return *this;
}

// ---------  Submatrices  -----------------------------------------------


Mat& Mat::submatrix(const Mat& M, int rl, int ru, int cl, int cu) {
   const int n__rows = ru - rl + 1;
   const int n__cols = cu - cl + 1;
   assert(n__rows >= 0 && n__cols >= 0);
   Mat_free();      // always free memory when creating a submatrix

   Mat_malloc(n__rows, 0);
   col = n__cols;
   for (int k=1, i=rl; i<=ru; k++, i++)
      row_ptr[k] = (Float*)(&M[i][cl]) - 1;   // cast discarding const

   return *this;
}

// ---------  Mat friends  -----------------------------------------------


std::istream& operator>>(std::istream& istr, Mat& A) {
   int r = 0, c = 0;
   assert(istr.good());
   istr >> r >> c;
   if (r != A.row || c != A.col) {
      A.Mat_free();
      A.Mat_malloc(r, c);
   }

   for (int i = 1; i <= r; i++)
      for (int j = 1; j <= c; j++)
	 istr >> A[i][j];

   return istr;
}

std::ostream& operator<<(std::ostream& ostr, const Mat& A) {
   int aw = ostr.width();
   if (aw > 0) aw--;
   ostr << std::setw(0);
   ostr << ' ' << std::setw(aw) << A.row << ' ' << std::setw(aw) << A.col << "\n\n";
   for (int i = 1; i <= A.row; i++) {
      for (int j = 1; j <= A.col; j++) {
	 ostr << ' ' << std::setw(aw) << A[i][j];
      }
      ostr << endl;
   }
 
   return ostr;
}

Mat inv(const Mat& A) {
   assert(A.rows() == A.cols());
   const int N = A.rows();
   Mat T(N, N);
   Vec c(N), b(N);
   SVD svd(A);
   svd.decompose();
   for (int j = 1; j <= N; j++)
      b[j] = 0;
   for (int i = 1; i <= N; i++) {
      b[i] = 1;
      c.solve(svd, b);
      T.set_col(i, c);
      b[i] = 0;
   }
   return T;
}

// ---------  Mat private functions  --------------------------------------


void Mat::Mat_malloc(int r, int c) {
   if( r == 0) {
      row = col = 0;
      mem = 0;
      return;
   }
   row = r;
   col = c;
   const int row_ptr_dim = ((sizeof(Float*)+sizeof(Float)-1)/sizeof(Float))*r;
   const int mem_dim = row_ptr_dim + r*c; 
   mem = new Float[mem_dim];
   row_ptr = (Float**)mem;
   Float* pf = mem + row_ptr_dim - 1;
   for (int i = 0; i < r; i++, pf += c)
      row_ptr[i] = pf;
   row_ptr--;
}

void Mat::Mat_free() {
   if (mem == 0) return;
   delete[] mem;
}

void Mat::Mat_copy(const Mat& A) {
   for (int i = 1; i <= A.row; i++)
      for (int j = 1; j <= A.col; j++)
	 (*this)[i][j] = A[i][j];
}

void Mat::Mat_init(int mi) {
   if (mem == 0)
      return;
   // if (mi == MAT_ZERO || mi == MAT_IDENTITY)

   for (int i=1; i<=row; i++) 
      for (int j=1; j<=col; j++)
	 (*this)[i][j] = 0;
   if (mi == MAT_ZERO)
      return;
   const int k = (row <= col) ? row : col;
   for (int d=1; d<=k; d++)
      (*this)[d][d] = 1;
}

// ---------  Vec methods  --------------------------------------------


Vec::~Vec() { if (mem_v != 0) delete[] ++mem_v; }

Vec::Vec(const Vec& v) {
   dim_v = v.dim_v;
   mem_v = new Float[dim_v] - 1;
   Float* old_m = v.mem_v;
   Float* new_m = mem_v;
   for (int i=0; i<dim_v; i++)
      *++new_m = *++old_m;    
}

Vec& Vec::operator=(const Vec& v) {
   if (dim_v != v.dim_v) {
      if (mem_v != 0) 
	 delete[] ++mem_v;  
      dim_v = v.dim_v;
      mem_v = new Float[dim_v] - 1;
   }
   Float* old_m = v.mem_v;
   Float* new_m = mem_v;
   for (int i=0; i<dim_v; i++)
      *++new_m = *++old_m;    
   return *this;
}

Vec& Vec::reset(int d) {
   if (mem_v != 0)
      delete[] ++mem_v;
   if (d > 0)
      mem_v = new Float[d] - 1;
   else
      mem_v = 0;
   dim_v = d;
   return *this;
}

Float Vec::L2_norm() const {
   Float s = 0;
   Float* p = mem_v + 1;
   for (int i=1; i<=dim(); i++, p++)
      s += (*p) * (*p);
   return sqrt(s);
}

Float Vec::L1_norm() const {
   Float s = 0;
   Float* p = mem_v + 1;
   for (int i=1; i<=dim(); i++, p++)
      s += fabs(*p);
   return s;
}

Float Vec::Lmax_norm() const {
   Float s = 0, as;
   Float* p = mem_v + 1;
   for (int i=1; i<=dim(); i++, p++) {
      as = fabs(s);
      if (as > s)
	 s = as;
   }
   return s;
}

Vec& Vec::add(const Vec& b) 
{
   assert(dim_v == b.dim_v);
   Float* arg_m = b.mem_v;
   Float* new_m = mem_v;
   for (int i=0; i<dim_v; i++)
      *++new_m += *++arg_m;
   return *this;
}

Vec& Vec::sub(const Vec& b) 
{
   assert(dim_v == b.dim_v);
   Float* arg_m = b.mem_v;
   Float* new_m = mem_v;
   for (int i=0; i<dim_v; i++)
      *++new_m -= *++arg_m;
   return *this;
}

Vec& Vec::mult(Float s) 
{
   Float* new_m = mem_v;
   for (int i=0; i<dim_v; i++)
      *++new_m *= s;
   return *this;
}

Vec& Vec::add(const Vec& a, const Vec& b) {
   assert(dim() == a.dim() && a.dim() == b.dim());   
   for (int i=1; i<=dim_v; i++)
      (*this)[i] = a[i] + b[i];
   return *this;
}

Vec& Vec::sub(const Vec& a, const Vec& b) {
   assert(dim() == a.dim() && a.dim() == b.dim());   
   for (int i=1; i<=dim_v; i++)
      (*this)[i] = a[i] - b[i];
   return *this;
}

Vec& Vec::add(const Vec& a, const Vec& b, Float q) {
   assert(dim() == a.dim() && a.dim() == b.dim());   
   for (int i=1; i<=dim_v; i++)
      (*this)[i] = a[i] + b[i]*q;
   return *this;
}

Vec& Vec::sub(const Vec& a, const Vec& b, Float q) {
   assert(dim() == a.dim() && a.dim() == b.dim());   
   for (int i=1; i<=dim_v; i++)
      (*this)[i] = a[i] - b[i]*q;
   return *this;
}

Vec& Vec::mult_Ab(const Mat& A, const Vec& b) {
   assert(A.col == b.dim_v);
   if (dim_v != A.row) {
      if (mem_v != 0) delete[] ++mem_v;
      dim_v = A.row;
      mem_v = new Float[dim_v] - 1;
   }
   Float* p_this = mem_v + 1;
   Float* pb;
   Float sum;
   for (int i = 1; i <= A.row; i++) {
      pb = b.mem_v;
      sum = 0;
      for (int j = 1; j <= A.col; j++)
	 sum += A[i][j] * *++pb;
      *p_this = sum;
      p_this++;
   }
   return *this;
}

Vec& Vec::mult_trA_b(const Mat& A, const Vec& b) {
   assert(A.row == b.dim_v);
   if (dim_v != A.col) {
      if (mem_v != 0) delete[] ++mem_v;
      dim_v = A.col;
      mem_v = new Float[dim_v] - 1;
   }
   Float* p_this = mem_v + 1;
   Float* pb;
   Float sum;
   for (int i = 1; i <= A.col; i++) {
      pb = b.mem_v;
      sum = 0;
      for (int j = 1; j <= A.row; j++)
	 sum += A[j][i] * *++pb;
      *p_this++ = sum;
   }
   return *this;
}

Vec& Vec::solve(SVD& svdec, const Vec& b) {
   svdec.svd();
   Vec tmp;
   tmp.mult_trA_b(svdec.U, b);
   for (int i=1; i<=svdec.inv_W.dim(); i++)
      tmp[i] *= svdec.inv_W[i];
   (*this).mult_Ab(svdec.V, tmp);
   return *this;
}

// ---------  Vec friends  --------------------------------------------


std::istream& operator>>(std::istream& istr, Vec& A) {
   int r = 0;
   istr >> r;
   if (r != A.dim_v) {
      if (A.mem_v != 0)
	 delete[] ++A.mem_v;
      A.mem_v = new Float[r] - 1;
      A.dim_v = r;
   }

   for (int i = 1; i <= r; i++)
      istr >> A[i];

   return istr;
}

std::ostream& operator<<(std::ostream& ostr, const Vec& A) {
   int aw = ostr.width();
   if (aw > 0) aw--;
   ostr << std::setw(0);
   ostr << ' ' << std::setw(aw) << A.dim_v << "\n\n";
   for (int i = 1; i <= A.dim_v; i++) {
      ostr << ' ' << std::setw(aw) << A[i] << "\n";
   }
   return ostr;
}

// --------- Singular Value Decompiosition methods ----------------------


inline static Float Fabs(Float x) { return (x >= (Float)0) ? x : -x ; }

void SVD::set_inv_W() {
   if (W_tol == 0) {   // if not defined set W_tol to 1000*comp_epsilon 

      Float  eps, eps_1, eps_min, eps_max, sum;
      const Float one = 1;

      eps_min = 0;
      eps_max = eps = 1e-5;
      do {
	 eps_1 = eps;
	 eps = (eps_min + eps_max) / 2;
	 sum = one + eps;
	 if (sum == one)
	    eps_min = eps;
	 else
	    eps_max = eps;
      } while (Fabs(eps - eps_1)/eps > 0.1);
      W_tol = 1000*eps;
   }
   for (int i=1; i<=inv_W.dim(); i++)
      if (Fabs(W[i]) >= W_tol)
	 inv_W[i] = 1/W[i];
      else
	 inv_W[i] = 0;
}

/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                                                   !!!
!!! The following source for SVD methods is based on  !!!
!!! the SVD functions written by C. R. Birchenhall    !!!
!!! (see the original copyright notice)               !!!
!!!                                                   !!!
!!! You can compare my derivative text with original  !!!
!!! code in the file matsvd.c, which should be        !!!
!!! distributed with this source.                     !!!
!!!                                                   !!!
!!!                  Ales Cepek  <cepek@fsv.cvut.cz>  !!!
!!!                                                   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/

/**************************************************/
/*      matsvd.c source for SVD functions         */
/**************************************************/


/**************************************************/
/*            MatClass Source File                */
/*       Copyright of C. R. Birchenhall           */
/*       University of Manchester, UK.            */
/*   MatClass is freeware. This file should be    */
/* made freely available to users of any software */
/* whose creation is wholly or partly dependent   */
/*                on this file.                   */
/**************************************************/


#include <math.h>


static Float at, bt, ct ;

inline Float PYTHAG( Float a, Float b )
{
  return
    ( ( at = Fabs(a) )  > ( bt = Fabs(b) ) )     ?
    ( ct = bt/at, at * sqrt( (Float)1 + ct * ct ) )   :
    ( bt ? ( ct = at/bt, bt * sqrt( (Float)1 + ct * ct ) ) : (Float)0 )
  ;
} // inline pythag



void SVD::svd()
{
   if (decomposed)
      return;

   int  flag, i, its, j, jj, k, l, nm ;
   Float   c, f, h, s, x, y, z ;
   Float   anorm=(Float)0, g=(Float)0, scale=(Float)0, r ;
   int  imin;
   
   Vec rv1(n) ;

   for ( i = 1 ; i <= n ; i++ ) {   // Householder reduction to bidiagonal form

      l = i+1;
      rv1[i] = scale * g ;
      g = s = scale = (Float)0 ;
      if ( i <= m ) {
         for ( k = i ; k <= m ; k++ )
            scale += Fabs( U[k][i] ) ;
         if ( scale != (Float)0 ) {
            for ( k = i ; k <= m ; k++ ) {
               U[k][i] /= scale;
               s += U[k][i] * U[k][i] ;
            } // for k

            f = U[i][i] ;
            g = sqrt(s) ;
            if ( f >= (Float)0 )
               g = - g ;
            h = f * g - s ;
            U[i][i] = f - g ;
            //??? if ( i != n ) {

               for ( j = l ; j <= n ; j++ ) {
                  for ( s = (Float)0 , k = i ; k <= m ; k++ )
                      s += U[k][i] * U[k][j] ;
                  f = s / h ;
                  for ( k = i ; k <= m ; k++ )
                      U[k][j] += f * U[k][i] ;
               } // for j

	    //??? } // if i != n

            for ( k = i ; k <= m ; k++ )
                U[k][i] *= scale ;
         } // if scale

      } // if i <= m

      W[i] = scale * g ;
      g = s = scale = (Float)0 ;
      if ( i <= m && i != n ) {
         for ( k = l ; k <= n ; k++ )
             scale += Fabs( U[i][k] ) ;
         if ( scale != (Float)0 ) {
            for ( k = l ; k <= n ; k++ ) {
               U[i][k] /= scale ;
               s += U[i][k] * U[i][k] ;
            } // for k

            f = U[i][l] ;
            g = sqrt(s) ;
            if ( f >= (Float)0 )
               g = -g ;
            h = f * g - s ;
            U[i][l] = f - g ;
            for ( k = l ; k <= n ; k++ )
                rv1[k] = U[i][k] / h ;
            //??? if ( i != m ) {

               for ( j = l ; j <= m ; j++ ) {
                  for ( s = (Float)0 , k = l ; k <= n ; k++ )
                      s += U[j][k] * U[i][k];
                  for ( k = l ; k <= n ; k++ )
                      U[j][k] += s * rv1[k];
               } // for j

	    //??? } // if i != m

            for ( k = l ; k <= n ; k++ )
                U[i][k] *= scale;
         } // if scale

      } // if i != m && i != n

      r = Fabs( W[i] ) + Fabs( rv1[i] ) ;
      if ( r > anorm )
         anorm = r ;
   } // for i

   for ( i = n ; i >= 1 ; i-- ) {
      //??? l = i + 1 ;

      if ( i < n ) {
         if ( g != (Float)0 ) {
            for ( j = l ; j <= n ; j++ )
               V[j][i] = ( U[i][j] / U[i][l] ) / g ;
               // double division to reduce underflow

            for ( j = l ; j <= n ; j++ ) {
               for ( s = (Float)0, k =l ; k <= n ; k++ )
                  s += U[i][k] * V[k][j] ;
               for ( k = l ; k <= n ; k++ )
                  V[k][j] += s * V[k][i] ;
            } // for j

         } // if g

         for ( j = l ; j <= n ; j++ )
            V[i][j] = V[j][i] = (Float)0 ;
      } // if i < n

      V[i][i] = (Float)1 ;
      g = rv1[i] ;
      l = i;
   } // for i

   imin = (m < n) ? m : n;
   for ( i = imin ; i >= 1 ; i-- ) {
      l = i + 1 ;
      g = W[i] ;
      //??? if ( i < n ) {

         for ( j = l ; j <= n ; j++ )
            U[i][j]=(Float)0 ;
      //???} // if

      if ( g != (Float)0 ) {
         g = (Float)1 / g ;
         //??? if ( i != n ) {

            for ( j = l ; j <= n ; j++ ) {
               for ( s = (Float)0 , k = l ; k <= m ; k++ )
                  s += U[k][i] * U[k][j];
               f = ( s / U[i][i] ) * g ;
               for ( k = i ; k <= m ; k++ )
                  U[k][j] += f * U[k][i] ;
            } // for j

	 //??? } // if i != n

         for ( j = i ; j <= m ; j++ )
            U[j][i] *= g ;
      } else {
         for ( j = i ; j <= m ; j++ )
            U[j][i] = (Float)0 ;
      } // else

      U[i][i] += (Float)1 ;
   } // for i


   for ( k = n ; k >= 1 ; k-- ) {
      for ( its = 1 ; its <= 30 ; its++ ) {
         flag = 1 ;
         for ( l = k ; l >= 1 ; l-- ) {
            nm = l - 1 ;
            if ( (Fabs( rv1[l] ) + anorm) == anorm) {
               flag = 0 ;
               break;
            } // if

            if ( (Fabs( W[nm] ) + anorm) == anorm )
               break;
         } // for l

         if ( flag ) {
            c = (Float)0 ;
            s = (Float)1;
            for ( i = l ; i <= k ; i++ ) {
               f = s * rv1[i] ;
	       rv1[i] = c * rv1[i];   // ???

               if ( (Fabs(f) + anorm) != anorm ) {
                  g = W[i] ;
                  h = PYTHAG(f,g) ;
                  W[i] = h ;
                  h = (Float)1 / h ;
                  c = g * h ;
                  s = ( -f * h ) ;
                  for ( j = 1 ; j <= m ; j++ ) {
                     y = U[j][nm] ;
                     z = U[j][i] ;
                     U[j][nm] = y * c + z * s ;
                     U[j][i] = z * c - y * s ;
                  } // for j

               } // if Fabs(f)

            } // for i

         } // if flag

         z = W[k] ;
         if ( l == k ) {
            if ( z < (Float)0 ) {
               W[k] = -z ;
               for ( j = 1 ; j <= n ; j++ )
                  V[j][k] = (-V[j][k]) ;
            } // if z < 0

            break;
         } // if l == k

		 //char * msgerr= "no convergence in 30 SVD iterations";
         if ( its == 30 ) return; //perror(msgerr);  
         x = W[l];
         nm = k - 1 ;
         y = W[nm] ;
         g = rv1[nm] ;
         h = rv1[k] ;
         f = ( (y-z)*(y+z) + (g-h)*(g+h) ) / ( (Float)2 * h * y ) ;
         g = PYTHAG(f,(Float)1) ;
         r = ( f >= (Float)0 ? g : - g ) ;
         f= ( (x-z)*(x+z) + h * ( ( y / ( f + r ) ) - h ) ) / x ;
         c = s = (Float)1 ;
         for ( j = l ; j <= nm ; j++ ) {
            i = j + 1 ;
            g = rv1[i] ;
            y = W[i] ;
            h = s * g ;
            g = c * g ;
            z = PYTHAG(f,h) ;
            rv1[j] = z ;
            c = f / z ;
            s = h / z ;
            f = x * c + g * s ;
            g = g * c - x * s ;
            h = y * s ;
            y = y * c ;
            for ( jj = 1 ; jj <= n ; jj++ ) {
               x = V[jj][j] ;
               z = V[jj][i];
               V[jj][j] = x * c + z * s ;
               V[jj][i]= z * c - x * s ;
            } // for jj

            z = PYTHAG(f,h) ;
            W[j] = z ;
            if (z) {
               z = (Float)1 / z ;
               c = f * z ;
               s = h * z ;
            } // if

            f = ( c * g ) + ( s * y ) ;
            x = ( c * y ) - ( s * g ) ;
            for ( jj = 1 ; jj <= m ; jj++ ) {
               y = U[jj][j] ;
               z = U[jj][i] ;
               U[jj][j] = y * c + z * s ;
               U[jj][i] = z * c - y * s ;
            } // for jj

         } // for j

         rv1[l] = (Float)0 ;
         rv1[k] = f ;
         W[k] = x ;
      } // for its

   } // for k


   decomposed = 1;
   set_inv_W();

} // svdcmp


/*
int svdBackSub( const Mat& a, const Mat& w, 
                  const Mat& v, const Mat& b, 
                  Mat& x, matError& error )
*/
/**************************************************************
   Assumes a, w and v are svd decomp of A. Solves Ax=b.
   Ignores components with zero singular values.
   Overwrites x with solution.
**************************************************************/
/*
{
   static char *mName = "svdBackSub" ;
   matFunc func( mName ) ; a.debugInfo( func ) ;
   error = NOERROR ;
   int nc = a.cols(), nr = a.rows() ;
   if ( nc != w.rows() || nc != v.rows() || 
        nc != v.cols() || nr != b.rows() ) {
      error = NEDIM ;
      return 0 ;
   } // if
   if ( b.cols() != 1 ) {
      error = NOTVC ;
      return 1 ;
   } // if
   // form tmp = w^-1 a'b, passing over zeros in w 
   Mat tmp( nc ) ;
   refMat m( a ) ;
   Float r ;
   int i ;
   for ( i = 1 ; i <= nc ; i++ ) {
      r = (Float)0 ;
      if ( W[i] != (Float)0 ) {
         m.refCol(i) ;
         r = m.inner(b) ;
	 r /= W[i] ;
      } // if
      tmp(i) = r ;
   } // if
   x.multOf( v, tmp ) ;
   return OK ;
} // svdBackSub
*/

SVD& SVD::reset(const Mat& A)
{
  const int R = A.rows();
  const int C = A.cols();
  if (R != U.rows() || C != U.cols()) {
    m = R;
    n = C;
    W.reset(C);
    V.reset(C, C);
    inv_W.reset(C);
  }
  U = A;
  decomposed = 0;
  
  return *this;
}

SVD& SVD::reset(const Mat& A, const Vec& w)
{
  reset(A);
  for (int i = 1; i <= A.rows(); i++) 
    for (int j = 1; j <= A.cols(); j++)
      U[i][j] *= w[i];

  return *this;
}

SVD& SVD::clear()
{
  m = 0;
  n = 0;
  W.reset();
  V.reset();
  inv_W.reset();
  U.reset();
  decomposed = 0;

  return *this;
}

Float SVD::cov_xx(int i, int j)
{
  assert(1 <= i && i <= n && 1 <= j && j <= n);
  // A = U*W*trans(V)

  // Covariance = V * inv_W * trans(inv_W) * trans(V);

  Float c = 0;   
  for (int k = 1; k <= n; k++)
    c += V[i][k] * inv_W[k] * inv_W[k] * V[j][k];
  return c;
}

Float SVD::cov_bb(int i, int j)
{
  assert(1 <= i && i <= m && 1 <= j && j <= m);
  // A = U*W*trans(V)

  // Covariance = U * trans(U);

  Float c = 0;   
  for (int k = 1; k <= n; k++)
    if (inv_W[k] != 0)            // !!!!!!!!!!!!!!!!!!!!!!!

      c += U[i][k] * U[j][k];
  return c;
}

Float SVD::cov_bx(int i, int j)
{
  assert(1 <= i && i <= m && 1 <= j && j <= n);
  // A = U*W*trans(V)

  // Covariance = U * trans(V);

  Float c = 0;   
  for (int k = 1; k <= n; k++)
    c += U[i][k] * inv_W[k] * V[j][k];
  return c;
}

SVD::SVD(const Mat& A) : U(A), m(A.rows()), n(A.cols()), W(n), V(n, n),
                       decomposed(0), W_tol(0), inv_W(n)
{}
