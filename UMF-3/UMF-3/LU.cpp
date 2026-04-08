#include "LU.h"
#include <cmath>
#include <stdexcept>

LU::LU( const SparseMatrix &matrix )
{
   n = matrix.n;
   di = matrix.di;
   ig = matrix.ig;
   jg = matrix.jg;
   ggl = matrix.ggl;
   ggu = matrix.ggu;

   for ( int i = 0; i < n; ++i )
   {
      for ( int idx = ig[i]; idx < ig[i + 1]; ++idx )
      {
         const int j = jg[idx];

         // Вычисление L(i,j): (A(i,j) - sum_{m<j} L(i,m) * U(m,j)) / D(j)
         double sumL = 0.0;
         int p = ig[i];
         int q = ig[j];
         while ( p < idx && q < ig[j + 1] )
         {
            const int colP = jg[p];
            const int colQ = jg[q];
            if ( colP == colQ )
            {
               sumL += ggl[p] * ggu[q];
               ++p;
               ++q;
            }
            else if ( colP < colQ )
            {
               ++p;
            }
            else
            {
               ++q;
            }
         }
         ggl[idx] = ( ggl[idx] - sumL ) / di[j];

         // Вычисление U(j,i): A(j,i) - sum_{m<j} L(j,m) * U(m,i)
         double sumU = 0.0;
         p = ig[i];
         q = ig[j];
         while ( p < idx && q < ig[j + 1] )
         {
            const int colP = jg[p];
            const int colQ = jg[q];
            if ( colP == colQ )
            {
               sumU += ggu[p] * ggl[q];
               ++p;
               ++q;
            }
            else if ( colP < colQ )
            {
               ++p;
            }
            else
            {
               ++q;
            }
         }
         ggu[idx] = ggu[idx] - sumU;
      }

      double sumD = 0.0;
      for ( int idx = ig[i]; idx < ig[i + 1]; ++idx )
      {
         sumD += ggl[idx] * ggu[idx];
      }
      di[i] -= sumD;
      if ( std::abs( di[i] ) < 1e-20 )
      {
         throw std::runtime_error( "LU factorization failed: near-zero diagonal pivot" );
      }
   }
}
