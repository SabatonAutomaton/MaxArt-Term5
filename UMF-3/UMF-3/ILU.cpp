#include "ILU.h"
#include <cmath>
#include <stdexcept>
#include <string> 

ILU::ILU( const SparseMatrix &matrix )
{
   n = matrix.n;
   di = matrix.di;
   ig = matrix.ig;
   jg = matrix.jg;
   ggl = matrix.ggl;
   ggu = matrix.ggu;

   for ( int i = 0; i < n; ++i )
   {
      // ── Обработка внедиагональных элементов строки i ──────────────────
      for ( int idx = ig[i]; idx < ig[i + 1]; ++idx )
      {
         const int j = jg[idx]; // j < i  (нижний треугольник)

         // L(i,j) = (A(i,j) - sum_{m: m<j, (i,m) и (j,m) в профиле} L(i,m)*U(m,j)) / U(j,j)
         double sumL = 0.0;
         int p = ig[i];   // бежим по строке i
         int q = ig[j];   // бежим по строке j
         // оба указателя идут по столбцам < j
         while ( p < idx && q < ig[j + 1] )
         {
            const int cp = jg[p];
            const int cq = jg[q];
            if ( cp == cq ) { sumL += ggl[p] * ggu[q]; ++p; ++q; }
            else if ( cp < cq ) { ++p; }
            else { ++q; }
         }
         ggl[idx] = ( ggl[idx] - sumL ) / di[j];   // di[j] = U(j,j)

         // U(j,i) = A(j,i) - sum_{m: m<j, (i,m) и (j,m) в профиле} L(j,m)*U(m,i)
         // Позиция U(j,i) в верхнем треугольнике — та же ячейка idx, но в ggu
         double sumU = 0.0;
         p = ig[i];
         q = ig[j];
         while ( p < idx && q < ig[j + 1] )
         {
            const int cp = jg[p];
            const int cq = jg[q];
            if ( cp == cq ) { sumU += ggu[p] * ggl[q]; ++p; ++q; }
            else if ( cp < cq ) { ++p; }
            else { ++q; }
         }
         ggu[idx] = ggu[idx] - sumU;
         // fill-in вне профиля молча отбрасывается — это и есть ILU(0)
      }

      // ── Диагональный элемент U(i,i) ───────────────────────────────────
      double sumD = 0.0;
      for ( int idx = ig[i]; idx < ig[i + 1]; ++idx )
         sumD += ggl[idx] * ggu[idx];
      di[i] -= sumD;

      if ( std::abs( di[i] ) < 1e-20 )
         throw std::runtime_error( "ILU(0): near-zero diagonal pivot at row " + std::to_string( i ) );
   }
}