#include "MSG.h"
#include "ILU.h"
#include <cmath>


std::vector<double> MSG::SolutionWithoutConditioning( SLE sle )
{
   const int n = sle.matrix.n;
   const double fNorm2 = SLEAssistant::DotProduct( sle.f, sle.f );
   if ( fNorm2 < 1e-30 )
   {
      std::cout << "Residual: 0 k: 0\n";
      return sle.x;
   }
   const double stopCrit = sle.eps * sle.eps * fNorm2;

   std::vector<double> r( n ), z( n ), p( n ), q( n );

   // r = f - A*x0
   SLEAssistant::MatrixVectorProduct( sle.matrix, sle.x, r );
   for ( int i = 0; i < n; ++i )
      r[i] = sle.f[i] - r[i];

   // без предобусловливателя: z = r
   p = r;
   z = r;

   double rz = SLEAssistant::DotProduct( r, r );
   double residual2 = rz;
   int k = 0;

   for ( ; k < sle.maxIter && residual2 > stopCrit; ++k )
   {
      SLEAssistant::MatrixVectorProduct( sle.matrix, p, q );

      const double pq = SLEAssistant::DotProduct( p, q );
      if ( std::abs( pq ) < 1e-30 ) break;
      const double alpha = rz / pq;

      for ( int i = 0; i < n; ++i )
      {
         sle.x[i] += alpha * p[i];
         r[i] -= alpha * q[i];
      }

      const double rzNew = SLEAssistant::DotProduct( r, r );
      if ( std::abs( rz ) < 1e-30 ) break;
      const double beta = rzNew / rz;
      rz = rzNew;

      for ( int i = 0; i < n; ++i )
         p[i] = r[i] + beta * p[i];

      residual2 = rz;
   }

   std::cout << "Residual: " << std::sqrt( residual2 / fNorm2 ) << " k: " << k << "\n";
   return sle.x;
}

std::vector<double> MSG::SolutionWithDiagonalConditioning( SLE sle )
{
   const int n = sle.matrix.n;
   const double fNorm2 = SLEAssistant::DotProduct( sle.f, sle.f );
   if ( fNorm2 < 1e-30 )
   {
      std::cout << "Residual: 0 k: 0\n";
      return sle.x;
   }
   const double stopCrit = sle.eps * sle.eps * fNorm2;

   std::vector<double> r( n ), z( n ), p( n ), q( n );

   // r = f - A*x0
   SLEAssistant::MatrixVectorProduct( sle.matrix, sle.x, r );
   for ( int i = 0; i < n; ++i )
      r[i] = sle.f[i] - r[i];

   // z = D^{-1} * r
   SLEAssistant::DiagonalPreconditionMSG( sle.matrix, r, z );
   p = z;

   double rz = SLEAssistant::DotProduct( r, z );
   double residual2 = SLEAssistant::DotProduct( r, r );
   int k = 0;

   for ( ; k < sle.maxIter && residual2 > stopCrit; ++k )
   {
      SLEAssistant::MatrixVectorProduct( sle.matrix, p, q );

      const double pq = SLEAssistant::DotProduct( p, q );
      if ( std::abs( pq ) < 1e-30 ) break;
      const double alpha = rz / pq;

      for ( int i = 0; i < n; ++i )
      {
         sle.x[i] += alpha * p[i];
         r[i] -= alpha * q[i];
      }

      // z = D^{-1} * r_new
      SLEAssistant::DiagonalPreconditionMSG( sle.matrix, r, z );

      const double rzNew = SLEAssistant::DotProduct( r, z );
      if ( std::abs( rz ) < 1e-30 ) break;
      const double beta = rzNew / rz;
      rz = rzNew;

      for ( int i = 0; i < n; ++i )
         p[i] = z[i] + beta * p[i];

      residual2 = SLEAssistant::DotProduct( r, r );
   }

   std::cout << "Residual: " << std::sqrt( residual2 / fNorm2 ) << " k: " << k << "\n";
   return sle.x;
}

