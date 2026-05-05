#include "BCGSTAB.h"
#include <cmath>

std::vector<double> BCGSTAB::SolutionWithDiagonalConditioning( SLE sle )
{
   const int n = sle.matrix.n;
   const double fNorm2 = SLEAssistant::DotProduct( sle.f, sle.f );
   if ( fNorm2 < 1e-30 )
   {
      std::cout << "Residual: 0 k: 0\n";
      return sle.x;
   }
   const double stopCrit = sle.eps * sle.eps * fNorm2;

   std::vector<double> r( n ), r0( n ), p( n ), v( n ), s( n ), t( n ), z( n ), y( n );

   // r = f - A*x0
   SLEAssistant::MatrixVectorProduct( sle.matrix, sle.x, r );
   for ( int i = 0; i < n; ++i )
      r[i] = sle.f[i] - r[i];

   r0 = r;
   p = r;

   double rho = 1.0, alpha = 1.0, omega = 1.0;
   double residual2 = SLEAssistant::DotProduct( r, r );
   int k = 0;

   for ( ; k < sle.maxIter && residual2 > stopCrit; ++k )
   {
      const double rhoNew = SLEAssistant::DotProduct( r0, r );
      if ( std::abs( rhoNew ) < 1e-30 ) break;

      const double beta = ( rhoNew / rho ) * ( alpha / omega );
      rho = rhoNew;

      // p = r + beta * (p - omega * v)
      for ( int i = 0; i < n; ++i )
         p[i] = r[i] + beta * ( p[i] - omega * v[i] );

      // z = M^{-1} * p (диагональное предобуславливание)
      SLEAssistant::DiagonalPreconditionMSG( sle.matrix, p, z );

      // v = A * z
      SLEAssistant::MatrixVectorProduct( sle.matrix, z, v );

      const double r0v = SLEAssistant::DotProduct( r0, v );
      if ( std::abs( r0v ) < 1e-30 ) break;
      alpha = rho / r0v;

      // s = r - alpha * v
      for ( int i = 0; i < n; ++i )
         s[i] = r[i] - alpha * v[i];

      // Проверка на ранний выход
      const double sNorm2 = SLEAssistant::DotProduct( s, s );
      if ( sNorm2 < stopCrit )
      {
         for ( int i = 0; i < n; ++i )
            sle.x[i] += alpha * z[i];
         residual2 = sNorm2;
         break;
      }

      // y = M^{-1} * s
      SLEAssistant::DiagonalPreconditionMSG( sle.matrix, s, y );

      // t = A * y
      SLEAssistant::MatrixVectorProduct( sle.matrix, y, t );

      const double tt = SLEAssistant::DotProduct( t, t );
      if ( std::abs( tt ) < 1e-30 ) break;
      omega = SLEAssistant::DotProduct( t, s ) / tt;

      // x = x + alpha * z + omega * y
      for ( int i = 0; i < n; ++i )
         sle.x[i] += alpha * z[i] + omega * y[i];

      // r = s - omega * t
      for ( int i = 0; i < n; ++i )
         r[i] = s[i] - omega * t[i];

      residual2 = SLEAssistant::DotProduct( r, r );

      if ( std::abs( omega ) < 1e-30 ) break;
   }

   std::cout << "Residual: " << std::sqrt( residual2 / fNorm2 ) << " k: " << k << "\n";
   return sle.x;
}
