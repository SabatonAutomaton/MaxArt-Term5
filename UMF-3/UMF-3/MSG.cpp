#include "MSG.h"
#include <cmath>

namespace
{
   std::vector<double> SolveBiCGStab( SLE sle, bool diagonalPreconditioner )
   {
      const int n = sle.matrix.n;
      const double tiny = 1e-30;

      std::vector<double> x = sle.x;
      std::vector<double> ax( n, 0.0 );
      std::vector<double> r( n, 0.0 );
      std::vector<double> rHat( n, 0.0 );
      std::vector<double> p( n, 0.0 );
      std::vector<double> pHat( n, 0.0 );
      std::vector<double> v( n, 0.0 );
      std::vector<double> s( n, 0.0 );
      std::vector<double> sHat( n, 0.0 );
      std::vector<double> t( n, 0.0 );

      SLEAssistant::MatrixVectorProduct( sle.matrix, x, ax );
      for ( int i = 0; i < n; ++i )
      {
         r[i] = sle.f[i] - ax[i];
      }
      rHat = r;

      const double fNorm2 = SLEAssistant::DotProduct( sle.f, sle.f );
      if ( fNorm2 < tiny )
      {
         std::cout << "Residual: 0 k: 0\n";
         return x;
      }

      const double stop2 = ( sle.eps * sle.eps ) * fNorm2;
      double residual2 = SLEAssistant::DotProduct( r, r );
      if ( residual2 <= stop2 )
      {
         std::cout << "Residual: " << std::sqrt( residual2 / fNorm2 ) << " k: 0\n";
         return x;
      }

      double rhoOld = 1.0;
      double alpha = 1.0;
      double omega = 1.0;
      int k = 0;

      for ( ; k < sle.maxIter; ++k )
      {
         const double rhoNew = SLEAssistant::DotProduct( rHat, r );
         if ( std::abs( rhoNew ) < tiny )
         {
            break;
         }

         if ( k == 0 )
         {
            p = r;
         }
         else
         {
            const double beta = ( rhoNew / rhoOld ) * ( alpha / omega );
            for ( int i = 0; i < n; ++i )
            {
               p[i] = r[i] + beta * ( p[i] - omega * v[i] );
            }
         }

         if ( diagonalPreconditioner )
         {
            SLEAssistant::DiagonalPreconditionMSG( sle.matrix, p, pHat );
         }
         else
         {
            pHat = p;
         }

         SLEAssistant::MatrixVectorProduct( sle.matrix, pHat, v );
         const double denomAlpha = SLEAssistant::DotProduct( rHat, v );
         if ( std::abs( denomAlpha ) < tiny )
         {
            break;
         }
         alpha = rhoNew / denomAlpha;

         for ( int i = 0; i < n; ++i )
         {
            s[i] = r[i] - alpha * v[i];
         }

         const double sNorm2 = SLEAssistant::DotProduct( s, s );
         if ( sNorm2 <= stop2 )
         {
            for ( int i = 0; i < n; ++i )
            {
               x[i] += alpha * pHat[i];
            }
            residual2 = sNorm2;
            ++k;
            break;
         }

         if ( diagonalPreconditioner )
         {
            SLEAssistant::DiagonalPreconditionMSG( sle.matrix, s, sHat );
         }
         else
         {
            sHat = s;
         }

         SLEAssistant::MatrixVectorProduct( sle.matrix, sHat, t );
         const double tt = SLEAssistant::DotProduct( t, t );
         if ( tt < tiny )
         {
            break;
         }
         omega = SLEAssistant::DotProduct( t, s ) / tt;
         if ( std::abs( omega ) < tiny )
         {
            break;
         }

         for ( int i = 0; i < n; ++i )
         {
            x[i] += alpha * pHat[i] + omega * sHat[i];
            r[i] = s[i] - omega * t[i];
         }

         residual2 = SLEAssistant::DotProduct( r, r );
         if ( residual2 <= stop2 )
         {
            ++k;
            break;
         }

         rhoOld = rhoNew;
      }

      std::cout << "Residual: " << std::sqrt( residual2 / fNorm2 ) << " k: " << k << "\n";
      return x;
   }
}

std::vector<double> MSG::SolutionWithoutConditioning( SLE sle )
{
   return SolveBiCGStab( sle, false );
}

std::vector<double> MSG::SolutionWithDiagonalConditioning( SLE sle )
{
   return SolveBiCGStab( sle, true );
}

std::vector<double> MSG::SolutionWithLUConditioning( SLE sle )
{
   LU     lu( sle.matrix );
   int    k = 0;
   double residual = 1;
   double a = 0;
   double b = 0;
   double fNorm = sqrt( SLEAssistant::DotProduct( sle.f, sle.f ) );
   double condition = ( sle.eps * fNorm ) * ( sle.eps * fNorm );

   std::vector<double> ax( sle.matrix.n, 0 );
   std::vector<double> z( sle.matrix.n, 0 );
   std::vector<double> az( sle.matrix.n, 0 );
   std::vector<double> rLU( sle.matrix.n, 0 );
   std::vector<double> r = sle.f;

   SLEAssistant::LUDirect( lu, r, z );
   SLEAssistant::LUReverse( lu, z, z );

   for ( ; k < sle.maxIter and residual > condition; ++k )
   {
      SLEAssistant::MatrixVectorProduct( sle.matrix, z, az );
      SLEAssistant::LUDirect( lu, r, rLU );
      SLEAssistant::LUReverse( lu, rLU, rLU );
      b = SLEAssistant::DotProduct( rLU, r );
      a = b / SLEAssistant::DotProduct( az, z );

      for ( int i = 0; i < sle.matrix.n; ++i )
      {
         sle.x[i] = sle.x[i] + a * z[i];
         r[i] = r[i] - a * az[i];
      }

      SLEAssistant::LUDirect( lu, r, rLU );
      SLEAssistant::LUReverse( lu, rLU, rLU );
      b = SLEAssistant::DotProduct( rLU, r ) / b;

      for ( int i = 0; i < sle.matrix.n; ++i )
         z[i] = rLU[i] + b * z[i];

      residual = SLEAssistant::DotProduct( r, r );
   }

   std::cout << "Residual: " << sqrt( residual ) / fNorm << " k: " << k << "\n";
   return sle.x;
}

std::vector<double> MSG::SolutionWithLLtConditioning( SLE sle )
{
   LLt    llt( sle.matrix );
   int    k = 0;
   double residual = 1;
   double a = 0;
   double b = 0;
   double fNorm = sqrt( SLEAssistant::DotProduct( sle.f, sle.f ) );
   double condition = ( sle.eps * fNorm ) * ( sle.eps * fNorm );

   std::vector<double> ax( sle.matrix.n, 0 );
   std::vector<double> z( sle.matrix.n, 0 );
   std::vector<double> az( sle.matrix.n, 0 );
   std::vector<double> rLLt( sle.matrix.n, 0 );
   std::vector<double> r = sle.f;

   SLEAssistant::LLtDirect( llt, r, z );
   SLEAssistant::LLtReverse( llt, z, z );

   bool nan = false;

   for ( ; k < sle.maxIter and residual > condition; ++k )
   {
      SLEAssistant::MatrixVectorProduct( sle.matrix, z, az );
      SLEAssistant::LLtDirect( llt, r, rLLt );
      SLEAssistant::LLtReverse( llt, rLLt, rLLt );
      b = SLEAssistant::DotProduct( rLLt, r );
      a = b / SLEAssistant::DotProduct( az, z );
      if ( std::isnan( a ) ) break;

      for ( int i = 0; i < sle.matrix.n; ++i )
      {
         sle.x[i] = sle.x[i] + a * z[i];
         r[i] = r[i] - a * az[i];
      }

      SLEAssistant::LLtDirect( llt, r, rLLt );
      SLEAssistant::LLtReverse( llt, rLLt, rLLt );
      b = SLEAssistant::DotProduct( rLLt, r ) / b;

      for ( int i = 0; i < sle.matrix.n; ++i )
         z[i] = rLLt[i] + b * z[i];

      residual = SLEAssistant::DotProduct( r, r );
   }

   std::cout << "Residual: " << sqrt( residual ) / fNorm << " k: " << k << "\n";
   return sle.x;
}
