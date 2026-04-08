#include "MSG.h"

std::vector<double> MSG::SolutionWithoutConditioning( SLE sle )
{
   int    k = 0;
   double residual = 1;
   double a = 0;
   double b = 0;
   double fNorm = sqrt( SLEAssistant::DotProduct( sle.f, sle.f ) );
   double condition = ( sle.eps * fNorm ) * ( sle.eps * fNorm );

   std::vector<double> ax( sle.matrix.n, 0 );
   std::vector<double> r( sle.matrix.n, 0 );
   std::vector<double> z( sle.matrix.n, 0 );
   std::vector<double> az( sle.matrix.n, 0 );

   SLEAssistant::MatrixVectorProduct( sle.matrix, sle.x, ax );

   for ( int i = 0; i < sle.matrix.n; ++i )
   {
      r[i] = sle.f[i] - ax[i];
      z[i] = r[i];
   }

   for ( ; k < sle.maxIter and residual > condition; ++k )
   {
      SLEAssistant::MatrixVectorProduct( sle.matrix, z, az );
      a = SLEAssistant::DotProduct( r, r ) / SLEAssistant::DotProduct( az, z );
      b = SLEAssistant::DotProduct( r, r );

      for ( int i = 0; i < sle.matrix.n; ++i )
      {
         sle.x[i] = sle.x[i] + a * z[i];
         r[i] = r[i] - a * az[i];
      }

      b = SLEAssistant::DotProduct( r, r ) / b;
      for ( int i = 0; i < sle.matrix.n; ++i )
      {
         z[i] = r[i] + b * z[i];
      }

      residual = SLEAssistant::DotProduct( r, r );
   }

   std::cout << "Residual: " << sqrt( residual ) / fNorm << " k: " << k << "\n";
   return sle.x;
}

std::vector<double> MSG::SolutionWithDiagonalConditioning( SLE sle )
{
   int    k = 0;
   double residual = 1;
   double a = 0;
   double b = 0;
   double fNorm = sqrt( SLEAssistant::DotProduct( sle.f, sle.f ) );
   double condition = ( sle.eps * fNorm ) * ( sle.eps * fNorm );

   std::vector<double> ax( sle.matrix.n, 0 );
   std::vector<double> r( sle.matrix.n, 0 );
   std::vector<double> z( sle.matrix.n, 0 );
   std::vector<double> az( sle.matrix.n, 0 );
   std::vector<double> rD( sle.matrix.n, 0 );
   SLEAssistant::MatrixVectorProduct( sle.matrix, sle.x, ax );

   for ( int i = 0; i < sle.matrix.n; ++i )
   {
      r[i] = sle.f[i] - ax[i];
      z[i] = r[i];
   }

   SLEAssistant::DiagonalPreconditionMSG( sle.matrix, z, z );
   SLEAssistant::DiagonalPreconditionMSG( sle.matrix, r, rD );

   for ( ; k < sle.maxIter and residual > condition; ++k )
   {
      SLEAssistant::MatrixVectorProduct( sle.matrix, z, az );
      a = SLEAssistant::DotProduct( rD, r ) / SLEAssistant::DotProduct( az, z );
      b = SLEAssistant::DotProduct( rD, r );

      for ( int i = 0; i < sle.matrix.n; ++i )
      {
         sle.x[i] = sle.x[i] + a * z[i];
         r[i] = r[i] - a * az[i];
      }
      SLEAssistant::DiagonalPreconditionMSG( sle.matrix, r, rD );

      b = SLEAssistant::DotProduct( rD, r ) / b;
      for ( int i = 0; i < sle.matrix.n; ++i )
         z[i] = rD[i] + b * z[i];

      residual = SLEAssistant::DotProduct( r, r );
   }

   std::cout << "Residual: " << sqrt( residual ) / fNorm << " k: " << k << "\n";
   return sle.x;
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