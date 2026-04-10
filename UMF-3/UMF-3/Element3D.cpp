#include "Element3D.h"
#include <cmath>
#include <stdexcept>

double Element3D::Det3( const std::array<std::array<double, 3>, 3> &m )
{
   return
      m[0][0] * ( m[1][1] * m[2][2] - m[1][2] * m[2][1] ) -
      m[0][1] * ( m[1][0] * m[2][2] - m[1][2] * m[2][0] ) +
      m[0][2] * ( m[1][0] * m[2][1] - m[1][1] * m[2][0] );
}

std::array<std::array<double, 3>, 3> Element3D::Inverse3( const std::array<std::array<double, 3>, 3> &m )
{
   const double det = Det3( m );
   if ( std::abs( det ) < 1e-20 )
   {
      throw std::runtime_error( "Degenerate hexahedron: Jacobian determinant is near zero" );
   }

   std::array<std::array<double, 3>, 3> inv {};
   inv[0][0] = ( m[1][1] * m[2][2] - m[1][2] * m[2][1] ) / det;
   inv[0][1] = ( m[0][2] * m[2][1] - m[0][1] * m[2][2] ) / det;
   inv[0][2] = ( m[0][1] * m[1][2] - m[0][2] * m[1][1] ) / det;

   inv[1][0] = ( m[1][2] * m[2][0] - m[1][0] * m[2][2] ) / det;
   inv[1][1] = ( m[0][0] * m[2][2] - m[0][2] * m[2][0] ) / det;
   inv[1][2] = ( m[0][2] * m[1][0] - m[0][0] * m[1][2] ) / det;

   inv[2][0] = ( m[1][0] * m[2][1] - m[1][1] * m[2][0] ) / det;
   inv[2][1] = ( m[0][1] * m[2][0] - m[0][0] * m[2][1] ) / det;
   inv[2][2] = ( m[0][0] * m[1][1] - m[0][1] * m[1][0] ) / det;
   return inv;
}

double Element3D::Approx( const std::vector<double> &val, const std::array<double, 3> &rst ) const
{
   if ( val.size( ) < 8 )
   {
      return 0.0;
   }
   double res = 0.0;
   for ( int i = 0; i < 8; ++i )
   {
      res += val[i] * basis.Basis( i, rst );
   }
   return res;
}

void Element3D::BuildLocalMatrix( )
{
   localMatrix.assign( 16, std::vector<double>( 16, 0.0 ) );

   const std::vector<std::vector<double>> quadratures = quad.gaussHex2( );
   for ( const auto &qv : quadratures )
   {
      const std::array<double, 3> q = { qv[0], qv[1], qv[2] };
      std::array<std::array<double, 3>, 3> J {};
      std::array<std::array<double, 3>, 8> gradRef {};
      for ( int aNode = 0; aNode < 8; ++aNode )
      {
         gradRef[aNode] = basis.GradBasis( aNode, q );
         J[0][0] += coords[aNode][0] * gradRef[aNode][0];
         J[0][1] += coords[aNode][0] * gradRef[aNode][1];
         J[0][2] += coords[aNode][0] * gradRef[aNode][2];
         J[1][0] += coords[aNode][1] * gradRef[aNode][0];
         J[1][1] += coords[aNode][1] * gradRef[aNode][1];
         J[1][2] += coords[aNode][1] * gradRef[aNode][2];
         J[2][0] += coords[aNode][2] * gradRef[aNode][0];
         J[2][1] += coords[aNode][2] * gradRef[aNode][1];
         J[2][2] += coords[aNode][2] * gradRef[aNode][2];
      }

      const double detJ = Det3( J );
      const auto invJ = Inverse3( J );
      const double wj = std::abs( detJ ) * qv[3];
      const double lambda = Approx( lambdaNodes, q );
      const double w = Approx( wNodes, q );
      const double sigma = Approx( sigmaNodes, q );
      const double xi = Approx( xiNodes, q );

      std::array<std::array<double, 3>, 8> gradPhys {};
      for ( int i = 0; i < 8; ++i )
      {
         for ( int c = 0; c < 3; ++c )
         {
            gradPhys[i][c] =
               invJ[0][c] * gradRef[i][0] +
               invJ[1][c] * gradRef[i][1] +
               invJ[2][c] * gradRef[i][2];
         }
      }

      for ( int i = 0; i < 8; ++i )
      {
         for ( int j = 0; j < 8; ++j )
         {
            const double dot =
               gradPhys[i][0] * gradPhys[j][0] +
               gradPhys[i][1] * gradPhys[j][1] +
               gradPhys[i][2] * gradPhys[j][2];
            const double psiI = basis.Basis( i, q );
            const double psiJ = basis.Basis( j, q );
            const double a = wj * ( lambda * dot - w * w * xi * psiI * psiJ );
            const double b = wj * w * sigma * psiI * psiJ;

            localMatrix[2 * i][2 * j] += a;
            localMatrix[2 * i][2 * j + 1] += -b;
            localMatrix[2 * i + 1][2 * j] += b;
            localMatrix[2 * i + 1][2 * j + 1] += a;
         }
      }
   }
}

void Element3D::BuildLocalB( )
{
   localB.assign( 16, 0.0 );

   const std::vector<std::vector<double>> quadratures = quad.gaussHex2( );
   for ( const auto &qv : quadratures )
   {
      const std::array<double, 3> q = { qv[0], qv[1], qv[2] };
      std::array<std::array<double, 3>, 3> J {};
      for ( int aNode = 0; aNode < 8; ++aNode )
      {
         const auto g = basis.GradBasis( aNode, q );
         J[0][0] += coords[aNode][0] * g[0];
         J[0][1] += coords[aNode][0] * g[1];
         J[0][2] += coords[aNode][0] * g[2];
         J[1][0] += coords[aNode][1] * g[0];
         J[1][1] += coords[aNode][1] * g[1];
         J[1][2] += coords[aNode][1] * g[2];
         J[2][0] += coords[aNode][2] * g[0];
         J[2][1] += coords[aNode][2] * g[1];
         J[2][2] += coords[aNode][2] * g[2];
      }

      const double detJ = Det3( J );
      const double wj = std::abs( detJ ) * qv[3];
      const double fs = Approx( fsNodes, q );
      const double fc = Approx( fcNodes, q );
      for ( int i = 0; i < 8; ++i )
      {
         const double Ni = basis.Basis( i, q );
         localB[2 * i] += wj * fs * Ni;
         localB[2 * i + 1] += wj * fc * Ni;
      }
   }
}
