#include "Basis3D.h"

double Basis3D::Basis1( const std::array<double, 3> &p ) const
{
   const double r = p[0], s = p[1], t = p[2];
   return ( 1.0 - r ) * ( 1.0 - s ) * ( 1.0 - t );
}

double Basis3D::Basis2( const std::array<double, 3> &p ) const
{
   const double r = p[0], s = p[1], t = p[2];
   return r * ( 1.0 - s ) * ( 1.0 - t );
}

double Basis3D::Basis3( const std::array<double, 3> &p ) const
{
   const double r = p[0], s = p[1], t = p[2];
   return r * s * ( 1.0 - t );
}

double Basis3D::Basis4( const std::array<double, 3> &p ) const
{
   const double r = p[0], s = p[1], t = p[2];
   return ( 1.0 - r ) * s * ( 1.0 - t );
}

double Basis3D::Basis5( const std::array<double, 3> &p ) const
{
   const double r = p[0], s = p[1], t = p[2];
   return ( 1.0 - r ) * ( 1.0 - s ) * t;
}

double Basis3D::Basis6( const std::array<double, 3> &p ) const
{
   const double r = p[0], s = p[1], t = p[2];
   return r * ( 1.0 - s ) * t;
}

double Basis3D::Basis7( const std::array<double, 3> &p ) const
{
   const double r = p[0], s = p[1], t = p[2];
   return r * s * t;
}

double Basis3D::Basis8( const std::array<double, 3> &p ) const
{
   const double r = p[0], s = p[1], t = p[2];
   return ( 1.0 - r ) * s * t;
}

std::array<double, 3> Basis3D::GradBasis1( const std::array<double, 3> &p ) const
{
   const double r = p[0], s = p[1], t = p[2];
   return { -( 1.0 - s ) * ( 1.0 - t ), -( 1.0 - r ) * ( 1.0 - t ), -( 1.0 - r ) * ( 1.0 - s ) };
}

std::array<double, 3> Basis3D::GradBasis2( const std::array<double, 3> &p ) const
{
   const double r = p[0], s = p[1], t = p[2];
   return { ( 1.0 - s ) * ( 1.0 - t ), -r * ( 1.0 - t ), -r * ( 1.0 - s ) };
}

std::array<double, 3> Basis3D::GradBasis3( const std::array<double, 3> &p ) const
{
   const double r = p[0], s = p[1], t = p[2];
   return { s * ( 1.0 - t ), r * ( 1.0 - t ), -r * s };
}

std::array<double, 3> Basis3D::GradBasis4( const std::array<double, 3> &p ) const
{
   const double r = p[0], s = p[1], t = p[2];
   return { -s * ( 1.0 - t ), ( 1.0 - r ) * ( 1.0 - t ), -( 1.0 - r ) * s };
}

std::array<double, 3> Basis3D::GradBasis5( const std::array<double, 3> &p ) const
{
   const double r = p[0], s = p[1], t = p[2];
   return { -( 1.0 - s ) * t, -( 1.0 - r ) * t, ( 1.0 - r ) * ( 1.0 - s ) };
}

std::array<double, 3> Basis3D::GradBasis6( const std::array<double, 3> &p ) const
{
   const double r = p[0], s = p[1], t = p[2];
   return { ( 1.0 - s ) * t, -r * t, r * ( 1.0 - s ) };
}

std::array<double, 3> Basis3D::GradBasis7( const std::array<double, 3> &p ) const
{
   const double r = p[0], s = p[1], t = p[2];
   return { s * t, r * t, r * s };
}

std::array<double, 3> Basis3D::GradBasis8( const std::array<double, 3> &p ) const
{
   const double r = p[0], s = p[1], t = p[2];
   return { -s * t, ( 1.0 - r ) * t, ( 1.0 - r ) * s };
}

double Basis3D::Basis( int i, const std::array<double, 3> &p ) const
{
   switch ( i )
   {
   case 0: return Basis1( p );
   case 1: return Basis2( p );
   case 2: return Basis3( p );
   case 3: return Basis4( p );
   case 4: return Basis5( p );
   case 5: return Basis6( p );
   case 6: return Basis7( p );
   case 7: return Basis8( p );
   default: return 0.0;
   }
}

std::array<double, 3> Basis3D::GradBasis( int i, const std::array<double, 3> &p ) const
{
   switch ( i )
   {
   case 0: return GradBasis1( p );
   case 1: return GradBasis2( p );
   case 2: return GradBasis3( p );
   case 3: return GradBasis4( p );
   case 4: return GradBasis5( p );
   case 5: return GradBasis6( p );
   case 6: return GradBasis7( p );
   case 7: return GradBasis8( p );
   default: return { 0.0, 0.0, 0.0 };
   }
}
