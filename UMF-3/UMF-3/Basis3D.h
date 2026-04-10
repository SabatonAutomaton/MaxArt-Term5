#pragma once
#include <array>

class Basis3D
{
public:
   double Basis1( const std::array<double, 3> &p ) const;
   double Basis2( const std::array<double, 3> &p ) const;
   double Basis3( const std::array<double, 3> &p ) const;
   double Basis4( const std::array<double, 3> &p ) const;
   double Basis5( const std::array<double, 3> &p ) const;
   double Basis6( const std::array<double, 3> &p ) const;
   double Basis7( const std::array<double, 3> &p ) const;
   double Basis8( const std::array<double, 3> &p ) const;

   std::array<double, 3> GradBasis1( const std::array<double, 3> &p ) const;
   std::array<double, 3> GradBasis2( const std::array<double, 3> &p ) const;
   std::array<double, 3> GradBasis3( const std::array<double, 3> &p ) const;
   std::array<double, 3> GradBasis4( const std::array<double, 3> &p ) const;
   std::array<double, 3> GradBasis5( const std::array<double, 3> &p ) const;
   std::array<double, 3> GradBasis6( const std::array<double, 3> &p ) const;
   std::array<double, 3> GradBasis7( const std::array<double, 3> &p ) const;
   std::array<double, 3> GradBasis8( const std::array<double, 3> &p ) const;

   double Basis( int i, const std::array<double, 3> &p ) const;
   std::array<double, 3> GradBasis( int i, const std::array<double, 3> &p ) const;
};
