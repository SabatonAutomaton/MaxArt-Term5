#include "Mesh.h"
#include <fstream>
#include <iostream>
#include <cmath>

void Mesh::InputMeshes( )
{
   double a_x, b_x, q_x;
   std::ifstream f( "mesh.txt" );
   if ( !f )
   {
      std::cerr << "File open error" << std::endl;
      return;
   }

   f >> uniformX >> n_x >> a_x >> b_x;
   meshX.resize( n_x );
   if ( uniformX )
      CalculateUniformMesh( meshX, a_x, b_x, n_x );
   else
   {
      f >> q_x;
      CalculateNonUniformMesh( meshX, a_x, b_x, n_x, q_x );
   }
}

void Mesh::CalculateUniformMesh( std::vector<double> &mesh, double a, double b, int n )
{
   double h = ( b - a ) / ( n - 1 );
   for ( int i = 0; i < n; i++ )
      mesh[i] = a + i * h;
}

void Mesh::CalculateNonUniformMesh( std::vector<double> &mesh, double a, double b, int n, double q )
{
   double h = ( b - a ) * ( 1 - q ) / ( 1 - pow( q, n - 1 ) );
   double hi = h;
   mesh[0] = a;
   double sum = 0;
   for ( int i = 1; i < n - 1; i++ )
   {
      sum += h * pow( q, i - 1 );
      mesh[i] = a + sum;
   }
   mesh[n - 1] = b;
}
