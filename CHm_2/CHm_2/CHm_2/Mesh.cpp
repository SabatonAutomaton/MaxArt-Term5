#include "Mesh.h"
#include <fstream>
#include <iostream>

void Mesh::InputMeshes( )
{
   double a_x, b_x, q_x, a_t, b_t, q_t;
   std::ifstream f( "mesh.txt" );
   if ( !f )
   {
      std::cerr << "File open error" << std::endl;
      return;
   }

   f >> uniformX >> n_x >> a_x >> b_x;
   meshX = new double[n_x];
   if ( uniformX )
      CalculateUniformMesh( meshX, a_x, b_x, n_x );
   else
   {
      f >> q_x;
      CalculateNonUniformMesh( meshX, a_x, b_x, n_x, q_x );
   }

   f >> uniformT >> n_t >> a_t >> b_t;
   meshT = new double[n_t];
   if ( uniformT )
      CalculateUniformMesh( meshT, a_t, b_t, n_t );
   else
   {
      f >> q_t;
      CalculateNonUniformMesh( meshT, a_t, b_t, n_t, q_t );
   }
}

void Mesh::CalculateUniformMesh( double *mesh, double a, double b, int n )
{
   double h = ( b - a ) / ( n - 1 );
   for ( int i = 0; i < n; i++ )
      mesh[i] = a + i * h;
}

void Mesh::CalculateNonUniformMesh( double *mesh, double a, double b, int n, double q )
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
   mesh[n-1] = b;
}
