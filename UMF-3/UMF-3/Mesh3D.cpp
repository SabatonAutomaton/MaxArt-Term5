#include "Mesh3D.h"
#include <fstream>
#include <iostream>

void Mesh3D::InputMeshes( )
{
   std::ifstream f( "mesh3d.txt" );
   if ( !f )
   {
      std::cerr << "File open error: mesh3d.txt\n";
      return;
   }

   double ax, bx, qx;
   double ay, by, qy;
   double az, bz, qz;

   f >> uniformX >> n_x >> ax >> bx;
   f >> uniformY >> n_y >> ay >> by;
   f >> uniformZ >> n_z >> az >> bz;

   meshX.resize( n_x );
   meshY.resize( n_y );
   meshZ.resize( n_z );

   if ( uniformX )
   {
      CalculateUniformMesh( meshX, ax, bx, n_x );
   }
   else
   {
      f >> qx;
      CalculateNonUniformMesh( meshX, ax, bx, n_x, qx );
   }

   if ( uniformY )
   {
      CalculateUniformMesh( meshY, ay, by, n_y );
   }
   else
   {
      f >> qy;
      CalculateNonUniformMesh( meshY, ay, by, n_y, qy );
   }

   if ( uniformZ )
   {
      CalculateUniformMesh( meshZ, az, bz, n_z );
   }
   else
   {
      f >> qz;
      CalculateNonUniformMesh( meshZ, az, bz, n_z, qz );
   }

   BuildStructuredHexMesh( );
}

void Mesh3D::BuildStructuredHexMesh( )
{
   nodes.clear( );
   nodes.reserve( n_x * n_y * n_z );
   for ( int k = 0; k < n_z; ++k )
   {
      for ( int j = 0; j < n_y; ++j )
      {
         for ( int i = 0; i < n_x; ++i )
         {
            nodes.push_back( { meshX[i], meshY[j], meshZ[k] } );
         }
      }
   }

   auto idx = [this]( int i, int j, int k ) -> int
      {
         return i + j * n_x + k * n_x * n_y;
      };

   hexes.clear( );
   for ( int k = 0; k < n_z - 1; ++k )
   {
      for ( int j = 0; j < n_y - 1; ++j )
      {
         for ( int i = 0; i < n_x - 1; ++i )
         {
            const int v000 = idx( i, j, k );
            const int v100 = idx( i + 1, j, k );
            const int v110 = idx( i + 1, j + 1, k );
            const int v010 = idx( i, j + 1, k );
            const int v001 = idx( i, j, k + 1 );
            const int v101 = idx( i + 1, j, k + 1 );
            const int v111 = idx( i + 1, j + 1, k + 1 );
            const int v011 = idx( i, j + 1, k + 1 );
            hexes.push_back( { { v000, v100, v110, v010, v001, v101, v111, v011 } } );
         }
      }
   }
}

