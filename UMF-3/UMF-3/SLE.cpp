#include "SLE.h"

void SLE::Input( )
{
   InputProperty( );
   InputSLE( );
}

void SLE::InputProperty( )
{
   std::fstream file( "C:\\NumMeth\\data\\property.txt", std::ios::in );
   if ( !file.is_open( ) ) throw std::exception( "The file \"property.txt\" cannot be opened!" );

   file >> matrix.n >> maxIter >> eps;

   file.close( );
}

void SLE::InputSLE( )
{
   std::ifstream file( "C:\\NumMeth\\data\\ig.txt" );
   if ( !file.is_open( ) ) throw std::exception( "The file \"ig.txt\" cannot be opened!" );

   matrix.ig.resize( matrix.n + 1 );
   for ( unsigned int i = 0; i <= matrix.n; ++i )
      file >> matrix.ig[i];

   file.close( );
   //---------------------------------------------------------------------------------
   file.open( "C:\\NumMeth\\data\\jg.txt" );
   if ( !file.is_open( ) ) throw std::exception( "The file \"jg.txt\" cannot be opened!" );

   size_t buf = matrix.ig[matrix.n];
   matrix.jg.resize( buf );
   for ( size_t i = 0; i < buf; ++i )
      file >> matrix.jg[i];

   file.close( );

   if ( matrix.ig[0] )
   {
      for ( int i = 0; i < matrix.ig.size( ); ++i )
      {
         matrix.ig[i]--;
      }
      for ( int j = 0; j < matrix.jg.size( ); ++j )
      {
         matrix.jg[j]--;
      }
   }

   //---------------------------------------------------------------------------------
   matrix.ggl.resize( buf );
   file.open( "C:\\NumMeth\\data\\gg.txt" );
   if ( !file.is_open( ) )throw std::exception( "The file \"gg.txt\" cannot be opened!" );
   for ( size_t i = 0; i < buf; ++i )
      file >> matrix.ggl[i];
   file.close( );
   //---------------------------------------------------------------------------------
   file.open( "C:\\NumMeth\\data\\di.txt" );
   if ( !file.is_open( ) ) throw std::exception( "The file \"di.txt\" cannot be opened!" );

   matrix.di.resize( matrix.n );
   for ( size_t i = 0; i < matrix.n; ++i )
      file >> matrix.di[i];

   file.close( );
   //---------------------------------------------------------------------------------
   file.open( "C:\\NumMeth\\data\\f.txt" );
   if ( !file.is_open( ) ) throw std::exception( "The file \"f.txt\" cannot be opened!" );

   f.resize( matrix.n );
   for ( size_t i = 0; i < matrix.n; ++i )
      file >> f[i];

   file.close( );

}