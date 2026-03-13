#include "Problem.h"
#include <iostream>
#include <fstream>

double Problem::lambda( double u )
{
   return 1;
}

double Problem::theta( double u )
{
   return 1;
}

double Problem::beta( double u )
{
   return 1;
}

double Problem::uBeta( double u )
{
   return 1;
}

void Problem::input( )
{
   mesh.InputMeshes( );
   std::ifstream f( "params.txt" );
   if ( !f )
   {
      std::cerr << "File open error" << std::endl;
      return;
   }
   f >> eps >> maxIter >> boundaryConditions[0] >> boundaryConditions[1];

   simpleIteration.al = new double *[2];//?
   simpleIteration.au = new double *[2];//?
   for ( int i = 0; i < 2; i++ )
   {
      simpleIteration.al[i] = new double[mesh.n_x];
      simpleIteration.au[i] = new double[mesh.n_x];
   }
   simpleIteration.di = new double [mesh.n_x];
}