#include "SimpleIteration.h"
#include <iostream>
#include <fstream>

double SimpleIteration::Condition( int number, double x )
{
   return functions.boundaryValue( number, x );
}

void SimpleIteration::Condition1( )
{
   if ( elements[0].cond == 1 )
   {
      matrix.di[0] = 1;
      b[0] = Condition( 1, 0 );;
      for ( int i = 0; i < 2; i++ )
      {
         matrix.ggu[1][1] = 0;
         matrix.ggu[2][0] = 0;
      }
   }
   if ( elements[elemCount - 1].cond == 1 )
   {
      matrix.di[n - 1] = 1;
      b[n - 1] = Condition( 1, mesh.meshX[elemCount] );;
      for ( int i = 0; i < 2; i++ )
      {
         matrix.ggl[n - 1][0] = 0;
         matrix.ggl[n - 1][1] = 0;
      }
   }
}

void SimpleIteration::Input( )
{
   mesh.InputMeshes( );
   std::ifstream file( "params.txt" );
   if ( !file )
   {
      std::cerr << "File open error" << std::endl;
      return;
   }
   file >> eps >> maxIter >> boundaryConditions[0] >> boundaryConditions[1];

   elemCount = mesh.n_x - 1;
   n = 3 * elemCount - elemCount + 1;
   matrix.n = n;
   matrix.ggl.resize( n );
   matrix.ggu.resize( n );
   for ( int i = 0; i < n; i++ )
   {
      matrix.ggl[i] = std::vector<double>( 2, 0.0 );
      matrix.ggu[i] = std::vector<double>( 2, 0.0 );
   }
   matrix.di = std::vector<double>( n, 0.0 );

   b.resize( n );
   lambdaNodes.resize( n );
   for ( int i = 0; i <= elemCount; i++ )
   {
      b[i * 2] = functions.f( mesh.meshX[i] );
      lambdaNodes[i * 2] = functions.lambda( mesh.meshX[i] );
      if ( i < elemCount )
      {
         b[i * 2 + 1] = functions.f( ( mesh.meshX[i] + mesh.meshX[i + 1] ) / 2.0 );
         lambdaNodes[i * 2 + 1] = functions.lambda( ( mesh.meshX[i] + mesh.meshX[i + 1] ) / 2.0 );
      }

   }
   b0 = b;


   std::vector<double> coords;
   elements.resize( elemCount );
   q = std::vector<double>( n, 0.0 );
   q[0] = Condition( 1, 0 );
   q[n - 1] = Condition( 1, mesh.meshX[elemCount] );
   std::vector<double> localB( 3, 0.0 );
   std::vector<double> localLambdaNodes( 3, 0.0 );
   for ( int i = 0; i < elemCount; i++ )
   {
      coords = { mesh.meshX[i], mesh.meshX[i + 1] };
      for ( int j = 0; j < 3; j++ )
      {
         if ( i == 0 )
         {
            localB[j] = b[3 * i + j];
            localLambdaNodes[j] = lambdaNodes[3 * i + j];
         }
         else
         {
            localB[j] = b[2 * i + j];
            localLambdaNodes[j] = lambdaNodes[2 * i + j];
         }
      }
      elements[i] = Element( coords, 0, q, localB, localLambdaNodes, &functions );
   }
   elements[0].cond = boundaryConditions[0];
   elements[elemCount - 1].cond = boundaryConditions[1];


}

void SimpleIteration::BuildMatrix( )
{
   matrix.di = std::vector<double>( n, 0.0 );
   for ( int i = 0; i < n; i++ )
   {
      matrix.ggl[i] = std::vector<double>( 2, 0.0 );
      matrix.ggu[i] = std::vector<double>( 2, 0.0 );
   }

   for ( int i = 0; i < elemCount; i++ )
   {
      elements[i].BuildLocalMatrix( );
      for ( int j = 0; j < 3; j++ )
      {
         if ( i == 0 )
            matrix.di[3 * i + j] += elements[i].localMatrix[j][j];
         else
            matrix.di[2 * i + j] += elements[i].localMatrix[j][j];
      }
      matrix.ggl[2 * i + 1][1] += elements[i].localMatrix[1][0];
      matrix.ggl[2 * i + 2][0] += elements[i].localMatrix[2][0];
      matrix.ggl[2 * i + 2][1] += elements[i].localMatrix[2][1];

      matrix.ggu[2 * i + 1][1] += elements[i].localMatrix[0][1];
      matrix.ggu[2 * i + 2][0] += elements[i].localMatrix[0][2];
      matrix.ggu[2 * i + 2][1] += elements[i].localMatrix[1][2];
   }
}

void SimpleIteration::BuildB( )
{
   b = std::vector<double>( n, 0.0 );
   for ( int i = 0; i < elemCount; i++ )
   {
      elements[i].BuildLocalB( );
      for ( int j = 0; j < 3; j++ )
      {
         if ( i == 0 )
            b[3 * i + j] += elements[i].localB[j];
         else
            b[2 * i + j] += elements[i].localB[j];
      }
   }
}

void SimpleIteration::SolveIter( )
{
   qOLd = q;
   BuildMatrix( );
   BuildB( );
   Condition1( );
   lu = LU( matrix, b );
   lu.calcLU( );
   lu.calcY( );
   lu.calcQ( );
   q = ops.AddVec( ops.MultVecScal( lu.q, w ), ops.MultVecScal( qOLd, 1 - w ) );
   std::vector<double> localQ( 3, 0.0 );
   for ( int i = 0; i < elemCount; i++ )
   {
      for ( int j = 0; j < 3; j++ )
      {
         if ( i == 0 )
         {
            localQ[j] = q[3 * i + j];
         }
         else
         {
            localQ[j] = q[2 * i + j];
         }
      }
      elements[i].q = localQ;
   }

}

void SimpleIteration::Solve( )
{
   int k = 0;
   SolveIter( );
   double residual = CalcResidual( );

   while ( residual > eps )
   {

      k++;
      residual = CalcResidual( );
      std::cout << k << "\t";
      for ( int i = 0; i < n; i++ )
      {
         std::cout << q[i] << "\t";
      }
      std::cout << "\tresidual: " << residual << std::endl;
      for ( int i = 0; i < elemCount; i++ )
         for ( int j = 0; j < 3; j++ )
         {
            if ( i == 0 )
               elements[i].q[j] = q[3 * i + j];
            else
               elements[i].q[j] = q[2 * i + j];
         }
      SolveIter( );
   }

}

double SimpleIteration::CalcResidual( )
{
   return ( ops.DotProduct( ops.AddVec( ops.MultMatVec( matrix, q ), ops.MultVecScal( b, -1 ) ), ops.AddVec( ops.MultMatVec( matrix, q ), ops.MultVecScal( b, -1 ) ) ) ) / ops.DotProduct( b, b );
}

