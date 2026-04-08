#include "Problem.h"
#include "LU.h"
#include "MSG.h"
#include "SLE.h"
#include "SLEAssistant.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <string>
#include <stdexcept>

void Problem::Run( )
{
   Input( );
   BuildMatrix( );
   BuildB( );
   Solve( );
   PrintResult( );
}

void Problem::Input( )
{
   mesh.InputMeshes( );
   func.SetDomain( mesh.meshX.front( ), mesh.meshX.back( ) );

   std::ifstream solverFile( "solver.txt" );
   if ( solverFile )
   {
      std::string mode;
      solverFile >> mode;
      std::transform( mode.begin( ), mode.end( ), mode.begin( ), []( unsigned char c ) { return static_cast<char>( std::toupper( c ) ); } );
      if ( mode == "MSG" )
      {
         solveMethod = SolveMethod::MSGDiagonal;
      }
      else
      {
         solveMethod = SolveMethod::LU;
      }
   }

   elemCount = mesh.n_x - 1;
   n = 2 * mesh.n_x;

   matrix = SparseMatrix( n );
   b.assign( n, 0.0 );
   q.assign( n, 0.0 );
   elements.resize( elemCount );

   std::vector<double> coords( 2, 0.0 );
   std::vector<double> localFs( 2, 0.0 );
   std::vector<double> localFc( 2, 0.0 );
   for ( int i = 0; i < elemCount; ++i )
   {
      coords = { mesh.meshX[i], mesh.meshX[i + 1] };

      for ( int j = 0; j < 2; ++j )
      {
         localFs[j] = func.fs( mesh.meshX[i + j] );
         localFc[j] = func.fc( mesh.meshX[i + j] );
      }

      elements[i] = Element( coords, localFs, localFc );
   }
}

void Problem::BuildMatrix( )
{
   std::vector<std::set<int>> rowProfile( n );
   for ( int e = 0; e < elemCount; ++e )
   {
      const std::vector<int> map = { 2 * e, 2 * e + 1, 2 * e + 2, 2 * e + 3 };
      for ( int i = 0; i < 4; ++i )
      {
         for ( int j = 0; j < 4; ++j )
         {
            if ( map[j] < map[i] )
            {
               rowProfile[map[i]].insert( map[j] );
            }
         }
      }
   }

   matrix.ig.assign( n + 1, 0 );
   for ( int i = 0; i < n; ++i )
   {
      matrix.ig[i + 1] = matrix.ig[i] + static_cast<int>( rowProfile[i].size( ) );
   }

   const int nnz = matrix.ig[n];
   matrix.jg.assign( nnz, 0 );
   matrix.ggl.assign( nnz, 0.0 );
   matrix.ggu.assign( nnz, 0.0 );
   matrix.di.assign( n, 0.0 );

   for ( int i = 0; i < n; ++i )
   {
      int k = matrix.ig[i];
      for ( const int col : rowProfile[i] )
      {
         matrix.jg[k++] = col;
      }
   }

   auto addToGlobal = [this]( int row, int col, double value )
      {
         if ( row == col )
         {
            matrix.di[row] += value;
            return;
         }

         if ( row > col )
         {
            for ( int k = matrix.ig[row]; k < matrix.ig[row + 1]; ++k )
            {
               if ( matrix.jg[k] == col )
               {
                  matrix.ggl[k] += value;
                  return;
               }
            }
         }
         else
         {
            for ( int k = matrix.ig[col]; k < matrix.ig[col + 1]; ++k )
            {
               if ( matrix.jg[k] == row )
               {
                  matrix.ggu[k] += value;
                  return;
               }
            }
         }

         throw std::runtime_error( "Sparse profile mismatch while assembling global matrix" );
      };

   for ( int e = 0; e < elemCount; ++e )
   {
      Element &element = elements[e];
      element.localMatrix.assign( 4, std::vector<double>( 4, 0.0 ) );
      element.BuildLocalMatrix( );

      const std::vector<int> map = { 2 * e, 2 * e + 1, 2 * e + 2, 2 * e + 3 };
      for ( int i = 0; i < 4; ++i )
      {
         for ( int j = 0; j < 4; ++j )
         {
            addToGlobal( map[i], map[j], element.localMatrix[i][j] );
         }
      }
   }
}

void Problem::BuildB( )
{
   std::fill( b.begin( ), b.end( ), 0.0 );
   for ( int e = 0; e < elemCount; ++e )
   {
      Element &element = elements[e];
      element.localB.assign( 4, 0.0 );
      element.BuildLocalB( );

      const std::vector<int> map = { 2 * e, 2 * e + 1, 2 * e + 2, 2 * e + 3 };
      for ( int i = 0; i < 4; ++i )
      {
         b[map[i]] += element.localB[i];
      }
   }
}

void Problem::ApplyDirichlet( int dof, double value )
{
   for ( int row = dof + 1; row < n; ++row )
   {
      for ( int k = matrix.ig[row]; k < matrix.ig[row + 1]; ++k )
      {
         if ( matrix.jg[k] == dof )
         {
            b[row] -= matrix.ggl[k] * value;
            matrix.ggl[k] = 0.0;
            matrix.ggu[k] = 0.0;
            break;
         }
      }
   }

   for ( int k = matrix.ig[dof]; k < matrix.ig[dof + 1]; ++k )
   {
      const int col = matrix.jg[k];
      b[col] -= matrix.ggu[k] * value;
      matrix.ggu[k] = 0.0;
      matrix.ggl[k] = 0.0;
   }

   matrix.di[dof] = 1.0;
   b[dof] = value;
}

void Problem::Solve( )
{
   const int leftNode = 0;
   const int rightNode = mesh.n_x - 1;
   const double dAs_dz_right = func.thetaS( mesh.meshX[rightNode] );
   const double dAc_dz_right = func.thetaC( mesh.meshX[rightNode] );

   // Условие Неймана на правой границе.
   b[2 * rightNode] += dAs_dz_right;
   b[2 * rightNode + 1] += dAc_dz_right;

   // Условие Дирихле на левой границе.
   ApplyDirichlet( 2 * leftNode, func.boundaryValue( 0, mesh.meshX[leftNode] ) );
   ApplyDirichlet( 2 * leftNode + 1, func.boundaryValue( 1, mesh.meshX[leftNode] ) );

   const auto isSymmetric = [this]( double tol )
      {
         if ( matrix.ggl.size( ) != matrix.ggu.size( ) )
         {
            return false;
         }
         for ( size_t k = 0; k < matrix.ggl.size( ); ++k )
         {
            if ( std::abs( matrix.ggl[k] - matrix.ggu[k] ) > tol )
            {
               return false;
            }
         }
         return true;
      };

   if ( solveMethod == SolveMethod::MSGDiagonal && isSymmetric( 1e-12 ) )
   {
      SLE sle;
      sle.matrix = matrix;
      sle.f = b;
      sle.x.assign( n, 0.0 );
      sle.maxIter = 10000;
      sle.eps = 1e-12;
      q = MSG::SolutionWithDiagonalConditioning( sle );
      std::cout << "Solver: MSG (diagonal preconditioning)\n";
   }
   else
   {
      if ( solveMethod == SolveMethod::MSGDiagonal )
      {
         std::cout << "Solver: MSG requested, but matrix is non-symmetric. Fallback to LU.\n";
      }
      LU lu( matrix );
      std::vector<double> y( n, 0.0 );
      q.assign( n, 0.0 );
      SLEAssistant::LUDirect( lu, b, y );
      SLEAssistant::LUReverse( lu, y, q );
      std::cout << "Solver: LU\n";
   }
}

void Problem::PrintResult( ) const
{
   std::ofstream out( "solution.txt" );
   std::ofstream checkOut( "verification.txt" );
   if ( !out )
   {
      std::cerr << "Cannot open solution.txt for writing\n";
      return;
   }

   out << std::scientific << std::setprecision( 12 );
   out << "# x u_s u_c\n";
   if ( checkOut )
   {
      checkOut << std::scientific << std::setprecision( 12 );
      checkOut << "# x us_num us_exact err_us uc_num uc_exact err_uc\n";
   }

   double maxErrUs = 0.0;
   double maxErrUc = 0.0;
   const double zLeft = mesh.meshX.front( );
   const double zRight = mesh.meshX.back( );
   const double len = ( std::abs( zRight - zLeft ) > 0.0 ) ? ( zRight - zLeft ) : 1.0;

   for ( int node = 0; node < mesh.n_x; ++node )
   {
      const double x = mesh.meshX[node];
      const double usNum = q[2 * node];
      const double ucNum = q[2 * node + 1];
      const double t = ( x - zLeft ) / len;
      const double usRef = 0.2 + 0.6 * t;
      const double ucRef = -0.4 + 0.8 * t;
      const double errUs = std::abs( usNum - usRef );
      const double errUc = std::abs( ucNum - ucRef );

      maxErrUs = std::max( maxErrUs, errUs );
      maxErrUc = std::max( maxErrUc, errUc );

      out << x << " " << usNum << " " << ucNum << "\n";
      if ( checkOut )
      {
         checkOut << x << " " << usNum << " " << usRef << " " << errUs
            << " " << ucNum << " " << ucRef << " " << errUc << "\n";
      }
   }

   std::cout << "Solution written to solution.txt\n";
   if ( checkOut )
   {
      std::cout << "Verification written to verification.txt\n";
   }
   std::cout << "Max |us-us_exact| = " << std::scientific << maxErrUs << "\n";
   std::cout << "Max |uc-uc_exact| = " << std::scientific << maxErrUc << "\n";
   std::cout << "First values:\n";
   for ( int node = 0; node < std::min( mesh.n_x, 5 ); ++node )
   {
      std::cout << "x=" << mesh.meshX[node]
         << " us=" << q[2 * node]
         << " uc=" << q[2 * node + 1] << "\n";
   }
}
