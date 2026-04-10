#include "Problem.h"
#include "LU.h"
#include "MSG.h"
#include "SLE.h"
#include "SLEAssistant.h"
#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <string>
#include <stdexcept>

void Problem::Run( int testCase )
{
   switch ( testCase )
   {
   case 1:
      dimension = Dimension::OneD;
      Input( );
      BuildMatrix( );
      BuildB( );
      Solve( );
      PrintResult( );
      break;
   case 2:
      dimension = Dimension::ThreeD;
      Run3DHarmonicTest( false );
      break;
   case 3:
      dimension = Dimension::ThreeD;
      Run3DHarmonicTest( true );
      break;
   default:
      throw std::runtime_error( "Unknown test case. Use 1 (1D), 2 (3D harmonic linear #1), 3 (3D harmonic linear #2)." );
   }
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
            if ( fixedDofs.empty( ) || !fixedDofs[row] )
            {
               b[row] -= matrix.ggl[k] * value;
            }
            matrix.ggl[k] = 0.0;
            matrix.ggu[k] = 0.0;
            break;
         }
      }
   }

   for ( int k = matrix.ig[dof]; k < matrix.ig[dof + 1]; ++k )
   {
      const int col = matrix.jg[k];
      if ( fixedDofs.empty( ) || !fixedDofs[col] )
      {
         b[col] -= matrix.ggu[k] * value;
      }
      matrix.ggu[k] = 0.0;
      matrix.ggl[k] = 0.0;
   }

   matrix.di[dof] = 1.0;
   b[dof] = value;
   if ( !fixedDofs.empty( ) )
   {
      fixedDofs[dof] = 1;
   }
}

void Problem::Solve( )
{
   fixedDofs.assign( n, 0 );
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

   if ( solveMethod == SolveMethod::MSGDiagonal )
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
}

void Problem::Run3DHarmonicTest( bool linearVariant2 )
{
   Input3D( linearVariant2 );
   BuildMatrix3D( );
   BuildB3D( );
   Solve3D( );
   PrintResult3D( );
}

void Problem::Input3D( bool linearVariant2 )
{
   use3DLinearVariant2 = linearVariant2;
   mesh3d.InputMeshes( );

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

   const int nNodes = static_cast<int>( mesh3d.nodes.size( ) );
   n = 2 * nNodes;
   matrix = SparseMatrix( n );
   b.assign( n, 0.0 );
   q3d.assign( n, 0.0 );
   usExact3D.assign( nNodes, 0.0 );
   ucExact3D.assign( nNodes, 0.0 );
   elements3D.resize( mesh3d.hexes.size( ) );
   const int variant = linearVariant2 ? 1 : 0;

   for ( int i = 0; i < nNodes; ++i )
   {
      const auto &p = mesh3d.nodes[i];
      usExact3D[i] = func.us3DExact( variant, p.x, p.y, p.z );
      ucExact3D[i] = func.uc3DExact( variant, p.x, p.y, p.z );
   }

   for ( size_t e = 0; e < mesh3d.hexes.size( ); ++e )
   {
      const auto &hex = mesh3d.hexes[e];
      std::array<Element3D::Point3, 8> coords {};
      std::vector<double> localFs( 8, 0.0 );
      std::vector<double> localFc( 8, 0.0 );
      std::vector<double> localLambda( 8, 0.0 );
      std::vector<double> localW( 8, 0.0 );
      std::vector<double> localSigma( 8, 0.0 );
      std::vector<double> localXi( 8, 0.0 );

      for ( int i = 0; i < 8; ++i )
      {
         const auto &p = mesh3d.nodes[hex.v[i]];
         coords[i] = { p.x, p.y, p.z };
         localLambda[i] = func.lambda3D( variant, p.x, p.y, p.z );
         localW[i] = func.w3D( variant, p.x, p.y, p.z );
         localSigma[i] = func.sigma3D( variant, p.x, p.y, p.z );
         localXi[i] = func.xi3D( variant, p.x, p.y, p.z );
         localFs[i] = func.fs3D( variant, p.x, p.y, p.z );
         localFc[i] = func.fc3D( variant, p.x, p.y, p.z );
      }

      elements3D[e] = Element3D( coords, localFs, localFc, localLambda, localW, localSigma, localXi );
   }
}

void Problem::BuildMatrix3D( )
{
   std::vector<std::set<int>> rowProfile( n );
   for ( const auto &hex : mesh3d.hexes )
   {
      std::array<int, 16> map {};
      for ( int i = 0; i < 8; ++i )
      {
         map[2 * i] = 2 * hex.v[i];
         map[2 * i + 1] = 2 * hex.v[i] + 1;
      }

      for ( int i = 0; i < 16; ++i )
      {
         for ( int j = 0; j < 16; ++j )
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

         throw std::runtime_error( "Sparse profile mismatch while assembling 3D global matrix" );
      };

   for ( size_t e = 0; e < mesh3d.hexes.size( ); ++e )
   {
      elements3D[e].BuildLocalMatrix( );
      std::array<int, 16> map {};
      for ( int i = 0; i < 8; ++i )
      {
         map[2 * i] = 2 * mesh3d.hexes[e].v[i];
         map[2 * i + 1] = 2 * mesh3d.hexes[e].v[i] + 1;
      }

      for ( int i = 0; i < 16; ++i )
      {
         for ( int j = 0; j < 16; ++j )
         {
            addToGlobal( map[i], map[j], elements3D[e].localMatrix[i][j] );
         }
      }
   }
}

void Problem::BuildB3D( )
{
   std::fill( b.begin( ), b.end( ), 0.0 );
   for ( size_t e = 0; e < mesh3d.hexes.size( ); ++e )
   {
      elements3D[e].BuildLocalB( );
      std::array<int, 16> map {};
      for ( int i = 0; i < 8; ++i )
      {
         map[2 * i] = 2 * mesh3d.hexes[e].v[i];
         map[2 * i + 1] = 2 * mesh3d.hexes[e].v[i] + 1;
      }

      for ( int i = 0; i < 16; ++i )
      {
         b[map[i]] += elements3D[e].localB[i];
      }
   }
}

void Problem::ApplyDirichlet3D( int dof, double value )
{
   ApplyDirichlet( dof, value );
}

void Problem::Solve3D( )
{
   fixedDofs.assign( n, 0 );
   const double xMin = mesh3d.meshX.front( ), xMax = mesh3d.meshX.back( );
   const double yMin = mesh3d.meshY.front( ), yMax = mesh3d.meshY.back( );
   const double zMin = mesh3d.meshZ.front( ), zMax = mesh3d.meshZ.back( );
   const int nNodes = static_cast<int>( mesh3d.nodes.size( ) );

   std::vector<char> isDir( n, 0 );
   std::vector<double> dirVal( n, 0.0 );

   for ( int node = 0; node < nNodes; ++node )
   {
      const auto &p = mesh3d.nodes[node];
      const bool isBoundary =
         std::abs( p.x - xMin ) < 1e-12 || std::abs( p.x - xMax ) < 1e-12 ||
         std::abs( p.y - yMin ) < 1e-12 || std::abs( p.y - yMax ) < 1e-12 ||
         std::abs( p.z - zMin ) < 1e-12 || std::abs( p.z - zMax ) < 1e-12;

      if ( !isBoundary )
      {
         continue;
      }

      isDir[2 * node] = 1;
      dirVal[2 * node] = usExact3D[node];
      isDir[2 * node + 1] = 1;
      dirVal[2 * node + 1] = ucExact3D[node];
   }

   for ( int row = 0; row < n; ++row )
   {
      for ( int k = matrix.ig[row]; k < matrix.ig[row + 1]; ++k )
      {
         const int col = matrix.jg[k];

         if ( isDir[col] && !isDir[row] )
         {
            b[row] -= matrix.ggl[k] * dirVal[col];
         }
         if ( isDir[row] && !isDir[col] )
         {
            b[col] -= matrix.ggu[k] * dirVal[row];
         }

         if ( isDir[row] || isDir[col] )
         {
            matrix.ggl[k] = 0.0;
            matrix.ggu[k] = 0.0;
         }
      }
   }

   for ( int i = 0; i < n; ++i )
   {
      if ( !isDir[i] )
      {
         continue;
      }
      matrix.di[i] = 1.0;
      b[i] = dirVal[i];
   }

   if ( solveMethod == SolveMethod::MSGDiagonal )
   {
      SLE sle;
      sle.matrix = matrix;
      sle.f = b;
      sle.x.assign( n, 0.0 );
      sle.maxIter = 10000;
      sle.eps = 1e-12;
      q3d = MSG::SolutionWithDiagonalConditioning( sle );
      std::cout << "Solver: MSG (diagonal preconditioning)\n";
   }
   else
   {
      LU lu( matrix );
      std::vector<double> y( n, 0.0 );
      q3d.assign( n, 0.0 );
      SLEAssistant::LUDirect( lu, b, y );
      SLEAssistant::LUReverse( lu, y, q3d );
      std::cout << "Solver: LU\n";
   }
}

void Problem::PrintResult3D( ) const
{
   std::ofstream out( "solution.txt" );
   std::ofstream checkOut( "verification.txt" );
   out << std::scientific << std::setprecision( 12 );
   out << "# x y z u_s u_c\n";
   if ( checkOut )
   {
      checkOut << std::scientific << std::setprecision( 12 );
      checkOut << "# x y z us_num us_exact err_us uc_num uc_exact err_uc\n";
   }

   double maxErrUs = 0.0;
   double maxErrUc = 0.0;
   const int nNodes = static_cast<int>( mesh3d.nodes.size( ) );
   for ( int i = 0; i < nNodes; ++i )
   {
      const auto &p = mesh3d.nodes[i];
      const double usNum = q3d[2 * i];
      const double ucNum = q3d[2 * i + 1];
      const double errUs = std::abs( usNum - usExact3D[i] );
      const double errUc = std::abs( ucNum - ucExact3D[i] );
      maxErrUs = std::max( maxErrUs, errUs );
      maxErrUc = std::max( maxErrUc, errUc );

      out << p.x << " " << p.y << " " << p.z << " " << usNum << " " << ucNum << "\n";
      if ( checkOut )
      {
         checkOut << p.x << " " << p.y << " " << p.z
            << " " << usNum << " " << usExact3D[i] << " " << errUs
            << " " << ucNum << " " << ucExact3D[i] << " " << errUc << "\n";
      }
   }

   std::cout << "Solution written to solution.txt\n";
   std::cout << "Verification written to verification.txt\n";
   if ( use3DLinearVariant2 )
   {
      std::cout << "3D harmonic test (linear manufactured #2)\n";
   }
   else
   {
      std::cout << "3D harmonic test (linear manufactured #1)\n";
   }
   std::cout << "Max |us-us_exact| = " << std::scientific << maxErrUs << "\n";
   std::cout << "Max |uc-uc_exact| = " << std::scientific << maxErrUc << "\n";
}
