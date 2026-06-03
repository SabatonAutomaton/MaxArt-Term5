#include "Problem.h"
#include "Basis3D.h"
#include "Quadratures.h"
#include "LU.h"
#include "LOS.h"
#include "BCGSTAB.h"
#include "SLE.h"
#include "SLEAssistant.h"
#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <set>
#include <string>
#include <stdexcept>
#include <chrono>

void Problem::Run( int testCase )
{
   if ( testCase == 1 )
   {
      dimension = Dimension::OneD;

      Input( );
      BuildMatrix( );
      BuildB( );
      Solve( );
      PrintResult( );

      return;
   }

   if ( testCase == 4 )
   {
      dimension = Dimension::ThreeD;

      for ( int p = 1; p <= 4; ++p )
      {
         std::cout << "\n====================\n";
         std::cout << "Polynomial degree p = "
            << p
            << "\n";
         std::cout << "====================\n";

         Run3DHarmonicTest( p + 9 );
      }

      return;
   }

   if ( testCase == 5 )
   {
      dimension = Dimension::ThreeD;
      RunCoefficientTests( );
      return;
   }

   if ( testCase >= 2 )
   {
      dimension = Dimension::ThreeD;

      Run3DHarmonicTest( testCase );

      return;
   }

   throw std::runtime_error(
      "Unknown test case."
   );
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
      std::transform( mode.begin( ), mode.end( ), mode.begin( ), []( unsigned char c ) { return static_cast< char >( std::toupper( c ) ); } );
      if ( mode == "LOS" )
      {
         solveMethod = SolveMethod::LOSDiagonal;
      }
      else if ( mode == "BCGSTAB" )
      {
         solveMethod = SolveMethod::BCGSTAB;
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
      matrix.ig[i + 1] = matrix.ig[i] + static_cast< int >( rowProfile[i].size( ) );
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

   auto start = std::chrono::steady_clock::now( );
   if ( solveMethod == SolveMethod::LOSDiagonal )
   {
      SLE sle;
      sle.matrix = matrix;
      sle.f = b;
      sle.x.assign( n, 0.0 );
      sle.maxIter = 10000;
      sle.eps = 1e-12;
      q = LOS::SolutionWithDiagonalConditioning( sle );
      std::cout << "Solver: LOS (diagonal preconditioning)\n";
   }
   else if ( solveMethod == SolveMethod::BCGSTAB )
   {
      SLE sle;
      sle.matrix = matrix;
      sle.f = b;
      sle.x.assign( n, 0.0 );
      sle.maxIter = 10000;
      sle.eps = 1e-12;
      q = BCGSTAB::SolutionWithDiagonalConditioning( sle );
      std::cout << "Solver: BCGSTAB (diagonal preconditioning)\n";
   }
   else
   {
      Matrix denseMatrix = ConvertToMatrix( matrix );
      LU lu( denseMatrix, b );
      lu.calcLU( );
      lu.calcY( );
      lu.calcQ( );
      q = lu.q;
      std::cout << "Solver: LU\n";
   }
   auto end = std::chrono::steady_clock::now( );
   auto duration = std::chrono::duration_cast< std::chrono::milliseconds >( end - start );

   std::cout << "Time: " << duration.count( ) << " ms" << std::endl;
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

void Problem::Run3DHarmonicTest( int testNumber )
{
   Input3D( testNumber );
   BuildMatrix3D( );
   BuildB3D( );
   Solve3D( );
   PrintResult3D( );
}

void Problem::Input3D( int testNumber )
{
   current3DTestNumber = testNumber;
   mesh3d.InputMeshes( );

   std::ifstream solverFile( "solver.txt" );
   std::ifstream coeffsFile( "coefficients.txt" );
   if ( solverFile )
   {
      std::string mode;
      solverFile >> mode;
      std::transform( mode.begin( ), mode.end( ), mode.begin( ), []( unsigned char c ) { return static_cast< char >( std::toupper( c ) ); } );
      if ( mode == "LOS" )
      {
         solveMethod = SolveMethod::LOSDiagonal;
      }
      else if ( mode == "BCGSTAB" )
      {
         solveMethod = SolveMethod::BCGSTAB;
      }
      else
      {
         solveMethod = SolveMethod::LU;
      }
   }

   if ( coeffsFile )
   {
      coeffsFile >> func.lambdaVal >> func.wVal >> func.sigmaVal >> func.xiVal;
   }

   const int nNodes = static_cast< int >( mesh3d.nodes.size( ) );
   n = 2 * nNodes;
   matrix = SparseMatrix( n );
   b.assign( n, 0.0 );
   q3d.assign( n, 0.0 );
   usExact3D.assign( nNodes, 0.0 );
   ucExact3D.assign( nNodes, 0.0 );
   elements3D.resize( mesh3d.hexes.size( ) );
   const int variant = current3DTestNumber;

   for ( int i = 0; i < nNodes; ++i )
   {
      const auto &p = mesh3d.nodes[i];
      usExact3D[i] = func.us3DExact( variant, p.x, p.y, p.z );
      ucExact3D[i] = func.uc3DExact( variant, p.x, p.y, p.z );
   }

   for ( size_t e = 0; e < mesh3d.hexes.size( ); ++e )
   {
      const auto &hex = mesh3d.hexes[e];
      std::array<Element3D::Point3, 8> coords{};
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
      std::array<int, 16> map{};
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
      matrix.ig[i + 1] = matrix.ig[i] + static_cast< int >( rowProfile[i].size( ) );
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
      std::array<int, 16> map{};
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
      std::array<int, 16> map{};
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
   const double xMin = mesh3d.meshX.front( ), xMax = mesh3d.meshX.back( );
   const double yMin = mesh3d.meshY.front( ), yMax = mesh3d.meshY.back( );
   const double zMin = mesh3d.meshZ.front( ), zMax = mesh3d.meshZ.back( );
   const int nNodes = static_cast< int >( mesh3d.nodes.size( ) );

   std::vector<char>   isDir( n, 0 );
   std::vector<double> dirVal( n, 0.0 );

   for ( int node = 0; node < nNodes; ++node )
   {
      const auto &p = mesh3d.nodes[node];
      const bool isBoundary =
         std::abs( p.x - xMin ) < 1e-12 || std::abs( p.x - xMax ) < 1e-12 ||
         std::abs( p.y - yMin ) < 1e-12 || std::abs( p.y - yMax ) < 1e-12 ||
         std::abs( p.z - zMin ) < 1e-12 || std::abs( p.z - zMax ) < 1e-12;
      if ( !isBoundary ) continue;

      isDir[2 * node] = 1;  dirVal[2 * node] = usExact3D[node];
      isDir[2 * node + 1] = 1;  dirVal[2 * node + 1] = ucExact3D[node];
   }

   std::vector<int> oldToNew( n, -1 );
   std::vector<int> newToOld;
   newToOld.reserve( n );
   for ( int i = 0; i < n; ++i )
   {
      if ( !isDir[i] )
      {
         oldToNew[i] = static_cast< int >( newToOld.size( ) );
         newToOld.push_back( i );
      }
   }
   const int nFree = static_cast< int >( newToOld.size( ) );

   std::vector<double> bFree( nFree, 0.0 );
   for ( int k = 0; k < nFree; ++k )
      bFree[k] = b[newToOld[k]];

   for ( int row = 0; row < n; ++row )
   {
      if ( isDir[row] ) continue;
      const int kRow = oldToNew[row];
      for ( int idx = matrix.ig[row]; idx < matrix.ig[row + 1]; ++idx )
      {
         const int col = matrix.jg[idx];
         if ( isDir[col] )
            bFree[kRow] -= matrix.ggl[idx] * dirVal[col];
      }
   }

   for ( int col = 0; col < n; ++col )
   {
      if ( !isDir[col] ) continue;
      for ( int idx = matrix.ig[col]; idx < matrix.ig[col + 1]; ++idx )
      {
         const int row = matrix.jg[idx];
         if ( isDir[row] ) continue;
         bFree[oldToNew[row]] -= matrix.ggu[idx] * dirVal[col];
      }
   }


   std::vector<int> igR( nFree + 1, 0 );
   for ( int k = 0; k < nFree; ++k )
   {
      const int row = newToOld[k];
      int cnt = 0;
      for ( int idx = matrix.ig[row]; idx < matrix.ig[row + 1]; ++idx )
      {
         const int col = matrix.jg[idx];
         if ( !isDir[col] )
            ++cnt;
      }
      igR[k + 1] = igR[k] + cnt;
   }
   const int nnzR = igR[nFree];

   SparseMatrix reduced( nFree );
   reduced.ig = igR;
   reduced.jg.resize( nnzR );
   reduced.ggl.resize( nnzR, 0.0 );
   reduced.ggu.resize( nnzR, 0.0 );
   reduced.di.assign( nFree, 0.0 );

   for ( int k = 0; k < nFree; ++k )
   {
      const int row = newToOld[k];
      int pos = igR[k];
      for ( int idx = matrix.ig[row]; idx < matrix.ig[row + 1]; ++idx )
      {
         const int col = matrix.jg[idx];
         if ( isDir[col] ) continue;
         reduced.jg[pos] = oldToNew[col];
         reduced.ggl[pos] = matrix.ggl[idx];
         reduced.ggu[pos] = matrix.ggu[idx];
         ++pos;
      }
      reduced.di[k] = matrix.di[row];
   }
   Matrix denseMatrix;
   std::vector<double> qFree( nFree, 0.0 );
   if ( solveMethod == SolveMethod::LU )
   {
      denseMatrix = ConvertToMatrix( reduced );
      std::cout << "Converting done!";
   }

   auto start = std::chrono::steady_clock::now( );

   if ( solveMethod == SolveMethod::LOSDiagonal )
   {
      SLE sle;
      sle.matrix = reduced;
      sle.f = bFree;
      sle.x.assign( nFree, 0.0 );
      sle.maxIter = 1000;
      sle.eps = 1e-15;
      qFree = LOS::SolutionWithDiagonalConditioning( sle );
      std::cout << "Solver: LOS (diagonal preconditioning)\n";
   }
   else if ( solveMethod == SolveMethod::BCGSTAB )
   {
      SLE sle;
      sle.matrix = reduced;
      sle.f = bFree;
      sle.x.assign( nFree, 0.0 );
      sle.maxIter = 1000;
      sle.eps = 1e-15;
      qFree = BCGSTAB::SolutionWithDiagonalConditioning( sle );
      std::cout << "Solver: BCGSTAB (diagonal preconditioning)\n";
   }
   else
   {
      LU lu( denseMatrix, bFree );
      lu.calcLU( );
      lu.calcY( );
      lu.calcQ( );
      qFree = lu.q;
      std::cout << "Solver: LU\n";
   }

   auto end = std::chrono::steady_clock::now( );
   std::cout << "Time: "
      << std::chrono::duration_cast< std::chrono::milliseconds >( end - start ).count( )
      << " ms\n";

   q3d.assign( n, 0.0 );
   for ( int k = 0; k < nFree; ++k )
      q3d[newToOld[k]] = qFree[k];
   for ( int i = 0; i < n; ++i )
      if ( isDir[i] )
         q3d[i] = dirVal[i];
}

void Problem::PrintResult3D( )
{
   std::ofstream out( "solution.txt" );
   std::ofstream checkOut( "verification.txt" );
   const double tEval = 1.0;
   out << std::scientific << std::setprecision( 12 );
   out << "# x y z u_s u_c u(t)\n";
   if ( checkOut )
   {
      checkOut << std::scientific << std::setprecision( 12 );
      checkOut << "# x y z us_num us_exact err_us uc_num uc_exact err_uc u_num u_exact err_u\n";
   }

   double maxErrUs = 0.0;
   double maxErrUc = 0.0;
   double maxErrU = 0.0;
   const int nNodes = static_cast< int >( mesh3d.nodes.size( ) );
   const int variant = current3DTestNumber;
   const L2Errors3D l2 = ComputeL2Errors3D( tEval );
   for ( int i = 0; i < nNodes; ++i )
   {
      const auto &p = mesh3d.nodes[i];
      const double usNum = q3d[2 * i];
      const double ucNum = q3d[2 * i + 1];
      const double usEx = usExact3D[i];
      const double ucEx = ucExact3D[i];
      const double errUs = std::abs( usNum - usEx );
      const double errUc = std::abs( ucNum - ucEx );
      maxErrUs = std::max( maxErrUs, errUs );
      maxErrUc = std::max( maxErrUc, errUc );

      const double w = func.w3D( variant, p.x, p.y, p.z );
      const double uNum = usNum * std::sin( w * tEval ) + ucNum * std::cos( w * tEval );
      const double uEx = usEx * std::sin( w * tEval ) + ucEx * std::cos( w * tEval );
      const double errU = std::abs( uNum - uEx );
      maxErrU = std::max( maxErrU, errU );

      out << p.x << " " << p.y << " " << p.z << " " << usNum << " " << ucNum << " " << uNum << "\n";
      if ( checkOut )
      {
         checkOut << p.x << " " << p.y << " " << p.z
            << " " << usNum << " " << usEx << " " << errUs
            << " " << ucNum << " " << ucEx << " " << errUc
            << " " << uNum << " " << uEx << " " << errU << "\n";
      }
   }

   std::cout << "Solution written to solution.txt\n";
   std::cout << "Verification written to verification.txt\n";
   std::cout << "3D harmonic test #" << current3DTestNumber << " (t = " << std::scientific << tEval << ")\n";
   std::cout << "Max |us-us_exact| = " << std::scientific << maxErrUs << "\n";
   std::cout << "Max |uc-uc_exact| = " << std::scientific << maxErrUc << "\n";
   std::cout << "Max |u-u_exact| = " << std::scientific << maxErrU << "\n";
   std::cout << "|us* - us|L2 = " << std::scientific << l2.us << "\n";
   std::cout << "|uc* - uc|L2 = " << std::scientific << l2.uc << "\n";
   std::cout << "|u* - u|L2 = " << std::scientific << l2.u << "\n";
}

Problem::L2Errors3D Problem::ComputeL2Errors3D( double t ) const
{
   double l2ErrUs2 = 0.0;
   double l2ErrUc2 = 0.0;
   double l2ErrU2 = 0.0;
   const int variant = current3DTestNumber;
   Basis3D basis3d;
   Quadratures quad;
   const auto quadratures = quad.gaussHex2( );
   auto det3 = []( const std::array<std::array<double, 3>, 3> &m ) -> double
      {
         return
            m[0][0] * ( m[1][1] * m[2][2] - m[1][2] * m[2][1] ) -
            m[0][1] * ( m[1][0] * m[2][2] - m[1][2] * m[2][0] ) +
            m[0][2] * ( m[1][0] * m[2][1] - m[1][1] * m[2][0] );
      };

   for ( size_t e = 0; e < mesh3d.hexes.size( ); ++e )
   {
      const auto &hex = mesh3d.hexes[e];
      for ( const auto &qv : quadratures )
      {
         const std::array<double, 3> rst = { qv[0], qv[1], qv[2] };
         std::array<double, 8> N{};
         std::array<std::array<double, 3>, 8> gradRef{};
         for ( int a = 0; a < 8; ++a )
         {
            N[a] = basis3d.Basis( a, rst );
            gradRef[a] = basis3d.GradBasis( a, rst );
         }

         std::array<std::array<double, 3>, 3> J{};
         double xq = 0.0, yq = 0.0, zq = 0.0;
         double usNumQ = 0.0, ucNumQ = 0.0;
         for ( int a = 0; a < 8; ++a )
         {
            const auto &p = mesh3d.nodes[hex.v[a]];
            J[0][0] += p.x * gradRef[a][0];
            J[0][1] += p.x * gradRef[a][1];
            J[0][2] += p.x * gradRef[a][2];
            J[1][0] += p.y * gradRef[a][0];
            J[1][1] += p.y * gradRef[a][1];
            J[1][2] += p.y * gradRef[a][2];
            J[2][0] += p.z * gradRef[a][0];
            J[2][1] += p.z * gradRef[a][1];
            J[2][2] += p.z * gradRef[a][2];

            xq += N[a] * p.x;
            yq += N[a] * p.y;
            zq += N[a] * p.z;
            usNumQ += N[a] * q3d[2 * hex.v[a]];
            ucNumQ += N[a] * q3d[2 * hex.v[a] + 1];
         }

         const double wj = std::abs( det3( J ) ) * qv[3];
         const double usExQ = func.us3DExact( variant, xq, yq, zq );
         const double ucExQ = func.uc3DExact( variant, xq, yq, zq );
         const double dUs = usNumQ - usExQ;
         const double dUc = ucNumQ - ucExQ;
         l2ErrUs2 += wj * dUs * dUs;
         l2ErrUc2 += wj * dUc * dUc;

         const double w = func.w3D( variant, xq, yq, zq );
         const double uNumQ = usNumQ * std::sin( w * t ) + ucNumQ * std::cos( w * t );
         const double uExQ = usExQ * std::sin( w * t ) + ucExQ * std::cos( w * t );
         const double dU = uNumQ - uExQ;
         l2ErrU2 += wj * dU * dU;
      }
   }

   const double l2ErrUs = std::sqrt( std::max( 0.0, l2ErrUs2 ) );
   const double l2ErrUc = std::sqrt( std::max( 0.0, l2ErrUc2 ) );
   const double l2ErrU = std::sqrt( std::max( 0.0, l2ErrU2 ) );
   return L2Errors3D{ l2ErrUs, l2ErrUc, l2ErrU };
}

Matrix Problem::ConvertToMatrix( const SparseMatrix &sparse ) const
{
   Matrix dense;
   dense.n = sparse.n;
   dense.di = sparse.di;

   // Вычисляем максимальную ширину ленты и сразу заполняем матрицу
   int maxBandwidth = 0;

   // Первый проход - только для определения ширины
   for ( int i = 0; i < sparse.n; ++i )
   {
      if ( sparse.ig[i + 1] > sparse.ig[i] )
      {
         const int j = sparse.jg[sparse.ig[i]];
         const int bandwidth = i - j;
         maxBandwidth = std::max( maxBandwidth, bandwidth );
      }
   }

   dense.ggl.assign( sparse.n, std::vector<double>( maxBandwidth, 0.0 ) );
   dense.ggu.assign( sparse.n, std::vector<double>( maxBandwidth, 0.0 ) );

   // Второй проход - заполнение
   for ( int i = 0; i < sparse.n; ++i )
   {
      for ( int k = sparse.ig[i]; k < sparse.ig[i + 1]; ++k )
      {
         const int j = sparse.jg[k];
         const int jl = maxBandwidth - ( i - j );

         dense.ggl[i][jl] = sparse.ggl[k];
         dense.ggu[i][jl] = sparse.ggu[k];
      }
   }

   return dense;
}

void Problem::RunCoefficientTests( )
{
   std::ifstream csvFile( "../../Lab_3/testsCoeffs.csv" );
   if ( !csvFile )
   {
      std::cerr << "Cannot open testsCoeffs.csv\n";
      return;
   }

   struct TestRow
   {
      int rowNum;
      double lambda;
      double omega;
      double sigma;
      double chi;
      bool hasLambda, hasOmega, hasSigma, hasChi;
   };

   std::vector<TestRow> tests;
   std::string line;

   // Пропускаем заголовок
   std::getline( csvFile, line );

   double currentLambda = 0, currentOmega = 0, currentSigma = 0, currentChi = 0;
   int rowNum = 1;

   while ( std::getline( csvFile, line ) )
   {
      if ( line.empty( ) ) continue;

      std::istringstream iss( line );
      std::string token;
      std::vector<std::string> tokens;

      while ( std::getline( iss, token, ';' ) )
      {
         tokens.push_back( token );
      }

      if ( tokens.size( ) < 4 ) continue;

      TestRow test;
      test.rowNum = rowNum++;
      test.hasLambda = !tokens[0].empty( );
      test.hasOmega = !tokens[1].empty( );
      test.hasSigma = !tokens[2].empty( );
      test.hasChi = !tokens[3].empty( );

      if ( test.hasLambda )
      {
         currentLambda = std::stod( tokens[0] );
      }
      if ( test.hasOmega )
      {
         currentOmega = std::stod( tokens[1] );
      }
      if ( test.hasSigma )
      {
         currentSigma = std::stod( tokens[2] );
      }
      if ( test.hasChi )
      {
         currentChi = std::stod( tokens[3] );
      }

      test.lambda = currentLambda;
      test.omega = currentOmega;
      test.sigma = currentSigma;
      test.chi = currentChi;

      tests.push_back( test );
   }
   csvFile.close( );

   // Открываем файл для записи результатов
   std::ofstream outCsv( "../../Lab_3/testsCoeffs_results.csv" );
   outCsv << "Row;lambda;omega;sigma;xi;time LOS, ms;k LOS;||u - u*||LOS;time BCGSTAB, ms;k BCGSTAB;||u - u*||BCGSTAB\n";
   outCsv << std::scientific << std::setprecision( 2 );

   for ( const auto &test : tests )
   {
      std::cout << "\n========================================\n";
      std::cout << "Test row " << test.rowNum << ": lambda=" << test.lambda
         << " omega=" << test.omega << " sigma=" << test.sigma << " xi=" << test.chi << "\n";
      std::cout << "========================================\n";

      // Записываем коэффициенты в файл
      std::ofstream coeffFile( "coefficients.txt" );
      coeffFile << std::scientific << std::setprecision( 16 );
      coeffFile << test.lambda << " " << test.omega << " " << test.sigma << " " << test.chi;
      coeffFile.close( );

      // Тест с LOS
      std::ofstream solverFile( "solver.txt" );
      solverFile << "LOS";
      solverFile.close( );

      std::cout << "\n--- Running with LOS solver ---\n";
      Input3D( 2 );
      BuildMatrix3D( );
      BuildB3D( );

      auto startLOS = std::chrono::steady_clock::now( );

      const double xMin = mesh3d.meshX.front( ), xMax = mesh3d.meshX.back( );
      const double yMin = mesh3d.meshY.front( ), yMax = mesh3d.meshY.back( );
      const double zMin = mesh3d.meshZ.front( ), zMax = mesh3d.meshZ.back( );
      const int nNodes = static_cast< int >( mesh3d.nodes.size( ) );

      std::vector<char> isDir( n, 0 );
      std::vector<double> dirVal( n, 0.0 );

      for ( int node = 0; node < nNodes; ++node )
      {
         const auto &p = mesh3d.nodes[node];
         const bool isBoundary =
            std::abs( p.x - xMin ) < 1e-12 || std::abs( p.x - xMax ) < 1e-12 ||
            std::abs( p.y - yMin ) < 1e-12 || std::abs( p.y - yMax ) < 1e-12 ||
            std::abs( p.z - zMin ) < 1e-12 || std::abs( p.z - zMax ) < 1e-12;
         if ( !isBoundary ) continue;

         isDir[2 * node] = 1;  dirVal[2 * node] = usExact3D[node];
         isDir[2 * node + 1] = 1;  dirVal[2 * node + 1] = ucExact3D[node];
      }

      std::vector<int> oldToNew( n, -1 );
      std::vector<int> newToOld;
      newToOld.reserve( n );
      for ( int i = 0; i < n; ++i )
      {
         if ( !isDir[i] )
         {
            oldToNew[i] = static_cast< int >( newToOld.size( ) );
            newToOld.push_back( i );
         }
      }
      const int nFree = static_cast< int >( newToOld.size( ) );

      std::vector<double> bFree( nFree, 0.0 );
      for ( int k = 0; k < nFree; ++k )
         bFree[k] = b[newToOld[k]];

      for ( int row = 0; row < n; ++row )
      {
         if ( isDir[row] ) continue;
         const int kRow = oldToNew[row];
         for ( int idx = matrix.ig[row]; idx < matrix.ig[row + 1]; ++idx )
         {
            const int col = matrix.jg[idx];
            if ( isDir[col] )
               bFree[kRow] -= matrix.ggl[idx] * dirVal[col];
         }
      }

      for ( int col = 0; col < n; ++col )
      {
         if ( !isDir[col] ) continue;
         for ( int idx = matrix.ig[col]; idx < matrix.ig[col + 1]; ++idx )
         {
            const int row = matrix.jg[idx];
            if ( isDir[row] ) continue;
            bFree[oldToNew[row]] -= matrix.ggu[idx] * dirVal[col];
         }
      }

      std::vector<int> igR( nFree + 1, 0 );
      for ( int k = 0; k < nFree; ++k )
      {
         const int row = newToOld[k];
         int cnt = 0;
         for ( int idx = matrix.ig[row]; idx < matrix.ig[row + 1]; ++idx )
         {
            const int col = matrix.jg[idx];
            if ( !isDir[col] )
               ++cnt;
         }
         igR[k + 1] = igR[k] + cnt;
      }
      const int nnzR = igR[nFree];

      SparseMatrix reduced( nFree );
      reduced.ig = igR;
      reduced.jg.resize( nnzR );
      reduced.ggl.resize( nnzR, 0.0 );
      reduced.ggu.resize( nnzR, 0.0 );
      reduced.di.assign( nFree, 0.0 );

      for ( int k = 0; k < nFree; ++k )
      {
         const int row = newToOld[k];
         int pos = igR[k];
         for ( int idx = matrix.ig[row]; idx < matrix.ig[row + 1]; ++idx )
         {
            const int col = matrix.jg[idx];
            if ( isDir[col] ) continue;
            reduced.jg[pos] = oldToNew[col];
            reduced.ggl[pos] = matrix.ggl[idx];
            reduced.ggu[pos] = matrix.ggu[idx];
            ++pos;
         }
         reduced.di[k] = matrix.di[row];
      }

      SLE sle;
      sle.matrix = reduced;
      sle.f = bFree;
      sle.x.assign( nFree, 0.0 );
      sle.maxIter = 1000;
      sle.eps = test.chi;
      std::vector<double> qFree = LOS::SolutionWithDiagonalConditioning( sle );

      auto endLOS = std::chrono::steady_clock::now( );
      long long timeLOS = std::chrono::duration_cast< std::chrono::milliseconds >( endLOS - startLOS ).count( );

      q3d.assign( n, 0.0 );
      for ( int k = 0; k < nFree; ++k )
         q3d[newToOld[k]] = qFree[k];
      for ( int i = 0; i < n; ++i )
         if ( isDir[i] )
            q3d[i] = dirVal[i];

      double maxErrLOS = 0.0;
      for ( int i = 0; i < nNodes; ++i )
      {
         const double usNum = q3d[2 * i];
         const double ucNum = q3d[2 * i + 1];
         const double usEx = usExact3D[i];
         const double ucEx = ucExact3D[i];
         maxErrLOS = std::max( maxErrLOS, std::abs( usNum - usEx ) );
         maxErrLOS = std::max( maxErrLOS, std::abs( ucNum - ucEx ) );
      }

      std::cout << "LOS Time: " << timeLOS << " ms\n";
      std::cout << "LOS Iterations: " << sle.iterCount << "\n";
      std::cout << "LOS Max Error: " << maxErrLOS << "\n";

      // Тест с BCGSTAB
      solverFile.open( "solver.txt" );
      solverFile << "BCGSTAB";
      solverFile.close( );

      std::cout << "\n--- Running with BCGSTAB solver ---\n";
      Input3D( 2 );
      BuildMatrix3D( );
      BuildB3D( );

      auto startBCGSTAB = std::chrono::steady_clock::now( );

      // Повторяем подготовку для BCGSTAB
      isDir.assign( n, 0 );
      dirVal.assign( n, 0.0 );

      for ( int node = 0; node < nNodes; ++node )
      {
         const auto &p = mesh3d.nodes[node];
         const bool isBoundary =
            std::abs( p.x - xMin ) < 1e-12 || std::abs( p.x - xMax ) < 1e-12 ||
            std::abs( p.y - yMin ) < 1e-12 || std::abs( p.y - yMax ) < 1e-12 ||
            std::abs( p.z - zMin ) < 1e-12 || std::abs( p.z - zMax ) < 1e-12;
         if ( !isBoundary ) continue;

         isDir[2 * node] = 1;  dirVal[2 * node] = usExact3D[node];
         isDir[2 * node + 1] = 1;  dirVal[2 * node + 1] = ucExact3D[node];
      }

      oldToNew.assign( n, -1 );
      newToOld.clear( );
      newToOld.reserve( n );
      for ( int i = 0; i < n; ++i )
      {
         if ( !isDir[i] )
         {
            oldToNew[i] = static_cast< int >( newToOld.size( ) );
            newToOld.push_back( i );
         }
      }

      bFree.assign( nFree, 0.0 );
      for ( int k = 0; k < nFree; ++k )
         bFree[k] = b[newToOld[k]];

      for ( int row = 0; row < n; ++row )
      {
         if ( isDir[row] ) continue;
         const int kRow = oldToNew[row];
         for ( int idx = matrix.ig[row]; idx < matrix.ig[row + 1]; ++idx )
         {
            const int col = matrix.jg[idx];
            if ( isDir[col] )
               bFree[kRow] -= matrix.ggl[idx] * dirVal[col];
         }
      }

      for ( int col = 0; col < n; ++col )
      {
         if ( !isDir[col] ) continue;
         for ( int idx = matrix.ig[col]; idx < matrix.ig[col + 1]; ++idx )
         {
            const int row = matrix.jg[idx];
            if ( isDir[row] ) continue;
            bFree[oldToNew[row]] -= matrix.ggu[idx] * dirVal[col];
         }
      }

      igR.assign( nFree + 1, 0 );
      for ( int k = 0; k < nFree; ++k )
      {
         const int row = newToOld[k];
         int cnt = 0;
         for ( int idx = matrix.ig[row]; idx < matrix.ig[row + 1]; ++idx )
         {
            const int col = matrix.jg[idx];
            if ( !isDir[col] )
               ++cnt;
         }
         igR[k + 1] = igR[k] + cnt;
      }

      reduced = SparseMatrix( nFree );
      reduced.ig = igR;
      reduced.jg.resize( nnzR );
      reduced.ggl.resize( nnzR, 0.0 );
      reduced.ggu.resize( nnzR, 0.0 );
      reduced.di.assign( nFree, 0.0 );

      for ( int k = 0; k < nFree; ++k )
      {
         const int row = newToOld[k];
         int pos = igR[k];
         for ( int idx = matrix.ig[row]; idx < matrix.ig[row + 1]; ++idx )
         {
            const int col = matrix.jg[idx];
            if ( isDir[col] ) continue;
            reduced.jg[pos] = oldToNew[col];
            reduced.ggl[pos] = matrix.ggl[idx];
            reduced.ggu[pos] = matrix.ggu[idx];
            ++pos;
         }
         reduced.di[k] = matrix.di[row];
      }

      SLE sleBCG;
      sleBCG.matrix = reduced;
      sleBCG.f = bFree;
      sleBCG.x.assign( nFree, 0.0 );
      sleBCG.maxIter = 10000;
      sleBCG.eps = test.chi;
      qFree = BCGSTAB::SolutionWithDiagonalConditioning( sleBCG );

      auto endBCGSTAB = std::chrono::steady_clock::now( );
      long long timeBCGSTAB = std::chrono::duration_cast< std::chrono::milliseconds >( endBCGSTAB - startBCGSTAB ).count( );

      q3d.assign( n, 0.0 );
      for ( int k = 0; k < nFree; ++k )
         q3d[newToOld[k]] = qFree[k];
      for ( int i = 0; i < n; ++i )
         if ( isDir[i] )
            q3d[i] = dirVal[i];

      double maxErrBCGSTAB = 0.0;
      for ( int i = 0; i < nNodes; ++i )
      {
         const double usNum = q3d[2 * i];
         const double ucNum = q3d[2 * i + 1];
         const double usEx = usExact3D[i];
         const double ucEx = ucExact3D[i];
         maxErrBCGSTAB = std::max( maxErrBCGSTAB, std::abs( usNum - usEx ) );
         maxErrBCGSTAB = std::max( maxErrBCGSTAB, std::abs( ucNum - ucEx ) );
      }

      std::cout << "BCGSTAB Time: " << timeBCGSTAB << " ms\n";
      std::cout << "BCGSTAB Iterations: " << sleBCG.iterCount << "\n";
      std::cout << "BCGSTAB Max Error: " << maxErrBCGSTAB << "\n";

      /* ЗАКОММЕНТИРОВАНО - LU слишком медленный
      // Тест с LU
      solverFile.open("solver.txt");
      solverFile << "LU";
      solverFile.close();

      std::cout << "\n--- Running with LU solver ---\n";
      Input3D(2);
      BuildMatrix3D();
      BuildB3D();

      auto startLU = std::chrono::steady_clock::now();

      // Повторяем подготовку для LU
      isDir.assign(n, 0);
      dirVal.assign(n, 0.0);

      for (int node = 0; node < nNodes; ++node)
      {
         const auto& p = mesh3d.nodes[node];
         const bool isBoundary =
            std::abs(p.x - xMin) < 1e-12 || std::abs(p.x - xMax) < 1e-12 ||
            std::abs(p.y - yMin) < 1e-12 || std::abs(p.y - yMax) < 1e-12 ||
            std::abs(p.z - zMin) < 1e-12 || std::abs(p.z - zMax) < 1e-12;
         if (!isBoundary) continue;

         isDir[2 * node] = 1;  dirVal[2 * node] = usExact3D[node];
         isDir[2 * node + 1] = 1;  dirVal[2 * node + 1] = ucExact3D[node];
      }

      oldToNew.assign(n, -1);
      newToOld.clear();
      newToOld.reserve(n);
      for (int i = 0; i < n; ++i)
      {
         if (!isDir[i])
         {
            oldToNew[i] = static_cast<int>(newToOld.size());
            newToOld.push_back(i);
         }
      }

      bFree.assign(nFree, 0.0);
      for (int k = 0; k < nFree; ++k)
         bFree[k] = b[newToOld[k]];

      for (int row = 0; row < n; ++row)
      {
         if (isDir[row]) continue;
         const int kRow = oldToNew[row];
         for (int idx = matrix.ig[row]; idx < matrix.ig[row + 1]; ++idx)
         {
            const int col = matrix.jg[idx];
            if (isDir[col])
               bFree[kRow] -= matrix.ggl[idx] * dirVal[col];
         }
      }

      for (int col = 0; col < n; ++col)
      {
         if (!isDir[col]) continue;
         for (int idx = matrix.ig[col]; idx < matrix.ig[col + 1]; ++idx)
         {
            const int row = matrix.jg[idx];
            if (isDir[row]) continue;
            bFree[oldToNew[row]] -= matrix.ggu[idx] * dirVal[col];
         }
      }

      igR.assign(nFree + 1, 0);
      for (int k = 0; k < nFree; ++k)
      {
         const int row = newToOld[k];
         int cnt = 0;
         for (int idx = matrix.ig[row]; idx < matrix.ig[row + 1]; ++idx)
         {
            const int col = matrix.jg[idx];
            if (!isDir[col])
               ++cnt;
         }
         igR[k + 1] = igR[k] + cnt;
      }

      reduced = SparseMatrix(nFree);
      reduced.ig = igR;
      reduced.jg.resize(nnzR);
      reduced.ggl.resize(nnzR, 0.0);
      reduced.ggu.resize(nnzR, 0.0);
      reduced.di.assign(nFree, 0.0);

      for (int k = 0; k < nFree; ++k)
      {
         const int row = newToOld[k];
         int pos = igR[k];
         for (int idx = matrix.ig[row]; idx < matrix.ig[row + 1]; ++idx)
         {
            const int col = matrix.jg[idx];
            if (isDir[col]) continue;
            reduced.jg[pos] = oldToNew[col];
            reduced.ggl[pos] = matrix.ggl[idx];
            reduced.ggu[pos] = matrix.ggu[idx];
            ++pos;
         }
         reduced.di[k] = matrix.di[row];
      }

      Matrix denseMatrix = ConvertToMatrix(reduced);
      LU lu(denseMatrix, bFree);
      lu.calcLU();
      lu.calcY();
      lu.calcQ();
      qFree = lu.q;

      auto endLU = std::chrono::steady_clock::now();
      long long timeLU = std::chrono::duration_cast<std::chrono::milliseconds>(endLU - startLU).count();

      q3d.assign(n, 0.0);
      for (int k = 0; k < nFree; ++k)
         q3d[newToOld[k]] = qFree[k];
      for (int i = 0; i < n; ++i)
         if (isDir[i])
            q3d[i] = dirVal[i];

      double maxErrLU = 0.0;
      for (int i = 0; i < nNodes; ++i)
      {
         const double usNum = q3d[2 * i];
         const double ucNum = q3d[2 * i + 1];
         const double usEx = usExact3D[i];
         const double ucEx = ucExact3D[i];
         maxErrLU = std::max(maxErrLU, std::abs(usNum - usEx));
         maxErrLU = std::max(maxErrLU, std::abs(ucNum - ucEx));
      }

      std::cout << "LU Time: " << timeLU << " ms\n";
      std::cout << "LU Max Error: " << maxErrLU << "\n";
      */

      // Записываем результаты (только LOS)
      outCsv << test.rowNum << ";"
             << test.lambda << ";"
             << test.omega << ";"
             << test.sigma << ";"
             << test.chi << ";"
             << timeLOS << ";"
             << sle.iterCount << ";"
             << maxErrLOS << ";"
             << timeBCGSTAB << ";"
             << sleBCG.iterCount << ";"
             << maxErrBCGSTAB << ";" << "\n";
   }

   outCsv.close( );
   std::cout << "\n========================================\n";
   std::cout << "All tests completed! Results saved to testsCoeffs_results.csv\n";
   std::cout << "========================================\n";
}

