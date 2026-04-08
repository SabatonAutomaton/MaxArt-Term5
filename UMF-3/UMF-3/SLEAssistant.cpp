#include "SLEAssistant.h"
#include <algorithm>

double SLEAssistant::DotProduct( const std::vector<double> &v1, const std::vector<double> &v2 )
{
   double sum = 0;
   for ( int i = 0; i < v1.size( ); ++i )
      sum += v1[i] * v2[i];
   return sum;
}

void SLEAssistant::MatrixVectorProduct( const SparseMatrix &m, const std::vector<double> &v, std::vector<double> &res )
{
   std::fill( res.begin( ), res.end( ), 0.0 );
   for ( int i = 0; i < m.n; ++i )
   {
      res[i] = m.di[i] * v[i];
      for ( int k = m.ig[i]; k < m.ig[i + 1]; k++ )
      {
         res[i] += m.ggl[k] * v[m.jg[k]];
         res[m.jg[k]] += m.ggu[k] * v[i];
      }
   }
}

void SLEAssistant::LUDirect( const LU &lu, const std::vector<double> &v, std::vector<double> &res )
{
   double sum = 0;
   res = v;
   for ( int i = 0; i < res.size( ); ++i )
   {
      sum = 0;
      for ( int k = lu.ig[i]; k < lu.ig[i + 1]; ++k )
         sum += lu.ggl[k] * res[lu.jg[k]];
      res[i] -= sum;
   }
}

void SLEAssistant::LUReverse( const LU &lu, const std::vector<double> &v, std::vector<double> &res )
{
   std::vector<double> z = v;
   for ( int i = z.size( ) - 1; i >= 0; --i )
   {
      res[i] = z[i] / lu.di[i];
      for ( int k = lu.ig[i]; k < lu.ig[i + 1]; ++k )
         z[lu.jg[k]] -= res[i] * lu.ggu[k];
   }
}

void SLEAssistant::DiagonalPreconditionMSG( SparseMatrix &M, std::vector<double> &v, std::vector<double> &res )
{
   for ( int i = 0; i < M.n; ++i )
      res[i] = v[i] / M.di[i];
}

void SLEAssistant::DiagonalPreconditionLOS( SparseMatrix &M, std::vector<double> &v, std::vector<double> &res )
{
   for ( int i = 0; i < M.n; ++i )
      res[i] = v[i] / sqrt( M.di[i] );
}

void SLEAssistant::LLtDirect( const LLt &llt, const std::vector<double> &v, std::vector<double> &res )
{
   res = v;
   for ( int i = 0; i < res.size( ); ++i )
   {
      for ( int k = llt.ig[i]; k < llt.ig[i + 1]; ++k )
         res[i] -= llt.ggl[k] * res[llt.jg[k]];
      res[i] /= llt.di[i];
   }
}

void SLEAssistant::LLtReverse( const LLt &llt, const std::vector<double> &v, std::vector<double> &res )
{
   res = v;
   for ( int i = res.size( ) - 1; i >= 0; --i )
   {
      res[i] /= llt.di[i];
      for ( int k = llt.ig[i]; k < llt.ig[i + 1]; ++k )
         res[llt.jg[k]] -= llt.ggl[k] * res[i];
   }
}
