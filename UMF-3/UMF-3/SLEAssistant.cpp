#include "SLEAssistant.h"
#include <algorithm>
#include <cmath>

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
      for ( int k = 0; k < lu.matrix.ggl[i].size( ); ++k )
         sum += lu.matrix.ggl[i][k] * res[i - k - 1];
      res[i] -= sum;
   }
}

void SLEAssistant::LUReverse( const LU &lu, const std::vector<double> &v, std::vector<double> &res )
{
   std::vector<double> z = v;
   for ( int i = z.size( ) - 1; i >= 0; --i )
   {
      res[i] = z[i] / lu.matrix.di[i];
      for ( int k = 0; k < lu.matrix.ggu[i].size( ); ++k )
         z[i - k - 1] -= res[i] * lu.matrix.ggu[i][k];
   }
}

void SLEAssistant::DiagonalPreconditionMSG( SparseMatrix &M, std::vector<double> &v, std::vector<double> &res )
{
   for ( int i = 0; i < M.n; ++i )
      res[i] = v[i] / M.di[i];
}

void SLEAssistant::DiagonalPreconditionLOS(SparseMatrix& M, std::vector<double>& v, std::vector<double>& res)
{
    for (int i = 0; i < M.n; ++i)
        res[i] = v[i] / M.di[i];
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

void SLEAssistant::ILUDirect( const ILU &ilu,
   const std::vector<double> &v,
   std::vector<double> &res )
{
   res = v;
   for ( int i = 0; i < ( int )res.size( ); ++i )
   {
      double sum = 0.0;
      for ( int k = ilu.ig[i]; k < ilu.ig[i + 1]; ++k )
         sum += ilu.ggl[k] * res[ilu.jg[k]];
      res[i] -= sum;
      // диагональ L = 1, делить не нужно
   }
}

void SLEAssistant::ILUReverse( const ILU &ilu,
   const std::vector<double> &v,
   std::vector<double> &res )
{
   std::vector<double> z = v;
   for ( int i = ( int )z.size( ) - 1; i >= 0; --i )
   {
      res[i] = z[i] / ilu.di[i];                          // U(i,i)
      for ( int k = ilu.ig[i]; k < ilu.ig[i + 1]; ++k )
         z[ilu.jg[k]] -= res[i] * ilu.ggu[k];            // обратная замена по U
   }
}
