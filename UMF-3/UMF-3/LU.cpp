#include "LU.h"
void LU::calcLU( ) {
   auto &ggl = matrix.ggl, &ggu = matrix.ggu;
   auto &di = matrix.di;
   const int ASize = di.size( );
   const int K = ggl.empty( ) ? 0 : static_cast<int>( ggl[0].size( ) );

   for ( int i = 0; i < ASize; i++ ) {
      double sumDi = 0;
      int j = i - K;
      for ( int jl = 0; jl < K; jl++, j++ ) {
         double sum1 = 0, sum2 = 0;
         if ( j < 0 ) continue;
         int mj = i - j;
         for ( int mi = 0; mi < jl; mi++, mj++ ) {
            sum1 += ggu[i][mi] * ggl[j][mj];
            sum2 += ggu[j][mj] * ggl[i][mi];
         }
         ggu[i][jl] -= sum1;
         ggl[i][jl] = ( ggl[i][jl] - sum2 ) / di[j];
         sumDi += ggu[i][jl] * ggl[i][jl];
      }
      di[i] -= sumDi;
   }
}

void LU::calcY( ) {
   auto &ggl = matrix.ggl;
   const int ASize = b.size( );
   const int K = ggl.empty( ) ? 0 : static_cast<int>( ggl[0].size( ) );
   y.resize( ASize );

   for ( int i = 0; i < ASize; i++ ) {
      double sum = 0;
      int j = i - K;
      for ( int jl = 0; jl < K; jl++, j++ ) {
         if ( j < 0 ) continue;
         sum += ggl[i][jl] * y[j];
      }
      y[i] = b[i] - sum;
   }
}

void LU::calcQ( ) {
   auto &ggu = matrix.ggu;
   auto &di = matrix.di;
   const int ASize = y.size( );
   const int K = ggu.empty( ) ? 0 : static_cast<int>( ggu[0].size( ) );

   for ( int i = ASize - 1; i >= 0; i-- ) {
      q[i] = y[i] / di[i];
      int j = i - K;
      for ( int jl = 0; jl < K; jl++, j++ ) {
         if ( j < 0 ) continue;
         y[j] -= ggu[i][jl] * q[i];
      }
   }
   y.clear( );
}