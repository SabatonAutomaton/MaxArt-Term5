#include "LLt.h"
#include <cmath>

LLt::LLt( const SparseMatrix &matrix )
{
   n = matrix.n;
   di = matrix.di;
   ig = matrix.ig;
   jg = matrix.jg;
   ggl = matrix.ggl;

   double sumL = 0.0, sumD = 0.0;
   for ( int k = 0; k < matrix.n; k++ )
   {
      sumD = 0;
      int i0 = ig[k];
      int i = i0;
      for ( ; i0 < ig[k + 1]; i0++ )
      {
         sumL = 0;
         int j0 = i, j1 = i0;
         for ( ; j0 < j1; j0++ )
         {
            int t0 = ig[jg[i0]];
            while ( t0 < ig[jg[i0] + 1] )
            {
               if ( jg[j0] == jg[t0] )
               {
                  sumL += ggl[j0] * ggl[t0];
               }
               t0++;
            }
         }
         ggl[i0] -= sumL;
         ggl[i0] /= di[jg[i0]];
         sumD += ggl[i0] * ggl[i0];
      }
      di[k] -= sumD;
      di[k] = sqrt( di[k] );
   }
}
