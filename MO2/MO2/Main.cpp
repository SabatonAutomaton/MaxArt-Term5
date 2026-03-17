#include "MSGFR.h"
#include "Broyden.h"
#include <iostream>
#include <iomanip>
int main( )
{
   MSGFR msgfr;
   Broyden broyden;
   msgfr.Init( 2 );
   double *res;
   res = msgfr.Solver( );
   for ( int i = 0; i < 2; i++ )
      std::cout<<std::setprecision(16) << res[i] << std::endl;
}