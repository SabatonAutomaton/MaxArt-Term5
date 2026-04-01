#include "MSGFR.h"
#include "Broyden.h"
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <iomanip>

namespace fs = std::filesystem;

fs::path FindWorkspaceRoot( )
{
   fs::path current = fs::current_path( );
   while ( !current.empty( ) )
   {
      if ( fs::exists( current / "MO2" ) && fs::exists( current / "MO_2_Graphiks" ) )
         return current;

      if ( current == current.root_path( ) )
         break;

      current = current.parent_path( );
   }

   return { };
}

void RunPlotter( int functionType )
{
   fs::path root = FindWorkspaceRoot( );
   if ( root.empty( ) )
   {
      std::cout << "Cannot find workspace root for Python plotter." << std::endl;
      return;
   }

   fs::path script = root / "MO_2_Graphiks" / "MO_2_Graphiks" / "MO_2_Graphiks.py";
   fs::path input = root / "MO2"/ "MO2" / "out2.txt";
   fs::path output = root / "MO_2_Graphiks" / "MO_2_Graphiks" / "result.svg";

   if ( !fs::exists( script ) )
   {
      std::cout << "Python plotter script not found: " << script << std::endl;
      return;
   }

   if ( !fs::exists( input ) )
   {
      std::cout << "Trajectory file not found: " << input << std::endl;
      return;
   }

   std::string args = "\"" + script.string( ) + "\" --input \"" + input.string( ) +
      "\" --output \"" + output.string( ) + "\" --function " + std::to_string( functionType );

   std::string pyCommand = "py " + args;
   int exitCode = std::system( pyCommand.c_str( ) );

   if ( exitCode != 0 )
   {
      std::string pythonCommand = "python " + args;
      exitCode = std::system( pythonCommand.c_str( ) );
   }

   if ( exitCode != 0 )
      std::cout << "Failed to run Python plotter." << std::endl;
}

int main( )
{
   MSGFR msgfr;
   Broyden broyden;
   int functionType = 1; // 1 - quadratic, 2 - Rosenbrock, 3 - target
   int method = 1;
   std::cout << "Function (1 - quadratic, 2 - Rosenbrock, 3 - target): ";
   std::cin >> functionType;
   std::cout << "\nMethod (1 - MSG, 2 - Broyden): ";
   std::cin >> method;
   msgfr.Init( 2 );
   broyden.Init( 2 );
   msgfr.SetFunctionType( functionType );
   broyden.SetFunctionType( functionType );
   double *res;
   switch ( method )
   {
      case 1:
         res = msgfr.Solver( );
         break;
      case 2:
         res = broyden.Solver( );
         break;
      default:
         res = { 0 };
         break;
         
   }
   
   for ( int i = 0; i < 2; i++ )
      std::cout<<std::setprecision(16) << res[i] << std::endl;

   RunPlotter( functionType );
}
