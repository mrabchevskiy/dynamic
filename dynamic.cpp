                                                                                                                              /*
 Copyright Mykola Rabchevskiy 2021.
 Distributed under the Boost Software License, Version 1.0.
 (See http://www.boost.org/LICENSE_1_0.txt)
 ______________________________________________________________________________

 Test aplication for `Dynamic` module
________________________________________________________________________________________________________________________________
                                                                                                                              */
#include <cstdio>
#include <cmath>

#include "polynomial.h"
#include "dynamic.h"

using namespace CoreAGI;

int main(){

  using Time = double;
  using Real = double;

  bool correct{ true };

  {
                                                                                                                              /*
    Test for polynomial reconstruction:
                                                                                                                              */
    constexpr unsigned L{ 11 };
                                                                                                                              /*
    Desired Dynamic value:
                                                                                                                              */
    auto f = Dynamic( L, Chebyshev4 );
    if( f.order() != 4 ) correct = false;
                                                                                                                              /*
    Original function to be approximated is a second order polynomial:
                                                                                                                              */
    auto u = [&]( const Time& t )->Real{ return ( 1.0*t - 2.0 )*t + 3.0; };  // : t^2 - 2t + 3
                                                                                                                              /*
    Make and add samples:
                                                                                                                              */
    for( auto k: RANGE{ L } ){
      const Real t{ 0.2*double( int( k )  - 5 ) };  // :time point
      f.update( t, u( t ) );                        // :set sample
    };
    if( f.length() != L ) correct = false;
                                                                                                                              /*
    Make approximation:
                                                                                                                              */
    auto[ Nr, Ne, Cn, dt ] = f.process();
                                                                                                                              /*
    Check:
                                                                                                                              */
    auto[ Q, To, Tt, Tx  ] = f.def();  // :Q is desired appoximation polynomial
    printf( "\n\n TEST: APPROXIMATION OF A POLYNOMIAL FUNCTION\n"  );
    printf( "\n   Number of rotations         %i",            Nr );
    printf( "\n   Number of used eigen values %i",            Ne );
    printf( "\n   Matrix condition number     %.2e",          Cn );
    printf( "\n   Elapsed time                %.2f microsec", dt );
    printf( "\n   Time range                  [ %.2f .. %.2f | .. %.2f ] sec", To, Tt, Tx );
    printf( "\n\n Approximation:\n" );
    printf( "\n  %2s %7s %7s %7s %7s \n", "#", "t  ", "orig", "proxy", "err" );
    Real Rsq{ 0.0 };
    for( auto k: RANGE{ L } ){
      const Real t    { 0.2*double( int( k )  - 5 ) }; // :time point
      const Real orig { u(t)                        }; // :original    value
      const Real proxy{ f(t)                        }; // :approximated value
      const Real err  { proxy - orig                }; // :error
      Rsq += err*err;
      printf( "\n  %2i %7.2f %7.2f %7.2f %7.2f", k+1, t, orig, proxy, err );
    }
    Rsq = sqrt( Rsq/Real( L ) );
    constexpr Real EPS{ 1.0e-6 };
    if( Rsq > EPS ) correct = false;
    printf( "\n\n Rsq %.3e  %s", Rsq, Rsq > EPS ? "unacceptable" : "acceptable" );
    printf( "\n\n Test result: %s\n", correct ? "CORRECT" : "FAILURE" );
  };

  {
    printf( "\n\n TEST: APPROXIMATION & EXTRAPOLATION OF THE POINT COORDINATES" );
                                                                                                                              /*
    Test for approximation of the arc trajectory
                                                                                                                              */
    constexpr unsigned L{ 11 };

    auto radians = []( const Real& degrees )->Real{
      constexpr Real FACTOR{ M_PI/180.0 };
      return FACTOR*degrees;
    };

    constexpr Real r{ 10.0           };  // :trajectory radius, m
    constexpr Real w{ radians( 9.0 ) };  // :angular velocity, rad/sec

    auto x = [&]( const Time& t )->Real{ return r*cos( w*t ); };
    auto y = [&]( const Time& t )->Real{ return r*sin( w*t ); };

    auto X{ Dynamic( L, Chebyshev6 ) };
    auto Y{ Dynamic( L, Chebyshev6 ) };

    for( auto i: RANGE{ L } ){
      const Time t{ Time( i ) };
      X.update( t, x( t ) );
      Y.update( t, y( t ) );
    }
   {
      auto[ Nr, Ne, Cn, dt ] = X.process();
      auto[ Q, To, Tt, Tx  ] = X.def();
      printf( "\n\n X approximation:\n"  );
      printf( "\n   Number of rotations         %i",            Nr );
      printf( "\n   Number of used eigen values %i",            Ne );
      printf( "\n   Matrix condition number     %.2e",          Cn );
      printf( "\n   Elapsed time                %.2f microsec", dt );
      printf( "\n   Time range                  [ %.2f .. %.2f | .. %.2f ] sec", To, Tt, Tx );
    }
    {
      auto[ Nr, Ne, Cn, dt ] = Y.process();
      auto[ Q, To, Tt, Tx  ] = Y.def();
      printf( "\n\n Y approximation:\n"  );
      printf( "\n   Number of rotations         %i",            Nr );
      printf( "\n   Number of used eigen values %i",            Ne );
      printf( "\n   Matrix condition number     %.2e",          Cn );
      printf( "\n   Elapsed time                %.2f microsec", dt );
      printf( "\n   Time range                  [ %.2f .. %.2f | .. %.2f ] sec", To, Tt, Tx );
    }
    printf( "\n\n POINT COORDINATES APPROXIMATION & EXTRAPOLATION:\n"  );
    printf( "\n   %2s %6s | %7s %7s %7s | %7s %7s %7s | %7s", "#", "time", "x  ", "y  ","r  ", "x  ", "y  ", "r  ", "dev  " );
    Real maxDeviation{ 0.0 };
    for( auto i: RANGE{ 15 } ){
      const Time t { Time( i )             };
      const Real xi{ x( t )                };
      const Real yi{ y( t )                };
      const Real ri{ sqrt( xi*xi + yi*yi ) };
      const Real Xi{ X( t )                };
      const Real Yi{ Y( t )                };
      const Real Ri{ hypot( Xi,    Yi    ) };
      const Real d { hypot( Xi-xi, Yi-yi ) };
      if( d > maxDeviation ) maxDeviation = d;
      printf( "\n   %2i %6.2f | %7.2f %7.2f %7.2f | %7.2f %7.2f %7.2f | %7.4f", i+1,t, xi,yi,ri, Xi,Yi,Ri, d );
      if( i > 10 ) printf( " extrapolated" );
    }
    const Real trajectoryLength    { r*radians( 120.0 )     };
    const Real acceptableDeviation { trajectoryLength/100.0 };
    printf( "\n\n   Trajectory length       %7.3f m", trajectoryLength    );
    printf(   "\n   Acceptable 1%% deviation %7.3f m", acceptableDeviation );
    printf(   "\n   Max deviation           %7.3f m", maxDeviation        );
    if( maxDeviation <= acceptableDeviation ){
      printf( "\n\n Test result: CORRECT\n" );
    } else {
      correct = false;
      printf( "\n\n Test result: FAILURE\n" );
    }
  }

  printf( "\n Verdict: %s\n", correct ? "CORRECT" : "FAILURE" );

	return correct ? EXIT_SUCCESS : EXIT_FAILURE;
}
