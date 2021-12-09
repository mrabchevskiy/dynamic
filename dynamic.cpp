                                                                                                                              /*
 Copyright Mykola Rabchevskiy 2021.
 Distributed under the Boost Software License, Version 1.0.
 (See http://www.boost.org/LICENSE_1_0.txt)
______________________________________________________________________________

 2021.11.16

 2021.12.08 Assignment text added

 Test application for `Dynamic` module
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

    if( f.defined()    ) correct = false;
    if( f.order() != 4 ) correct = false;
                                                                                                                              /*
    Original function to be approximated is a second order polynomial:
                                                                                                                              */
    auto u = [&]( const Time& t )->Real{ return ( 1.0*t - 2.0 )*t + 3.0; };  // : t^2 - 2t + 3
                                                                                                                              /*
    Set single sample; it should define constant value:
                                                                                                                              */
    {
      printf( "\n\n TEST FOT CONSTANT VALUE DEFINED BY THE SINGLE SAMPLE\n" );
      constexpr Real CONSTANT_VALUE{ 3.14 };
      constexpr Time SOME_TIME     { 2.72 };
      constexpr Time ANOTHER_TIME  { 1.00 };
      f.update( SOME_TIME, CONSTANT_VALUE );
      const auto V{ f( ANOTHER_TIME ) };
      if( fabs( V - CONSTANT_VALUE ) > 1.0e-6 ){
        correct = false;
        printf( "\n Test result: failed; expeced %.3f but got %.3f\n", CONSTANT_VALUE, V );
      } else {
        printf( "\n Test result: correct\n" );
      }
    }

    f.clear();
    if( f.length() != 0 ){ correct = false; printf( "\n\n Clearing test failed" ); }
                                                                                                                            /*
    Add samples:
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
                                                                                                                              /*
    Test the copy costructor and assignment operator:
                                                                                                                              */
    auto g{ f };
    printf( "\n\n Constructed by copy:\n" );
    printf( "\n  %2s %7s %7s %7s %7s \n", "#", "t  ", "orig", "proxy", "err" );
    for( auto k: RANGE{ L } ){
      const Real t    { 0.2*double( int( k )  - 5 ) }; // :time point
      const Real orig { u(t)                        }; // :original    value
      const Real proxy{ g(t)                        }; // :approximated value
      const Real err  { proxy - orig                }; // :error
      Rsq += err*err;
      printf( "\n  %2i %7.2f %7.2f %7.2f %7.2f", k+1, t, orig, proxy, err );
    }

    auto h = f;
    {
      printf( "\n\n TEST: RE-APPROXIMATION WITH ASSIGNMENT BUT WITHOUT RESETTING\n"  );
                                                                                                                            /*
      Add samples:
                                                                                                                              */
      for( auto k: RANGE{ L } ){
        const Real t{ 0.2*double( int( k )  - 5 ) };  // :time point
        h.update( t, u( t ) );                        // :set sample
        if( h.length() != L ) correct = false;
      };
                                                                                                                              /*
      Make approximation:
                                                                                                                              */
      auto[ Nr, Ne, Cn, dt ] = h.process();
      auto[ Q, To, Tt, Tx  ] = f.def();  // :Q is desired appoximation polynomial
      printf( "\n   Number of rotations         %i",            Nr );
      printf( "\n   Number of used eigen values %i",            Ne );
      printf( "\n   Matrix condition number     %.2e",          Cn );
      printf( "\n   Elapsed time                %.2f microsec", dt );
      printf( "\n   Time range                  [ %.2f .. %.2f | .. %.2f ] sec", To, Tt, Tx );
      printf( "\n\n Approximation:\n" );
      printf(   "\n  %2s %7s %7s %7s %7s \n", "#", "t  ", "orig", "proxy", "err" );
      for( auto k: RANGE{ L } ){
        const Real t    { 0.2*double( int( k )  - 5 ) }; // :time point
        const Real orig { u(t)                        }; // :original    value
        const Real proxy{ h(t)                        }; // :approximated value
        const Real err  { proxy - orig                }; // :error
        Rsq += err*err;
        printf( "\n  %2i %7.2f %7.2f %7.2f %7.2f", k+1, t, orig, proxy, err );
      }
    }
    printf( "\n\n TEST FOR ASSIGNMENT " );
    try {
      f = h;
      printf( " [ok]\n" );
    } catch(...){
      printf( " [failed]\n" );
      correct = false;
    }
  }

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
    RangePoint rangePoint{};
    for( auto i: RANGE{ 15 } ){
      const Time t { Time( i )             };
      const Real xi{ x( t )                };
      const Real yi{ y( t )                };
      const Real ri{ sqrt( xi*xi + yi*yi ) };
      const Real Xi{ X( t, &rangePoint )   };
      const Real Yi{ Y( t )                };
      const Real Ri{ hypot( Xi,    Yi    ) };
      const Real d { hypot( Xi-xi, Yi-yi ) };
      if( d > maxDeviation ) maxDeviation = d;
      printf( "\n   %2i %6.2f | %7.2f %7.2f %7.2f | %7.2f %7.2f %7.2f | %7.4f", i+1,t, xi,yi,ri, Xi,Yi,Ri, d );
      switch( rangePoint ){
        case RangePoint::INSIDE   : printf( " inside"    ); break;
        case RangePoint::BACKWARD : printf( " backward"  ); break;
        case RangePoint::FORWARD  : printf( " forward"   ); break;
        default                   : printf( " undefined" ); break;
      }//switch
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
