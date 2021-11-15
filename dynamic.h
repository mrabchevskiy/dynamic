                                                                                                                              /*
 Copyright Mykola Rabchevskiy 2021.
 Distributed under the Boost Software License, Version 1.0.
 (See http://www.boost.org/LICENSE_1_0.txt)
________________________________________________________________________________________________________________________________

  2021.11.11 Initial version
________________________________________________________________________________________________________________________________
                                                                                                                              */
#ifndef DYNAMIC_H_INCLUDED
#define DYNAMIC_H_INCLUDED

#include <cstring> // :memset

#include <atomic>
#include <mutex>
#include <thread>

#include "eigen.h"
#include "range.h"
#include "timer.h"

namespace CoreAGI {

  template< unsigned N, typename Real = double > class Dynamic {

    using Time = double;

    struct Sample {
      Time t;
      Real v;
      Sample( const Time& t, const Real& v ): t{  t },  v{  v  }{}
      Sample(                              ): t{ 0.0 }, v{ 0.0 }{}
      bool operator!= ( const Sample& S ) const { return t != S.t; }
    };

    const unsigned                    CAPACITY; // :queue capacity
    const PolynomialBasis< N, Real >& F;        // :basis
    Sample*                           S;        // :queue of samples
    unsigned                          len;      // :actual number of samples
    Polynomial< N, Real >             P;        // :approximation polynomial
    Time                              To;       // :extrapolation horizon
    Time                              Tt;       // :extrapolation horizon
    Time                              Tx;       // :extrapolation horizon
    Time                              T_;       // :size of the full time range [ To .. Tx ]
    mutable std::mutex                mutexP;   // :protects P, To, Tt, Tx, T_
    mutable std::mutex                mutexQ;   // :protects S, len

  public:

    std::atomic< bool > mutant;

    Dynamic( unsigned capacity, const PolynomialBasis< N, Real >& basis ):
      CAPACITY{ capacity               },
      F       { basis                  },
      S       { new Sample[ CAPACITY ] },
      len     { 0                      },
      P{}, To{}, Tt{}, Tx{}, T_{}, mutexP{}, mutexQ{}, mutant{ false }
    {}

    constexpr unsigned order() const { return N; }

    unsigned length() const {
      const std::lock_guard< std::mutex > lock( mutexQ );
      const unsigned L{ len };
      return L;
    }

    std::tuple< Polynomial< N, Real >, Time, Time, Time > def() const {
      const std::lock_guard< std::mutex > lock( mutexP );
      return std::make_tuple( P, To, Tt, Tx );
    }

    void clear(){
      const std::lock_guard< std::mutex > lock( mutexQ );
      len = 0;
      mutant.store( true );
    }

    unsigned update( const Time& t, const Real& v ){
      unsigned L{ 0 };
      {                                                                                                                       /*
        Lock queue:
                                                                                                                              */
        const std::lock_guard< std::mutex > lock( mutexQ );
                                                                                                                              /*
        Update queue:
                                                                                                                              */
        if( len < CAPACITY ){
          S[ len++ ] = Sample{ t, v };
        } else {
          assert( len == CAPACITY );
          for( auto i: RANGE{ 1u, CAPACITY } ) S[ i-1 ] = S[ i ];  // :shift
          S[ CAPACITY-1 ] = Sample{ t, v };
        }
        L = len;
      }
      mutant.store( true );
      return L;
    }//update

    std::tuple<
      unsigned, // :rotation number
      unsigned, // :number of used eigen values
      Real,     // :matrix condition number
      Time      // :elapsed time, microsec
    > process(){
                                                                                                                              /*
      (Re)Calculate approximation:
                                                                                                                              */
      constexpr Real FACTOR{ 0.5 };

      assert( mutant.load() );
      Time     T[ CAPACITY ];
      Real     Y[ CAPACITY ];
      unsigned L;
      {
                                                                                                                              /*
        Lock samples S and copy data into T and V:
                                                                                                                              */
        const std::lock_guard< std::mutex > lock( mutexQ );
        assert( len > 0 );
        for( auto i: RANGE{ len } ) T[i] = S[i].t, Y[i] = S[i].v;
                                                                                                                              /*
        Remember cureent length, original one can be changed any time:
                                                                                                                              */
        L = len;
      }
                                                                                                                              /*
      Local utility values:
                                                                                                                              */
      const Time& to{ T[   0   ]              };
      const Time& tt{ T[ len-1 ]              };
      const Time  tx{ tt + FACTOR*( tt - to ) };
      const Time  t_{ tx - to                 };
                                                                                                                              /*
      Function for mapping time range [ to..tt ] to the [ -1 .. +1 ] range
                                                                                                                              */
      auto U = [&]( Time t )->Real{ return 2.0*( t - to )/t_ - 1.0; };
                                                                                                                              /*
      Approximation:
                                                                                                                              */
      Polynomial< N > p;  // :desired polynomial
      Time            dt; // :elapsed time
      unsigned        nr; // :number of rotation in the Jacoby
      unsigned        nc; // :actual number of used eigen vectors
      Real            cn; // :condition number

      {

        constexpr Real COND{ 1.0e6 };
                                                                                                                              /*
        Convert time to dimensionless X:[ -1 .. 1 ]:
                                                                                                                              */
        Real X[ 32 ];
        for( auto k: RANGE{ L } ) X[k] = U( T[k] );
                                                                                                                              /*
        Compose problem `AC = B`:
                                                                                                                              */
        CoreAGI::Timer timer;
        CoreAGI::Eigen< N, Real > E;
        Real B[N]; memset( B, 0, N*sizeof( Real ) );
        for( auto k: RANGE{ L } ){
          const Real& Xk{ X[k] };
          for( auto i: RANGE{ N } ) for( auto j: RANGE{ i+1 } ) E.add( i, j, F[i]( Xk )*F[j]( Xk ) );
          for( auto i: RANGE{ N } ) B[i] += F[i]( Xk )*Y[k];
        }//for k
                                                                                                                              /*
        Solve problem:
                                                                                                                              */
        Real C[ N ]; memset( C, 0, N*sizeof( Real ) );
        nc = E.linearSystem( C, B, COND );
        nr = E.rotationNumber();
        dt = timer.elapsed( Timer::MICROSEC );
        cn = E.eigenValue( 0 )/E.eigenValue( nc-1 );
                                                                                                                              /*
        Compose desired polynomial as linear combination of elements of polynomial basis:
                                                                                                                              */
        p = F( C );
      }
                                                                                                                              /*
      Lock and update C[*], To, Tt, Tx:
                                                                                                                              */
      {
        const std::lock_guard< std::mutex > lock( mutexP );
        P  = p;
        To = to;
        Tt = tt;
        Tx = tx;
        T_ = t_;
        mutant.store( false );
      }
      return std::make_tuple( nr, nc, cn, dt );
    }//process

    Real operator() ( const Time& t, int* note = nullptr ){
                                                                                                                              /*
      Calculate approximated/extrapolated value.
      If `note` pointer defined, it value asigned:
        0 when t in range [ To, Tx ]
        1 when t > Tx
       -1 when t < To
                                                                                                                              */
      const std::lock_guard< std::mutex > lock( mutexP );
      if( note ) *note = t > Tx ? 1 : ( t < To ? -1 : 0 );
      return P( 2.0*( t  - To )/T_ - 1.0 ); // :mapping t:[ To, Tx ] => x:[ -1, 1 ]
    }//operator()

  };//class Dynamic

}//CoreAGI

#endif // DYNAMIC_H_INCLUDED
