                                                                                                                              /*
 Copyright Mykola Rabchevskiy 2021.
 Distributed under the Boost Software License, Version 1.0.
 (See http://www.boost.org/LICENSE_1_0.txt)
________________________________________________________________________________________________________________________________

  2021.11.16 Initial version

  2021.12.08 Assignment operator fixed
________________________________________________________________________________________________________________________________
                                                                                                                              */#ifndef DYNAMIC_H_INCLUDED
#define DYNAMIC_H_INCLUDED

#include <cstring> // :memset

#include <atomic>
#include <mutex>
#include <thread>

#include "eigen.h"
#include "range.h"
#include "timer.h"

namespace CoreAGI {

  enum class RangePoint: unsigned { UNDEFINED = 0, BACKWARD, INSIDE, FORWARD };

  template< unsigned N, typename Real = double > class Dynamic {

    using Time = double;

    struct Sample {
      Time t;
      Real v;
      Sample( const Time& t, const Real& v ): t{  t },  v{  v  }{}
      Sample(                              ): t{ 0.0 }, v{ 0.0 }{}
      bool operator!= ( const Sample& S ) const { return t != S.t; }
    };

//  const unsigned                    CAPACITY; // :queue capacity
    unsigned                          CAPACITY; // :queue capacity
    const PolynomialBasis< N, Real >& F;        // :basis
    Sample*                           S;        // :queue of samples
    unsigned                          pos;      // :sample incl position
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
      pos     { 0                      },
      len     { 0                      },
      P{}, To{}, Tt{}, Tx{}, T_{}, mutexP{}, mutexQ{}, mutant{ false }
    {
      P.undef(); assert( not defined() );
    }

    Dynamic( const Dynamic& D ):
      CAPACITY{ D.CAPACITY             },
      F       { D.F                    },
      S       { new Sample[ CAPACITY ] },
      pos     { D.pos                  },
      len     { D.len                  },
      P       { D.P                    },
      To      { D.To                   },
      Tt      { D.Tt                   },
      Tx      { D.Tx                   },
      T_      { D.T_                   },
      mutexP{}, mutexQ{}, mutant{}
    {
      for( auto i: RANGE{ len } ) S[i] = D.S[i];
      mutant.store( D.mutant.load() );
    }

    Dynamic& operator= ( const Dynamic& D ){                                                                   // [m] 2021.12.08
                                                                                                                              /*
      Assignment can be done only when both sides use the same functional basis:
                                                                                                                              */
      if( &F != &D.F ) throw std::invalid_argument( "Functional basises must be identical" );
      delete[] S;
      CAPACITY = D.CAPACITY;
      S        = new Sample[ CAPACITY ];
      pos      = D.pos;
      len      = D.len;
      P        = D.P;
      To       = D.To;
      Tt       = D.Tt;
      Tx       = D.Tx;
      T_       = D.T_;
      mutant.store( D.mutant.load() );
      for( auto i: RANGE{ len } ) S[i] = D.S[i];
      return *this;
    }

   ~Dynamic(){
      delete[] S;
    }

    constexpr bool defined() const { return P.defined(); }

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
      pos = 0;
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
//        if( len < CAPACITY ){
//          S[ len++ ] = Sample{ t, v };
//        } else {
//          assert( len == CAPACITY );
//          for( auto i: RANGE{ 1u, CAPACITY } ) S[ i-1 ] = S[ i ];  // :shift
//          S[ CAPACITY-1 ] = Sample{ t, v };
//        }
        if( len < CAPACITY ){
                                                                                                                              /*
          `len` and `pos` values are the same:
                                                                                                                              */
          S[ pos++ ] = Sample{ t, v };
          len++;
        } else {
                                                                                                                              /*
          `len` not changed, `pos` changed cyclically:
                                                                                                                              */
          assert( len == CAPACITY );
          if( ++pos >= CAPACITY ) pos = 0;
          S[ pos ] = Sample{ t, v };
        }
        L = len;
      }
      if( len == 1 ){
                                                                                                                              /*
        Assign approximation polynomial that actually represents constant:
                                                                                                                              */
        {
          const std::lock_guard< std::mutex > lock( mutexP );
          P = v;
        }
        mutant.store( false );
      } else {
        mutant.store( true );
      }
      return L;
    }//update

    std::tuple<
      unsigned, // :rotation number
      unsigned, // :number of used eigen values
      Real,     // :matrix condition number
      Time      // :elapsed time, microsec
    > process(){

      if( not mutant.load() ) return std::make_tuple( 0, 0, 0.0, 0.0 ); // :no changes, nothing to do
                                                                                                                              /*
      (Re)Calculate approximation:
                                                                                                                              */
      constexpr Real FACTOR{ 0.5 };


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

    Real operator() ( const Time& t, RangePoint* note = nullptr ){
                                                                                                                              /*
      Calculate approximated/extrapolated value.
      If `note` pointer defined, it value asigned:
        0 when t in range [ To, Tx ]
        1 when t > Tx
       -1 when t < To
                                                                                                                              */
      const std::lock_guard< std::mutex > lock( mutexP );
      auto value = P( 2.0*( t  - To )/T_ - 1.0 ); // :mapping t:[ To, Tx ] => x:[ -1, 1 ]
      if( note ){
        if( std::isnan( value ) ) *note = RangePoint::UNDEFINED;
        else *note = t > Tx ? RangePoint::FORWARD : ( t < To ? RangePoint::BACKWARD : RangePoint::INSIDE );
      }
      return value;
    }//operator()

  };//class Dynamic

}//CoreAGI

#endif // DYNAMIC_H_INCLUDED
