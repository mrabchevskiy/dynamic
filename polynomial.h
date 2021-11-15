                                                                                                                              /*
 Copyright Mykola Rabchevskiy 2021.
 Distributed under the Boost Software License, Version 1.0.
 (See http://www.boost.org/LICENSE_1_0.txt)
________________________________________________________________________________________________________________________________

  2021.11.11 Initial version
________________________________________________________________________________________________________________________________
                                                                                                                              */
#ifndef APPROXIMATION_H_INCLUDED
#define APPROXIMATION_H_INCLUDED

#include <cassert>
#include <cmath>
#include <functional>
#include <span>
#include <initializer_list>

#include "range.h"

namespace CoreAGI {

  template< unsigned L, typename Real = double > class Polynomial {

    Real C[ L ]; // :polynomial coefficients

  public:

    constexpr Polynomial(): C{}{ for( auto i: RANGE{ L } ) C[i] = 0.0; }

    constexpr unsigned order() const { return L; }

    Polynomial( const Polynomial& P ): C{}{ for( auto i: RANGE{ L } ) C[i] = P.C[i]; }

    explicit constexpr Polynomial( std::initializer_list< Real > coeff ): C{}{
      assert( coeff.size() == L );
      for( unsigned i = 0; const auto& Ci: coeff ) C[ i++ ] = Ci;
    }

    constexpr Polynomial& operator = ( const Polynomial& P ){
      for( auto i: RANGE{ L } ) C[i] = P.C[i];
      return *this;
    }

    constexpr Polynomial operator* ( Real factor ) const {
      Polynomial P;
      for( auto i: RANGE{ L } ) P.C[i] = C[i]*factor;
      return P;
    }

    constexpr Polynomial operator+ ( const Polynomial& Q ) const {
      Polynomial P;
      for( auto i: RANGE{ L } ) P.C[i] = C[i] + Q.C[i];
      return P;
    }

    constexpr Polynomial& operator+= ( const Polynomial& Q ){
      for( auto i: RANGE{ L } ) C[i] += Q.C[i];
      return *this;
    }
                                                                                                                              /*
    Polynomial coefficient by index:
                                                                                                                              */
    constexpr Real operator[]( const unsigned& i ) const { assert( i < L ); return C[i]; }
                                                                                                                              /*
    Polynomial value:
                                                                                                                              */
    constexpr Real operator()( const Real& x ) const {
      Real y{ 0.0 };
      for( auto i: RANGE{ L } ) y = y*x + C[i];
      return y;
    }

  };//Polynomial


  template< unsigned N, typename Real = double > class PolynomialBasis {

    Polynomial< N, Real > f[ N ];

  public:

    constexpr unsigned size() const { return N; }

    explicit constexpr PolynomialBasis( const std::initializer_list< const Polynomial< N, Real > > G ){
      assert( G.size() == N );
      for( unsigned i = 0; const auto& Gi: G ) f[ i++ ] = Gi;
    }

    Polynomial< N, Real > operator()( const std::initializer_list< Real > coeff ) const {
      assert( coeff.size() == N );
      Polynomial< N, Real > P;
      for( unsigned i = 0; const auto& Ci: coeff ) P += f[ i++ ]*Ci;
      return P;
    }

    Polynomial< N, Real > operator()( const Real* coeff ) const {
      Polynomial< N, Real > P;
      for( auto i: RANGE{ N } ) P += f[i]*coeff[i];
      return P;
    }

    const Polynomial< N, Real > operator[] ( unsigned i ) const { return f[i]; }

    constexpr PolynomialBasis& operator= ( const PolynomialBasis< N, Real >& basis ){
      for( auto i: RANGE{ N } ) f[i] = basis.f[i];
      return *this;
    }

  };//PolynomialBasis
                                                                                                                              /*
  Functional basises composed of Chebyshev polynomials:
                                                                                                                              */
  using Real = double;

  constexpr PolynomialBasis< 2, Real > Chebyshev2 {
                                                                                                                              /*
                     x^1  x^0
                     ---  ---                                                                                                         */
    Polynomial< 2 >{ 0.0, 1.0 },
    Polynomial< 2 >{ 1.0, 0.0 }
  };

  constexpr PolynomialBasis< 3, Real > Chebyshev3 {
                                                                                                                              /*
                     x^2  x^1  x^0
                     ---  ---  ---                                                                                                         */
    Polynomial< 3 >{ 0.0, 0.0, 1.0 },
    Polynomial< 3 >{ 0.0, 1.0, 0.0 },
    Polynomial< 3 >{ 2.0, 0.0,-1.0 }
  };

  constexpr PolynomialBasis< 4, Real > Chebyshev4 {
                                                                                                                              /*
                     x^3  x^2  x^1  x^0
                     ---  ---  ---  ---                                                                                                         */
    Polynomial< 4 >{ 0.0, 0.0, 0.0, 1.0 },
    Polynomial< 4 >{ 0.0, 0.0, 1.0, 0.0 },
    Polynomial< 4 >{ 0.0, 2.0, 0.0,-1.0 },
    Polynomial< 4 >{ 4.0, 0.0,-3.0, 0.0 }
  };

  constexpr PolynomialBasis< 5, Real > Chebyshev5 {
                                                                                                                              /*
                     x^4  x^3  x^2  x^1  x^0
                     ---  ---  ---  ---  ---                                                                                                         */
    Polynomial< 5 >{ 0.0, 0.0, 0.0, 0.0, 1.0 },
    Polynomial< 5 >{ 0.0, 0.0, 0.0, 1.0, 0.0 },
    Polynomial< 5 >{ 0.0, 0.0, 2.0, 0.0,-1.0 },
    Polynomial< 5 >{ 0.0, 4.0, 0.0,-3.0, 0.0 },
    Polynomial< 5 >{ 8.0, 4.0,-8.0, 0.0, 1.0 }
  };

  constexpr PolynomialBasis< 6, Real > Chebyshev6 {
                                                                                                                              /*
                      x^5   x^4    x^3   x^2   x^1   x^0
                      ---   ---    ---   ---   ---   ---                                                                                                         */
    Polynomial< 6 >{  0.0,  0.0,   0.0,  0.0,  0.0,  1.0 },
    Polynomial< 6 >{  0.0,  0.0,   0.0,  0.0,  1.0,  0.0 },
    Polynomial< 6 >{  0.0,  0.0,   0.0,  2.0,  0.0, -1.0 },
    Polynomial< 6 >{  0.0,  0.0,   4.0,  0.0, -3.0,  0.0 },
    Polynomial< 6 >{  0.0,  8.0,   4.0, -8.0,  0.0,  1.0 },
    Polynomial< 6 >{ 16.0,  0.0, -20.0,  0.0,  5.0,  1.0 }
  };

}//CoreAGI

#endif
