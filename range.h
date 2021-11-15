#ifndef RANGE_H_INCLUDED
#define RANGE_H_INCLUDED

#include <type_traits>

template< typename Elem = unsigned > class RANGE {
  static_assert( std::is_integral< Elem >::value, "Integral type required" );
  struct Iter {
    Elem i;
    constexpr Iter( Elem i ): i{ i }{}
    constexpr Elem const& operator*  (               ) const {      return i;        }
    constexpr bool        operator== ( const Iter& x ) const {      return i == x.i; }
    constexpr bool        operator!= ( const Iter& x ) const {      return i <  x.i; }
    constexpr Iter&       operator++ (               )       { i++; return *this;    }
  };
  Iter       I;
  Iter const N;
public:
  constexpr RANGE(            Elem upto ): I{ 0    }, N{ upto }{ assert( upto >= 0    ); }
  constexpr RANGE( Elem from, Elem upto ): I{ from }, N{ upto }{ assert( upto >= from ); }
  constexpr Iter begin(){ return I; }
  constexpr Iter end  (){ return N; }
};//RANGE

#endif // RANGE_H_INCLUDED
