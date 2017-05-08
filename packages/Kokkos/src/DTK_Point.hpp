#ifndef DTK_POINT_HPP
#define DTK_POINT_HPP

#include <Kokkos_Macros.hpp>

namespace DataTransferKit
{
class Point
{
  public:
    KOKKOS_INLINE_FUNCTION
    const double &operator()( unsigned int const i ) const
    {
        return _coords[i];
    }

    KOKKOS_INLINE_FUNCTION
    const double &operator[]( unsigned int i ) const { return _coords[i]; }

    KOKKOS_INLINE_FUNCTION
    double &operator()( unsigned int const i ) { return _coords[i]; }

    KOKKOS_INLINE_FUNCTION
    double &operator[]( unsigned int i ) { return _coords[i]; }

    // This should be private but if we make public we can use the list
    // initializer constructor.
    double _coords[3];
};
}

#endif
