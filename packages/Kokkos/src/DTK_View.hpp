//---------------------------------------------------------------------------//
/*
  Copyright (c) 2012, Stuart R. Slattery
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

  *: Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

  *: Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

  *: Neither the name of the University of Wisconsin - Madison nor the
  names of its contributors may be used to endorse or promote products
  derived from this software without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
//---------------------------------------------------------------------------//
/*!
 * \file DTK_View.hpp
 * \brief Flat view wrapper for Kokkos::View objects.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_VIEW_HPP
#define DTK_VIEW_HPP

#include <Kokkos_Core.hpp>
#include <Kokkos_DynRankView.hpp>

#include <cassert>
#include <type_traits>

namespace DataTransferKit
{

//---------------------------------------------------------------------------//
/*!
 * \class View
 *
 * \brief Basic unmanaged wrapper around Kokkos::View data designed for
 * handling user input in function callbacks.
 *
 * The View exposes the raw pointer under any Kokkos::View so it may be
 * treated as a flat array and so it can be accessed directly by the
 * application in any manner on and off device. If the Kokkos::View lives on
 * the device, the DTK view will provide a raw pointer to that device data and
 * the length of the view under that pointer.
 *
 * The view is copyable to device and provides a device-accesible element
 * accessor to facilitate user implementation.
 *
 * The view is intended only to be used during the local scope of a
 * function. It does not track the lifetime of the Kokkos::View it is wrapping
 * nor does it manage the memory of that Kokkos:View.
 *
 * In C++ the view data can be accessed directly via the [] operator. The raw
 * pointer to the view is provided for accessing view data in C and Fortran.
 */
template <class SC>
class View
{
  public:
    //@{
    //! Type aliases.
    using Scalar = SC;
    //@}

    // Empty constructor.
    KOKKOS_INLINE_FUNCTION
    View()
        : _size( 0 )
        , _data( nullptr )
    { /* ... */
    }

    // Kokkos::View constructor.
    template <class KokkosViewType>
    KOKKOS_INLINE_FUNCTION
    View( KokkosViewType kokkos_view,
          typename std::enable_if<
              Kokkos::is_view<KokkosViewType>::value ||
                  Kokkos::Experimental::is_dyn_rank_view<KokkosViewType>::value,
              void *>::type = nullptr )
        : _size( kokkos_view.size() )
        , _data( kokkos_view.data() )
    {
        // Make sure the Kokkos view value type and the DTK view scalar type is
        // the same.
        static_assert(
            std::is_same<typename KokkosViewType::value_type, SC>::value,
            "Kokkos View value type and DTK View Scalar type do not match" );

        // Make sure the Kokkos view is left layout. The assures that for
        // multi-dimensional views that the data will be column-major.
        static_assert( std::is_same<typename KokkosViewType::array_layout,
                                    Kokkos::LayoutLeft>::value,
                       "Kokkos View layout must be LayoutLeft" );
    }

    // Get size of the view.
    KOKKOS_INLINE_FUNCTION
    size_t size() const { return _size; }

    // Get the raw pointer to the view data.
    KOKKOS_INLINE_FUNCTION
    Scalar *data() { return _data; }

    // Get the const raw pointer to the view data.
    KOKKOS_INLINE_FUNCTION
    const Scalar *data() const { return _data; }

    // Access an element in the view.
    KOKKOS_FORCEINLINE_FUNCTION
    Scalar &operator[]( const size_t i )
    {
        assert( _data );
        assert( i < _size );
        return _data[i];
    }

    // Access a const element in the view.
    KOKKOS_FORCEINLINE_FUNCTION
    Scalar &operator[]( const size_t i ) const
    {
        assert( _data );
        assert( i < _size );
        return _data[i];
    }

  private:
    // Size of the view.
    // Note: The use of size_t matches the return type of
    // Kokkos::View<>::size().
    size_t _size;

    // View raw pointer.
    Scalar *_data;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_VIEW_HPP

//---------------------------------------------------------------------------//
// end DTK_View.hpp
//---------------------------------------------------------------------------//
