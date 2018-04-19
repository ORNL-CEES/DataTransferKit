/****************************************************************************
 * Copyright (c) 2012-2018 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 *                                                                          *
 * SPDX-License-Identifier: BSD-3-Clause                                    *
 ****************************************************************************/
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
 * The view is copyable to device and provides a device-accessible element
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

#ifdef KOKKOS_ENABLE_CUDA
        static_assert( std::is_same<typename KokkosViewType::memory_space,
                                    Kokkos::CudaSpace>::value == false,
                       "DTK for CUDA currently does not support CudaSpace "
                       "memory. Please use CudaUVM memory space instead." );
#endif
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

} // namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_VIEW_HPP

//---------------------------------------------------------------------------//
// end DTK_View.hpp
//---------------------------------------------------------------------------//
