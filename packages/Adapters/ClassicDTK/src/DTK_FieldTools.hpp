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
 * \file DTK_FieldTools.hpp
 * \author Stuart R. Slattery
 * \brief FieldTools definition
 */
//---------------------------------------------------------------------------//

#ifndef DTK_FIELDTOOLS_HPP
#define DTK_FIELDTOOLS_HPP

#include "DTK_BoundingBox.hpp"
#include "DTK_FieldTraits.hpp"

#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_RCP.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \class FieldTools
 * \brief A stateless class for operating on objects that have field traits.
 *
 * Field tools are meant to provide access to typical field/vector
 * operations. In additition, operations intended for coordinate fields are
 * also provided.
 */
//---------------------------------------------------------------------------//
template <class FieldType>
class FieldTools
{
  public:
    //@{
    //! Typedefs.
    typedef FieldType field_type;
    typedef FieldTraits<FieldType> FT;
    typedef typename FT::value_type value_type;
    typedef typename FT::size_type size_type;
    typedef typename FT::iterator iterator;
    typedef typename FT::const_iterator const_iterator;
    typedef Teuchos::Comm<int> CommType;
    typedef Teuchos::RCP<const CommType> RCP_Comm;
    //@}

    //! Constructor.
    FieldTools() { /* ... */}

    //@{
    // Dimension iterators.
    // Get the local size of the dimensions.
    static size_type dimSize( const FieldType &field );

    // Get an iterator to the beginning of a dimension.
    static iterator dimBegin( FieldType &field, const int dim );

    // Get a const iterator to the beginning of a dimension.
    static const_iterator dimBegin( const FieldType &field, const int dim );

    // Get an iterator to the end of a dimension.
    static iterator dimEnd( FieldType &field, const int dim );

    // Get a const iterator to the end of a dimension.
    static const_iterator dimEnd( const FieldType &field, const int dim );
    //@}

    //@{
    // View methods.
    // Get a const view of the field. The ArrayRCP object will not manage the
    // memory.
    static Teuchos::ArrayRCP<const value_type> view( const FieldType &field );

    // Get a non-const view of the field. The ArrayRCP object will not manage
    // the memory.
    static Teuchos::ArrayRCP<value_type> nonConstView( const FieldType &field );

    // Get a deep-copy of the field. The arrayRCP object will manage the
    // memory.
    static Teuchos::ArrayRCP<value_type> copy( const FieldType &field );

    // Get a const view of a dimension of the field. The ArrayRCP object will
    // not manage the memory.
    static Teuchos::ArrayRCP<const value_type> dimView( const FieldType &field,
                                                        const int dim );

    // Get a non-const view of a dimension of the field. The ArrayRCP object
    // will not manage the memory.
    static Teuchos::ArrayRCP<value_type>
    dimNonConstView( const FieldType &field, const int dim );
    //@}

    //@{
    // General global mathematical operations.
    // Fill a field with a scalar.
    static void putScalar( FieldType &field, const value_type &scalar );

    // Fill a field with a different scalar in each dimension.
    static void putScalar( FieldType &field,
                           const Teuchos::ArrayView<value_type> &scalars );

    // Scale a each dimension field by a single value.
    static void scale( FieldType &field, const value_type &scalar );

    // Scale a field by different value for each dimension.
    static void scale( FieldType &field,
                       const Teuchos::ArrayView<value_type> &scalars );

    // Compute the global infinity norm for each field dimension.
    static void normInf( const FieldType &field, const RCP_Comm &comm,
                         Teuchos::Array<value_type> &norms );

    // Compute the global L1 norm for each field dimension.
    static void norm1( const FieldType &field, const RCP_Comm &comm,
                       Teuchos::Array<value_type> &norms );

    // Compute the global L2 norm for each field dimension.
    static void norm2( const FieldType &field, const RCP_Comm &comm,
                       Teuchos::Array<value_type> &norms );

    // Compute the global q-norm for each field dimension.
    static void normQ( const FieldType &field, const RCP_Comm &comm,
                       const int q, Teuchos::Array<value_type> &norms );

    // Compute the global average value for each field dimension.
    static void average( const FieldType &field, const RCP_Comm &comm,
                         Teuchos::Array<value_type> &averages );

    // Get the global size of the field.
    static size_type globalSize( const FieldType &field, const RCP_Comm &comm );
    //@}

    //@{
    // Coordinate field operations.
    // Get the local bounding box for a field of coordinates.
    static BoundingBox coordLocalBoundingBox( const FieldType &field );

    // Get the global bounding box for a field of coordinates.
    static BoundingBox coordGlobalBoundingBox( const FieldType &field,
                                               const RCP_Comm &comm );
    //@}
};

} // end namepsace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_FieldTools_def.hpp"

#endif // end DTK_FIELDTOOLS_HPP

//---------------------------------------------------------------------------//
// end DTK_FieldTools.hpp
//---------------------------------------------------------------------------//
