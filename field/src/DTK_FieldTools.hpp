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

#include "DTK_FieldTraits.hpp"
#include <DTK_BoundingBox.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Comm.hpp>

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
template<class Field>
class FieldTools
{
  public:

    //@{
    //! Typedefs. 
    typedef Field                           field_type;
    typedef FieldTraits<Field>              FT;
    typedef typename FT::value_type         value_type;
    typedef typename FT::size_type          size_type;
    typedef Teuchos::Comm<int>              CommType;
    typedef Teuchos::RCP<const CommType>    RCP_Comm;
    //@}

    //! Constructor.
    FieldTools()
    { /* ... */ }

    //! Destructor.
    ~FieldTools()
    { /* ... */ }

    //@{
    //! View methods.
    // Get a const view of the field. The ArrayRCP object will not manage the
    // memory. 
    static Teuchos::ArrayRCP<const value_type> view( const Field& field );

    // Get a non-const view of the field. The ArrayRCP object will not manage
    // the memory. 
    static Teuchos::ArrayRCP<value_type> nonConstView( const Field& field );
    //@}


    //@{
    //! General global mathematical operations.
    // Fill a field with a scalar.
    static void putScalar( Field& field, const value_type& scalar );

    // Fill a field with a different scalar in each dimension.
    static void putScalar( Field& field, 
			   const Teuchos::ArrayView<value_type>& scalars );

    // Scale a each dimension field by a single value.
    static void scale( Field& field, const value_type& scalar );

    // Scale a field by different value for each dimension.
    static void scale( Field& field, 
		       const Teuchos::ArrayView<value_type>& scalars );

    // Compute the infinity norm for each field dimension.
    static void normInf( const Field& field,
			 const RCP_Comm& comm,
			 Teuchos::Array<value_type>& norms );

    // Compute the L1 norm for each field dimension.
    static void norm1( const Field& field, const RCP_Comm& comm,
		       Teuchos::Array<value_type>& norms );


    // Compute the L2 norm for each field dimension.
    static void norm2( const Field& field, const RCP_Comm& comm, 
		       Teuchos::Array<value_type>& norms );

    // Compute the q-norm for each field dimension.
    static void normQ( const Field& field,  const int& q, const RCP_Comm& comm,
		       Teuchos::Array<value_type>& norms );

    // Compute the average value for each field dimension.
    static void average( const Field& field, const RCP_Comm& comm, 
			 Teuchos::Array<value_type>& averages );

    // Get the global length of the field.
    static size_type globalLength( const Field& field, const RCP_Comm& comm );
    //@}

    //@{
    //! Coordinate field operations.
    // Get the local bounding box for a field of coordinates.
    static BoundingBox coordLocalBoundingBox( const Field& field );

    // Get the global bounding box for a field of coordinates.
    static BoundingBox coordGlobalBoundingBox( const Field& field,
					       const RCP_Comm& comm );
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

