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
 * \file DTK_FieldTools_def.hpp
 * \author Stuart R. Slattery
 * \brief FieldTools definition.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_FIELDTOOLS_DEF_HPP
#define DTK_FIELDTOOLS_DEF_HPP

#include <algorithm>
#include <iterator>

#include <DTK_Exception.hpp>

#include <Teuchos_Tuple.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_ScalarTraits.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Get a const view of the field. The ArrayRCP object will not manage
 * the memory. 
 */
template<class Field>
Teuchos::ArrayRCP<const typename FieldTools<Field>::value_type>
FieldTools<Field>::view( const Field& field )
{
    size_type field_size = std::distance( FT::begin( field ), 
					  FT::end( field ) );
    return Teuchos::ArrayRCP<const value_type>(
	&*FT::begin( field ), 0, field_size, false );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get a non-const view of the field. The ArrayRCP object will not
 * manage the memory. 
 */
template<class Field>
Teuchos::ArrayRCP<typename FieldTools<Field>::value_type>
FieldTools<Field>::nonConstView( const Field& field )
{
    size_type field_size = std::distance( FT::begin( field ), 
					  FT::end( field ) );
    return Teuchos::ArrayRCP<value_type>(
	(value_type*) &*FT::begin( field ), 0, field_size, false );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the local bounding box for a field of coordinates.
 */
template<class Field>
BoundingBox FieldTools<Field>::coordLocalBoundingBox( const Field& field )
{
    double x_min = -Teuchos::ScalarTraits<double>::rmax();
    double y_min = -Teuchos::ScalarTraits<double>::rmax();
    double z_min = -Teuchos::ScalarTraits<double>::rmax();

    double x_max = Teuchos::ScalarTraits<double>::rmax();
    double y_max = Teuchos::ScalarTraits<double>::rmax();
    double z_max = Teuchos::ScalarTraits<double>::rmax();

    std::size_t dim = FT::dim( field );

    typename FT::size_type num_elements = 
	std::distance( FT::begin( field ), FT::end( field ) ) / dim;


    if ( dim > 0 )
    {
	x_min = *std::min_element( 
	    FT::begin( field ),
	    FT::begin( field ) + num_elements );
	x_max = *std::max_element( 
	    FT::begin( field ),
	    FT::begin( field ) + num_elements );
    }
    if ( dim > 1 )
    {
	y_min = *std::min_element( 
	    FT::begin( field ) + num_elements,
	    FT::begin( field ) + 2*num_elements );
	y_max = *std::max_element( 
	    FT::begin( field ) + num_elements,
	    FT::begin( field ) + 2*num_elements );
    }
    if ( dim > 2 )
    {
	z_min = *std::min_element( 
	    FT::begin( field ) + 2*num_elements,
	    FT::begin( field ) + 3*num_elements );
	z_max = *std::max_element( 
	    FT::begin( field ) + 2*num_elements,
	    FT::begin( field ) + 3*num_elements );
    }
    if ( dim > 3 )
    {
	throw InvariantException( 
	    "Points with greater than 3 dimensions not supported" );
    }

    return BoundingBox( x_min, y_min, z_min, x_max, y_max, z_max );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the global bounding box for a field of coordinates.
 */
template<class Field>
BoundingBox FieldTools<Field>::coordGlobalBoundingBox( const Field& field,
						       const RCP_Comm& comm )
{
    BoundingBox local_box = coordLocalBoundingBox( field );
    Teuchos::Tuple<double,6> local_bounds = local_box.getBounds();

    double global_x_min, global_y_min, global_z_min;
    double global_x_max, global_y_max, global_z_max;

    Teuchos::reduceAll<int,double>( *comm, 
				    Teuchos::REDUCE_MIN,
				    1,
				    &local_bounds[0],
				    &global_x_min );

    Teuchos::reduceAll<int,double>( *comm, 
				    Teuchos::REDUCE_MIN,
				    1,
				    &local_bounds[1],
				    &global_y_min );

    Teuchos::reduceAll<int,double>( *comm, 
				    Teuchos::REDUCE_MIN,
				    1,
				    &local_bounds[2],
				    &global_z_min );

    Teuchos::reduceAll<int,double>( *comm, 
				    Teuchos::REDUCE_MAX,
				    1,
				    &local_bounds[3],
				    &global_x_max );

    Teuchos::reduceAll<int,double>( *comm, 
				    Teuchos::REDUCE_MAX,
				    1,
				    &local_bounds[4],
				    &global_y_max );

    Teuchos::reduceAll<int,double>( *comm, 
				    Teuchos::REDUCE_MAX,
				    1,
				    &local_bounds[5],
				    &global_z_max );

    return BoundingBox( global_x_min, global_y_min, global_z_min,
			global_x_max, global_y_max, global_z_max );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_FIELDTOOLS_DEF_HPP

//---------------------------------------------------------------------------//
// end DTK_FieldTools_def.hpp
//---------------------------------------------------------------------------//

