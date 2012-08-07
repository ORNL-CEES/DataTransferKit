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

#include "DTK_Assertion.hpp"

#include <Teuchos_Tuple.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_Ptr.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Get the local size of the dimensions.
template<class Field>
typename FieldTools<Field>::size_type
FieldTools<Field>::dimSize( const Field& field )
{
    int field_dim = FT::dim( field );
    testPrecondition( field_dim > 0 );
    size_type field_size = FT::size( field );
    return field_size / field_dim;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get an iterator to the beginning of a dimension.
 */
template<class Field>
typename FieldTools<Field>::iterator 
FieldTools<Field>::dimBegin( Field& field, const int dim )
{
    int field_dim = FT::dim( field );
    testPrecondition( field_dim > 0 );
    testPrecondition( dim >= 0 && dim < field_dim );
    return FT::begin(field) + dim*dimSize(field);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get a const iterator to the beginning of a dimension.
 */
template<class Field>
typename FieldTools<Field>::const_iterator 
FieldTools<Field>::dimBegin( const Field& field, const int dim )
{
    int field_dim = FT::dim( field );
    testPrecondition( field_dim > 0 );
    testPrecondition( dim >= 0 && dim < field_dim );
    return FT::begin(field) + dim*dimSize(field);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get an iterator to the end of a dimension.
 */
template<class Field>
typename FieldTools<Field>::iterator 
FieldTools<Field>::dimEnd( Field& field, const int dim )
{
    int field_dim = FT::dim( field );
    testPrecondition( field_dim > 0 );
    testPrecondition( dim >= 0 && dim < field_dim );
    return FT::begin(field) + (dim+1)*dimSize(field);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get a const iterator to the end of a dimension.
 */
template<class Field>
typename FieldTools<Field>::const_iterator 
FieldTools<Field>::dimEnd( const Field& field, const int dim )
{
    int field_dim = FT::dim( field );
    testPrecondition( field_dim > 0 );
    testPrecondition( dim >= 0 && dim < field_dim );
    return FT::begin(field) + (dim+1)*dimSize(field);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get a const view of the field. The ArrayRCP object will not manage
 * the memory. 
 */
template<class Field>
Teuchos::ArrayRCP<const typename FieldTools<Field>::value_type>
FieldTools<Field>::view( const Field& field )
{
    if ( FT::empty(field) )
    {
	return Teuchos::ArrayRCP<value_type>(0,0);
    }
    else
    {
	return Teuchos::ArrayRCP<const value_type>(
	    &*FT::begin( field ), 0, FT::size( field ), false );
    }
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
    if ( FT::empty(field) )
    {
	return Teuchos::ArrayRCP<value_type>(0,0);
    }
    else
    {
	return Teuchos::ArrayRCP<value_type>(
	    (value_type*) &*FT::begin( field ), 0, FT::size( field ), false );
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Fill a field with a scalar.
 */
template<class Field>
void FieldTools<Field>::putScalar( Field& field, const value_type& scalar )
{
    std::fill( FT::begin( field ), FT::end( field ), scalar );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Fill a field with a different scalar in each dimension.
 */
template<class Field>
void FieldTools<Field>::putScalar( 
    Field& field, const Teuchos::ArrayView<value_type>& scalars )
{
    int field_dim = FT::dim( field );
    testInvariant( field_dim == (int) scalars.size() );
    for ( int d = 0; d < field_dim; ++d )
    {
	std::fill( dimBegin( field, d ), dimEnd( field, d ), scalars[d] );
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Scale a field by a single value.
 */
template<class Field>
void FieldTools<Field>::scale( Field& field, const value_type& scalar )
{
    iterator field_iterator;
    for ( field_iterator = FT::begin( field ); 
	  field_iterator != FT::end( field );
	  ++field_iterator )
    {
	*field_iterator *= scalar;
    }
}

//---------------------------------------------------------------------------//
/*
 * \brief Scale a field by different value for each dimension.
 */
template<class Field>
void FieldTools<Field>::scale( Field& field, 
			       const Teuchos::ArrayView<value_type>& scalars )
{
    int field_dim = FT::dim( field );
    testPrecondition( field_dim == (int) scalars.size() );
    iterator dim_iterator;
    for ( int d = 0; d < field_dim; ++d )
    {
	for ( dim_iterator = dimBegin( field, d );
	      dim_iterator != dimEnd( field, d );
	      ++dim_iterator )
	{
	    *dim_iterator *= scalars[d];
	}
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compute the global infinity norm for each field dimension.
 */
template<class Field>
void FieldTools<Field>::normInf( const Field& field, const RCP_Comm& comm,
				 Teuchos::Array<value_type>& norms )
{
    int field_dim = FT::dim( field );
    norms.resize( field_dim );
    value_type local_max, local_min, local_norm;
    for ( int d = 0; d < field_dim; ++d )
    {
	local_max = *std::max_element( dimBegin( field, d ),
				       dimEnd( field, d ) );

	local_min = *std::min_element( dimBegin( field, d ),
				       dimEnd( field, d ) );

	local_norm = std::max( std::abs(local_max), std::abs(local_min) );

	Teuchos::reduceAll<int,value_type>( 
	    *comm,
	    Teuchos::REDUCE_MAX,
	    local_norm,
	    Teuchos::Ptr<value_type>( &norms[d] ) );
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compute the global L1 norm for each field dimension.
 */
template<class Field>
void FieldTools<Field>::norm1( const Field& field, const RCP_Comm& comm,
			       Teuchos::Array<value_type>& norms )
{
    int field_dim = FT::dim( field );
    norms.resize( field_dim );
    const_iterator dim_iterator;
    value_type local_norm;
    for ( int d = 0; d < field_dim; ++d )
    {
	local_norm = 0.0;

	for ( dim_iterator = dimBegin( field, d );
	      dim_iterator != dimEnd( field, d );
	      ++dim_iterator )
	{
	    local_norm += std::abs( *dim_iterator );
	}

	Teuchos::reduceAll<int,value_type>( 
	    *comm,
	    Teuchos::REDUCE_SUM,
	    local_norm,
	    Teuchos::Ptr<value_type>( &norms[d] ) );
    } 
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compute the global L2 norm for each field dimension.
 */
template<class Field>
void FieldTools<Field>::norm2( const Field& field, const RCP_Comm& comm,
			       Teuchos::Array<value_type>& norms )
{
    int field_dim = FT::dim( field );
    norms.resize( field_dim );
    const_iterator dim_iterator;
    value_type local_norm;
    for ( int d = 0; d < field_dim; ++d )
    {
	local_norm = 0.0;

	for ( dim_iterator = dimBegin( field, d );
	      dim_iterator != dimEnd( field, d );
	      ++dim_iterator )
	{
	    local_norm += std::abs(*dim_iterator)*std::abs(*dim_iterator);
	}

	Teuchos::reduceAll<int,value_type>( 
	    *comm,
	    Teuchos::REDUCE_SUM,
	    local_norm,
	    Teuchos::Ptr<value_type>( &norms[d] ) );

	norms[d] = std::pow( norms[d], 1.0/2.0 );
    } 
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compute the global q-norm for each field dimension.
 */
template<class Field>
void FieldTools<Field>::normQ( const Field& field, const RCP_Comm& comm, 
			       const int& q,
			       Teuchos::Array<value_type>& norms )
{
    testPrecondition( q > 0 );
    int field_dim = FT::dim( field );
    norms.resize( field_dim );
    const_iterator dim_iterator;
    value_type local_norm, element_product;
    for ( int d = 0; d < field_dim; ++d )
    {
	local_norm = 0.0;

	for ( dim_iterator = dimBegin( field, d );
	      dim_iterator != dimEnd( field, d );
	      ++dim_iterator )
	{
	    element_product = std::abs( *dim_iterator );
	    for ( int n = 0; n < (q-1); ++n )
	    {
		element_product *= std::abs(*dim_iterator);
	    }

	    local_norm += element_product;
	}

	Teuchos::reduceAll<int,value_type>( 
	    *comm,
	    Teuchos::REDUCE_SUM,
	    local_norm,
	    Teuchos::Ptr<value_type>( &norms[d] ) );

	norms[d] = std::pow( norms[d], 1.0/q );
    } 
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compute the global average value for each field dimension.
 */
template<class Field>
void FieldTools<Field>::average( const Field& field, const RCP_Comm& comm,
				 Teuchos::Array<value_type>& averages )
{
    size_type global_length = 
	FieldTools<Field>::globalSize( field, comm );
    testPrecondition( global_length > 0 );

    int field_dim = FT::dim( field );
    averages.resize( field_dim );
    size_type dim_length = global_length / field_dim;
    const_iterator dim_iterator;
    value_type local_sum;
    for ( int d = 0; d < field_dim; ++d )
    {
	local_sum = 0.0;

	for ( dim_iterator = dimBegin( field, d );
	      dim_iterator != dimEnd( field, d );
	      ++dim_iterator )
	{
	    local_sum += *dim_iterator;
	}

	Teuchos::reduceAll<int,value_type>( 
	    *comm,
	    Teuchos::REDUCE_SUM,
	    local_sum,
	    Teuchos::Ptr<value_type>( &averages[d] ) );

	averages[d] /= dim_length;
    } 
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the global size of the field (all dimensions, all instances).
 */
template<class Field>
typename FieldTools<Field>::size_type 
FieldTools<Field>::globalSize( const Field& field, 
			       const RCP_Comm& comm )
{
    size_type local_size = FT::size( field );
    size_type global_size = 0;

    Teuchos::reduceAll<int,size_type>( 
	*comm,
	Teuchos::REDUCE_SUM,
	local_size,
	Teuchos::Ptr<size_type>( &global_size ) );

    return global_size;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the local bounding box for a field of coordinates.
 */
template<class Field>
BoundingBox FieldTools<Field>::coordLocalBoundingBox( const Field& field )
{
    int dim = FT::dim( field );
    testPrecondition( 0 <= dim && dim <= 3 );

    double huge_val = Teuchos::ScalarTraits<double>::rmax();
    double x_min = -huge_val;
    double y_min = -huge_val;
    double z_min = -huge_val;

    double x_max = huge_val;
    double y_max = huge_val;
    double z_max = huge_val;

    if ( dim > 0 )
    {
	x_min = *std::min_element( dimBegin( field, 0 ), dimEnd( field, 0 ) );
	x_max = *std::max_element( dimBegin( field, 0 ), dimEnd( field, 0 ) );
    }
    if ( dim > 1 )
    {
	y_min = *std::min_element( dimBegin( field, 1 ), dimEnd( field, 1 ) );
	y_max = *std::max_element( dimBegin( field, 1 ), dimEnd( field, 1 ) );
    }
    if ( dim > 2 )
    {
	z_min = *std::min_element( dimBegin( field, 2 ), dimEnd( field, 2 ) );
	z_max = *std::max_element( dimBegin( field, 2 ), dimEnd( field, 2 ) );
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
				    local_bounds[0],
				    Teuchos::Ptr<double>( &global_x_min ) );

    Teuchos::reduceAll<int,double>( *comm, 
				    Teuchos::REDUCE_MIN,
				    local_bounds[1],
				    Teuchos::Ptr<double>( &global_y_min ) );

    Teuchos::reduceAll<int,double>( *comm, 
				    Teuchos::REDUCE_MIN,
				    local_bounds[2],
				    Teuchos::Ptr<double>( &global_z_min ) );

    Teuchos::reduceAll<int,double>( *comm, 
				    Teuchos::REDUCE_MAX,
				    local_bounds[3],
				    Teuchos::Ptr<double>( &global_x_max ) );

    Teuchos::reduceAll<int,double>( *comm, 
				    Teuchos::REDUCE_MAX,
				    local_bounds[4],
				    Teuchos::Ptr<double>( &global_y_max ) );

    Teuchos::reduceAll<int,double>( *comm, 
				    Teuchos::REDUCE_MAX,
				    local_bounds[5],
				    Teuchos::Ptr<double>( &global_z_max ) );

    return BoundingBox( global_x_min, global_y_min, global_z_min,
			global_x_max, global_y_max, global_z_max );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_FIELDTOOLS_DEF_HPP

//---------------------------------------------------------------------------//
// end DTK_FieldTools_def.hpp
//---------------------------------------------------------------------------//

