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

#include "DTK_DBC.hpp"

#include <Teuchos_Tuple.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_Ptr.hpp>
#include <Teuchos_as.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Get the local size of the dimensions.
 *
 * \param field The field to get the dimension size of.
 *
 * \return The size of the dimensions in the field.
 */
template<class Field>
typename FieldTools<Field>::size_type
FieldTools<Field>::dimSize( const Field& field )
{
    DTK_REQUIRE( FT::dim( field ) > 0 );
    size_type field_size = FT::size( field );
    return field_size / FT::dim( field );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get an iterator to the beginning of a dimension.
 *
 * \param field The field to get an iterator for.
 *
 * \param dim The dimension to get an iterator for.
 *
 * \return The field iterator the beginning given of the dimension
 */
template<class Field>
typename FieldTools<Field>::iterator 
FieldTools<Field>::dimBegin( Field& field, const int dim )
{
    DTK_REQUIRE( FT::dim( field ) > 0 );
    DTK_REQUIRE( dim >= 0 && dim < FT::dim( field ) );
    return FT::begin(field) + dim*dimSize(field);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get a const iterator to the beginning of a dimension.
 *
 * \param field The field to get an iterator for.
 *
 * \param dim The dimension to get an iterator for.
 *
 * \return The const field iterator the beginning of the given dimension
 */
template<class Field>
typename FieldTools<Field>::const_iterator 
FieldTools<Field>::dimBegin( const Field& field, const int dim )
{
    DTK_REQUIRE( FT::dim( field ) > 0 );
    DTK_REQUIRE( dim >= 0 && dim < FT::dim( field ) );
    return FT::begin(field) + dim*dimSize(field);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get an iterator to the end of a dimension.
 *
 * \param field The field to get an iterator for.
 *
 * \param dim The dimension to get an iterator for.
 *
 * \return The field iterator the end of the given dimension
 */
template<class Field>
typename FieldTools<Field>::iterator 
FieldTools<Field>::dimEnd( Field& field, const int dim )
{
    DTK_REQUIRE( FT::dim( field ) > 0 );
    DTK_REQUIRE( dim >= 0 && dim < FT::dim( field ) );
    return FT::begin(field) + (dim+1)*dimSize(field);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get a const iterator to the end of a dimension.
 *
 * \param field The field to get an iterator for.
 *
 * \param dim The dimension to get an iterator for.
 *
 * \return The const field iterator the end of the given dimension
 */
template<class Field>
typename FieldTools<Field>::const_iterator 
FieldTools<Field>::dimEnd( const Field& field, const int dim )
{
    DTK_REQUIRE( FT::dim( field ) > 0 );
    DTK_REQUIRE( dim >= 0 && dim < FT::dim( field ) );
    return FT::begin(field) + (dim+1)*dimSize(field);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get a const view of the field. The ArrayRCP object will not manage
 * the memory. 
 *
 * \param field The field to get a view of.
 *
 * \return A const view of the field.
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
 *
 * \param field The field to get a view of.
 *
 * \return A non-const view of the field.
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
	    const_cast<value_type*>(&*FT::begin( field )), 0, 
	    FT::size( field ), false );
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get a deep-copy of the field. The arrayRCP object will manage the
 * memory.
 *
 * \param field The field to copy.
 *
 * \return A copy of the field.
 */
template<class Field>
Teuchos::ArrayRCP<typename FieldTools<Field>::value_type>
FieldTools<Field>::copy( const Field& field )
{
    if ( FT::empty(field) )
    {
	return Teuchos::ArrayRCP<value_type>(0,0);
    }
    else
    {
	Teuchos::ArrayRCP<value_type> field_copy( FT::size(field) );
	std::copy( FT::begin(field), FT::end(field), field_copy.begin() );
	return field_copy;
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get a const view of a dimension of the field. The ArrayRCP object
 * will not manage the memory.
 *
 * \param field The field to get a view of.
 *
 * \param dim The dimension to get a view of.
 *
 * \return A const view of the field.
 */
template<class Field>
Teuchos::ArrayRCP<const typename FieldTools<Field>::value_type>
FieldTools<Field>::dimView( const Field& field, const int dim )
{
    DTK_REQUIRE( FT::dim( field ) > 0 );
    DTK_REQUIRE( dim >= 0 && dim < FT::dim( field ) );

    if ( FT::empty(field) )
    {
	return Teuchos::ArrayRCP<value_type>(0,0);
    }
    else
    {
	size_type dim_size = dimSize( field );
	Teuchos::ArrayRCP<const value_type> field_view = view( field );
	return field_view.persistingView( dim*dim_size, dim_size );
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get a non-const view of a dimension of the field. The ArrayRCP
 * object will not manage the memory.
 *
 * \param field The field to get a view of.
 *
 * \param dim The dimension to get a view of.
 *
 * \return A non-const view of the field.
 */
template<class Field>
Teuchos::ArrayRCP<typename FieldTools<Field>::value_type>
FieldTools<Field>::dimNonConstView( const Field& field, const int dim )
{
    DTK_REQUIRE( FT::dim( field ) > 0 );
    DTK_REQUIRE( dim >= 0 && dim < FT::dim( field ) );

    if ( FT::empty(field) )
    {
	return Teuchos::ArrayRCP<value_type>(0,0);
    }
    else
    {
	size_type dim_size = dimSize( field );
	Teuchos::ArrayRCP<const value_type> field_view = view( field );
	return field_view.persistingView( dim*dim_size, dim_size );
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Fill a field with a scalar.
 * 
 * \param field The field to fill.
 *
 * \param scalar The scalar to fill the field with.
 */
template<class Field>
void FieldTools<Field>::putScalar( Field& field, const value_type& scalar )
{
    std::fill( FT::begin( field ), FT::end( field ), scalar );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Fill a field with a different scalar in each dimension.
 * 
 * \param field The field to fill.
 *
 * \param scalars The array of scalars to fill the field with. This array must
 * be of the same length as the field dimension. Each entry in the array
 * correlates to the scalar for a particular field dimension.
 */
template<class Field>
void FieldTools<Field>::putScalar( 
    Field& field, const Teuchos::ArrayView<value_type>& scalars )
{
    DTK_CHECK( FT::dim( field ) == Teuchos::as<int>(scalars.size()) );
    for ( int d = 0; d < FT::dim( field ); ++d )
    {
	std::fill( dimBegin( field, d ), dimEnd( field, d ), scalars[d] );
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Scale a field by a single value.
 * 
 * \param field The field to scale.
 *
 * \param The scalar to scale the field with.
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
 * 
 * \param field The field to scale.
 *
 * \param scalars The array of scalars to scale the field with. This array must
 * be of the same length as the field dimension. Each entry in the array
 * correlates to the scalar for a particular field dimension.
 */
template<class Field>
void FieldTools<Field>::scale( Field& field, 
			       const Teuchos::ArrayView<value_type>& scalars )
{
    DTK_REQUIRE( FT::dim( field ) == Teuchos::as<int>(scalars.size()) );
    iterator dim_iterator;
    for ( int d = 0; d < FT::dim( field ); ++d )
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
 *
 * \param field The field to compute the norm of.
 *
 * \param comm The communicator over which the field is defined.
 *
 * \param norms The norms for each dimension in the field. This array will be
 * of the same length as the field dimension.
 */
template<class Field>
void FieldTools<Field>::normInf( const Field& field, const RCP_Comm& comm,
				 Teuchos::Array<value_type>& norms )
{
    norms.resize( FT::dim( field ) );
    value_type local_max, local_min, local_norm;
    for ( int d = 0; d < FT::dim( field ); ++d )
    {
	norms[d] = 0.0;

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
 *
 * \param field The field to compute the norm of.
 *
 * \param comm The communicator over which the field is defined.
 *
 * \param norms The norms for each dimension in the field. This array will be
 * of the same length as the field dimension.
 */
template<class Field>
void FieldTools<Field>::norm1( const Field& field, const RCP_Comm& comm,
			       Teuchos::Array<value_type>& norms )
{
    norms.resize( FT::dim( field ) );
    const_iterator dim_iterator;
    value_type local_norm;
    for ( int d = 0; d < FT::dim( field ); ++d )
    {
	norms[d] = 0.0;
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
 *
 * \param field The field to compute the norm of.
 *
 * \param comm The communicator over which the field is defined.
 *
 * \param norms The norms for each dimension in the field. This array will be
 * of the same length as the field dimension.
 */
template<class Field>
void FieldTools<Field>::norm2( const Field& field, const RCP_Comm& comm,
			       Teuchos::Array<value_type>& norms )
{
    norms.resize( FT::dim( field ) );
    const_iterator dim_iterator;
    value_type local_norm;
    for ( int d = 0; d < FT::dim( field ); ++d )
    {
	norms[d] = 0.0;
	local_norm = 0.0;

	for ( dim_iterator = dimBegin( field, d );
	      dim_iterator != dimEnd( field, d );
	      ++dim_iterator )
	{
	    local_norm += std::abs( (*dim_iterator)*(*dim_iterator) );
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
 *
 * \param field The field to compute the norm of.
 *
 * \param comm The communicator over which the field is defined.
 * 
 * \param q The integer norm to compute.
 *
 * \param norms The norms for each dimension in the field. This array will be
 * of the same length as the field dimension.
 */
template<class Field>
void FieldTools<Field>::normQ( const Field& field, const RCP_Comm& comm, 
			       const int q,
			       Teuchos::Array<value_type>& norms )
{
    DTK_REQUIRE( q > 0 );
    norms.resize( FT::dim( field ) );
    const_iterator dim_iterator;
    value_type local_norm, element_product;
    for ( int d = 0; d < FT::dim( field ); ++d )
    {
	norms[d] = 0.0;
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
 *
 * \param field The field to compute the average of.
 *
 * \param comm The communicator over which the field is defined.
 *
 * \param averages The averages for each dimension in the field. This array
 * will bep of the same length as the field dimension.
 */
template<class Field>
void FieldTools<Field>::average( const Field& field, const RCP_Comm& comm,
				 Teuchos::Array<value_type>& averages )
{
    size_type global_length = 
	FieldTools<Field>::globalSize( field, comm );
    DTK_REQUIRE( global_length > 0 );

    averages.resize( FT::dim( field ) );
    size_type dim_length = global_length / FT::dim( field );
    const_iterator dim_iterator;
    value_type local_sum;
    for ( int d = 0; d < FT::dim( field ); ++d )
    {
	averages[d] = 0.0;
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
 * 
 * \param field The field to get the global size of.
 *
 * \param comm The communicator over which the field is defined.
 *
 * \return The global size of the field. This is equivalent to all global
 * elements in the field in all dimensions.
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
 * \brief Get the local bounding box for a field of coordinates. Here the
 * field is directly interpreted as a coordinate field.
 *
 * \param field The coordinate field to compute the local box for.
 *
 * \return The local bounding box of the coordinate field.
 */
template<class Field>
BoundingBox FieldTools<Field>::coordLocalBoundingBox( const Field& field )
{
    DTK_REQUIRE( !FT::empty(field) );
    int dim = FT::dim( field );
    DTK_REQUIRE( 0 <= dim && dim <= 3 );

    double huge_val = Teuchos::ScalarTraits<double>::rmax();
    double x_min = -huge_val;
    double y_min = -huge_val;
    double z_min = -huge_val;

    double x_max = huge_val;
    double y_max = huge_val;
    double z_max = huge_val;

    if ( dim > 0 )
    {
	x_min = *std::min_element( dimBegin(field, 0), dimEnd(field, 0) );
	x_max = *std::max_element( dimBegin(field, 0), dimEnd(field, 0) );
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
 * \brief Get the global bounding box for a field of coordinates. Here the
 * field is directly interpreted as a coordinate field.
 *
 * \param field The coordinate field to compute the global box for.
 *
 * \param comm The communicator the field is defined over.
 *
 * \return The global bounding box of the coordinate field.
 */
template<class Field>
BoundingBox FieldTools<Field>::coordGlobalBoundingBox( const Field& field,
						       const RCP_Comm& comm )
{
    Teuchos::Tuple<double,6> local_bounds =
	Teuchos::tuple( Teuchos::ScalarTraits<double>::rmax(),
			Teuchos::ScalarTraits<double>::rmax(),
			Teuchos::ScalarTraits<double>::rmax(),
			-Teuchos::ScalarTraits<double>::rmax(),
			-Teuchos::ScalarTraits<double>::rmax(),
			-Teuchos::ScalarTraits<double>::rmax() );
    if ( !FT::empty(field) )
    {
	BoundingBox local_box = coordLocalBoundingBox( field );
	local_bounds = local_box.getBounds();
    }

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

