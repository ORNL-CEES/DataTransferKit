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
 * \brief DTK_POD_PointCloudEntityIterator.cpp
 * \author Stuart R. Slattery
 * \brief Entity iterator interface.
 */
//---------------------------------------------------------------------------//

#include "DTK_DBC.hpp"
#include "DTK_POD_PointCloudEntityIterator.hpp"
#include <DTK_POD_PointCloudEntity.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Default constructor.
POD_PointCloudEntityIterator::POD_PointCloudEntityIterator()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Constructor.
POD_PointCloudEntityIterator::POD_PointCloudEntityIterator(
    const double* cloud_coords,
    const EntityId* global_ids,                              
    const unsigned num_points,
    const int space_dim,
    const DataLayout layout,
    const int my_rank,
    const PredicateFunction& predicate )
    : d_cloud_coords( cloud_coords )
    , d_global_ids( global_ids )
    , d_num_points( num_points )
    , d_space_dim( space_dim )
    , d_layout( layout )
    , d_my_rank( my_rank )
    , d_current_lid( 0 )
{
    this->b_predicate = predicate;
}

//---------------------------------------------------------------------------//
// Copy constructor.
POD_PointCloudEntityIterator::POD_PointCloudEntityIterator( 
    const POD_PointCloudEntityIterator& rhs )
    : d_cloud_coords( rhs.d_cloud_coords )
    , d_global_ids( rhs.d_global_ids )
    , d_num_points( rhs.d_num_points )
    , d_space_dim( rhs.d_space_dim )
    , d_layout( rhs.d_layout )
    , d_my_rank( rhs.d_my_rank )
    , d_current_lid( rhs.d_current_lid )
{
    this->b_predicate = rhs.b_predicate;
}

//---------------------------------------------------------------------------//
// Assignment operator.
POD_PointCloudEntityIterator& POD_PointCloudEntityIterator::operator=( 
    const POD_PointCloudEntityIterator& rhs )
{
    this->b_predicate = rhs.b_predicate;
    if ( &rhs == this )
    {
	return *this;
    }
    d_cloud_coords = rhs.d_cloud_coords;
    d_global_ids = rhs.d_global_ids;
    d_num_points = rhs.d_num_points;
    d_space_dim = rhs.d_space_dim;
    d_layout = rhs.d_layout;
    d_my_rank = rhs.d_my_rank;
    d_current_lid = rhs.d_current_lid;
    return *this;
}

//---------------------------------------------------------------------------//
// Pre-increment operator.
EntityIterator& POD_PointCloudEntityIterator::operator++()
{
    ++d_current_lid;
    return *this;
}

//---------------------------------------------------------------------------//
// Dereference operator.
Entity& POD_PointCloudEntityIterator::operator*(void)
{
    this->operator->();
    return d_current_entity;
}

//---------------------------------------------------------------------------//
// Dereference operator.
Entity* POD_PointCloudEntityIterator::operator->(void)
{
    DTK_CHECK( d_current_lid < d_num_points );
    d_current_entity = 
	POD_PointCloudEntity( d_cloud_coords, d_num_points, d_space_dim,
                          d_layout, d_global_ids[d_current_lid],
                          d_current_lid, d_my_rank );
    return &d_current_entity;
}

//---------------------------------------------------------------------------//
// Equal comparison operator.
bool POD_PointCloudEntityIterator::operator==( 
    const EntityIterator& rhs ) const
{ 
    const POD_PointCloudEntityIterator* rhs_it = 
	static_cast<const POD_PointCloudEntityIterator*>(&rhs);
    const POD_PointCloudEntityIterator* rhs_it_impl = 
	static_cast<const POD_PointCloudEntityIterator*>(rhs_it->b_iterator_impl.get());
    return ( rhs_it_impl->d_current_lid == d_current_lid );
}

//---------------------------------------------------------------------------//
// Not equal comparison operator.
bool POD_PointCloudEntityIterator::operator!=( 
    const EntityIterator& rhs ) const
{
    const POD_PointCloudEntityIterator* rhs_it = 
	static_cast<const POD_PointCloudEntityIterator*>(&rhs);
    const POD_PointCloudEntityIterator* rhs_it_impl = 
	static_cast<const POD_PointCloudEntityIterator*>(rhs_it->b_iterator_impl.get());
    return ( rhs_it_impl->d_current_lid != d_current_lid );
}

//---------------------------------------------------------------------------//
// An iterator assigned to the beginning.
EntityIterator POD_PointCloudEntityIterator::begin() const
{ 
    return POD_PointCloudEntityIterator( 
	d_cloud_coords, d_global_ids, d_num_points, d_space_dim,
        d_layout, d_my_rank, this->b_predicate );
}

//---------------------------------------------------------------------------//
// An iterator assigned to the end.
EntityIterator POD_PointCloudEntityIterator::end() const
{
    POD_PointCloudEntityIterator end_it( 
	d_cloud_coords, d_global_ids, d_num_points, d_space_dim,
        d_layout, d_my_rank, this->b_predicate );
    end_it.d_current_lid = end_it.d_num_points;
    return end_it;
}

//---------------------------------------------------------------------------//
// Create a clone of the iterator. We need this for the copy constructor
// and assignment operator to pass along the underlying implementation.
std::unique_ptr<EntityIterator> POD_PointCloudEntityIterator::clone() const
{
    return std::unique_ptr<EntityIterator>( new POD_PointCloudEntityIterator(*this) );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_POD_PointCloudEntityIterator.cpp
//---------------------------------------------------------------------------//
