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
 * \brief DTK_IntegrationPointSet.cpp
 * \author Stuart R. Slattery
 * \brief Integration point set.
 */
//---------------------------------------------------------------------------//

#include "DTK_IntegrationPointSet.hpp"
#include "DTK_DBC.hpp"

#include <Teuchos_CommHelpers.hpp>

#include <Tpetra_Map.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// IntegrationPointSetIterator implementation.
//---------------------------------------------------------------------------//
// Default constructor.
IntegrationPointSetIterator::IntegrationPointSetIterator() { /* ... */}

//---------------------------------------------------------------------------//
// Constructor.
IntegrationPointSetIterator::IntegrationPointSetIterator(
    Teuchos::RCP<Teuchos::Array<IntegrationPoint>> points )
    : d_points( points )
    , d_points_it( d_points->begin() )
{ /* ... */
}

//---------------------------------------------------------------------------//
// Copy constructor.
IntegrationPointSetIterator::IntegrationPointSetIterator(
    const IntegrationPointSetIterator &rhs )
    : d_points( rhs.d_points )
    , d_points_it( rhs.d_points_it )
{ /* ... */
}

//---------------------------------------------------------------------------//
// Assignment operator.
IntegrationPointSetIterator &IntegrationPointSetIterator::
operator=( const IntegrationPointSetIterator &rhs )
{
    if ( &rhs == this )
    {
        return *this;
    }
    d_points = rhs.d_points;
    d_points_it = rhs.d_points_it;
    return *this;
}

//---------------------------------------------------------------------------//
// Pre-increment operator.
EntityIterator &IntegrationPointSetIterator::operator++()
{
    ++d_points_it;
    return *this;
}

//---------------------------------------------------------------------------//
// Dereference operator.
Entity &IntegrationPointSetIterator::operator*( void )
{
    this->operator->();
    return d_current_entity;
}

//---------------------------------------------------------------------------//
// Dereference operator.
Entity *IntegrationPointSetIterator::operator->( void )
{
    d_current_entity =
        IntegrationPointEntity( Teuchos::ptrFromRef( *d_points_it ) );
    return &d_current_entity;
}

//---------------------------------------------------------------------------//
// Equal comparison operator.
bool IntegrationPointSetIterator::operator==( const EntityIterator &rhs ) const
{
    const IntegrationPointSetIterator *rhs_vec =
        static_cast<const IntegrationPointSetIterator *>( &rhs );
    const IntegrationPointSetIterator *rhs_vec_impl =
        static_cast<const IntegrationPointSetIterator *>(
            rhs_vec->b_iterator_impl.get() );
    return ( rhs_vec_impl->d_points_it == d_points_it );
}

//---------------------------------------------------------------------------//
// Not equal comparison operator.
bool IntegrationPointSetIterator::operator!=( const EntityIterator &rhs ) const
{
    const IntegrationPointSetIterator *rhs_vec =
        static_cast<const IntegrationPointSetIterator *>( &rhs );
    const IntegrationPointSetIterator *rhs_vec_impl =
        static_cast<const IntegrationPointSetIterator *>(
            rhs_vec->b_iterator_impl.get() );
    return ( rhs_vec_impl->d_points_it != d_points_it );
}

//---------------------------------------------------------------------------//
// An iterator assigned to the beginning.
EntityIterator IntegrationPointSetIterator::begin() const
{
    return IntegrationPointSetIterator( d_points );
}

//---------------------------------------------------------------------------//
// An iterator assigned to the end.
EntityIterator IntegrationPointSetIterator::end() const
{
    IntegrationPointSetIterator end_it( d_points );
    end_it.d_points_it = d_points->end();
    return end_it;
}

//---------------------------------------------------------------------------//
// Create a clone of the iterator. We need this for the copy constructor
// and assignment operator to pass along the underlying implementation.
std::unique_ptr<EntityIterator> IntegrationPointSetIterator::clone() const
{
    return std::unique_ptr<EntityIterator>(
        new IntegrationPointSetIterator( *this ) );
}

//---------------------------------------------------------------------------//
// IntegrationPointSet Implementation
//---------------------------------------------------------------------------//
// Constructor.
IntegrationPointSet::IntegrationPointSet(
    const Teuchos::RCP<const Teuchos::Comm<int>> &comm )
    : d_comm( comm )
{ /* ... */
}

//---------------------------------------------------------------------------//
// Add an integration point to the set.
void IntegrationPointSet::addPoint( const IntegrationPoint &ip )
{
    d_points.push_back( ip );
}

//---------------------------------------------------------------------------//
// Finalize the point set to construct global ids.
void IntegrationPointSet::finalize()
{
    // Build a globally contiguous ordering of point global ids.
    EntityId num_local_ip = d_points.size();

    EntityId invalid = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
    Teuchos::RCP<const Tpetra::Map<int, EntityId>> map =
        Tpetra::createContigMap<int, EntityId>( invalid, num_local_ip, d_comm );

    // Get the starting id for this node.
    DTK_CHECK( map->isContiguous() );
    DTK_CHECK( map->getMaxGlobalIndex() - map->getMinGlobalIndex() + 1 ==
               num_local_ip );
    DTK_CHECK( map->getNodeNumElements() == num_local_ip );
    d_start_gid = map->getMinGlobalIndex();

    // Assign global ids to the points.
    Teuchos::ArrayView<const EntityId> ip_gids = map->getNodeElementList();
    for ( EntityId p = 0; p < num_local_ip; ++p )
    {
        d_points[p].d_gid = ip_gids[p];
    }
}

//---------------------------------------------------------------------------//
// Get an integration point with the given global id.
const IntegrationPoint &
IntegrationPointSet::getPoint( const EntityId ip_id ) const
{
    DTK_REQUIRE( ip_id >= d_start_gid );
    DTK_REQUIRE( ip_id - d_start_gid <
                 Teuchos::as<EntityId>( d_points.size() ) );
    return d_points[ip_id - d_start_gid];
}

//---------------------------------------------------------------------------//
// Get an entity iterator over the integration points.
EntityIterator IntegrationPointSet::entityIterator() const
{
    return IntegrationPointSetIterator( Teuchos::rcpFromRef( d_points ) );
}

//---------------------------------------------------------------------------//
// Get the global maximum support size for all integration points.
int IntegrationPointSet::globalMaxSupportSize() const
{
    int local_max = 0;
    if ( d_points.size() > 0 )
    {
        auto comp = []( const IntegrationPoint &a, const IntegrationPoint &b ) {
            return ( a.d_owner_support_ids.size() <
                     b.d_owner_support_ids.size() );
        };
        local_max = std::max_element( d_points.begin(), d_points.end(), comp )
                        ->d_owner_support_ids.size();
    }
    int global_max = 0;
    Teuchos::reduceAll( *d_comm, Teuchos::REDUCE_MAX, local_max,
                        Teuchos::ptrFromRef( global_max ) );
    return global_max;
}

//---------------------------------------------------------------------------//
// Get the centroid of the integration point.
void IntegrationPointSet::centroid(
    const Entity &entity, const Teuchos::ArrayView<double> &centroid ) const
{
    const IntegrationPoint &point = getPoint( entity.id() );
    centroid.assign( point.d_physical_coordinates );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_IntegrationPointSet.cpp
//---------------------------------------------------------------------------//
