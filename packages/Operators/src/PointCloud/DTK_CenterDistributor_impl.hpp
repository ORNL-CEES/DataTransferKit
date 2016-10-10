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
 * \file   DTK_CenterDistributor_impl.hpp
 * \author Stuart R. Slattery
 * \brief  Global acenter distributor.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_CENTERDISTRIBUTOR_IMPL_HPP
#define DTK_CENTERDISTRIBUTOR_IMPL_HPP

#include <algorithm>
#include <limits>

#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_Array.hpp>

#include "DTK_DBC.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<int DIM>
CenterDistributor<DIM>::CenterDistributor(
        const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
        const Teuchos::ArrayView<const double>& source_centers,
        const Teuchos::ArrayView<const double>& target_centers,
        const double radius,
        Teuchos::Array<double>& target_decomp_source_centers )
    : d_distributor( new Tpetra::Distributor(comm) )
{
    DTK_REQUIRE( 0 == source_centers.size() % DIM );
    DTK_REQUIRE( 0 == target_centers.size() % DIM );

    // Build the import/export data.
    Teuchos::Array<int> export_procs;
    {
        // Compute the radius to expand the local domain with.
        double radius_tol = 1.0e-2;
        double radius_expand = radius * ( 1.0 + radius_tol );

        // Gather the bounding domains for each target proc.
        CloudDomain<DIM> local_target_domain =
            localCloudDomain( target_centers );
        local_target_domain.expand( radius_expand );
        Teuchos::Array<CloudDomain<DIM> > global_target_domains(
            comm->getSize() );
        Teuchos::gatherAll<int,CloudDomain<DIM> >(
            *comm,
            1,
            &local_target_domain,
            global_target_domains.size(),
            global_target_domains.getRawPtr() );

        // Get those that are neighbors to this source proc.
        CloudDomain<DIM> local_source_domain =
            localCloudDomain( source_centers );
        Teuchos::Array<CloudDomain<DIM> > neighbor_target_domains;
        Teuchos::Array<int> neighbor_ranks;
        for ( unsigned i = 0; i < global_target_domains.size(); ++i )
        {
            if ( local_source_domain.checkForIntersection(
                     global_target_domains[i]) )
            {
                neighbor_target_domains.push_back(global_target_domains[i]);
                neighbor_ranks.push_back(i);
            }
        }
        global_target_domains.clear();

        // Find the procs to which the sources will be sent.
        Teuchos::ArrayView<const double> source_point;
        for ( unsigned source_id = 0;
              source_id < source_centers.size() / DIM;
              ++source_id )
        {
            source_point = source_centers.view( DIM*source_id, DIM );
            for ( unsigned b = 0; b < neighbor_target_domains.size(); ++b )
            {
                if ( neighbor_target_domains[b].pointInDomain(source_point) )
                {
                    export_procs.push_back( neighbor_ranks[b] );
                    d_export_ids.push_back( source_id );
                }
            }
        }
    }
    DTK_CHECK( d_export_ids.size() == export_procs.size() );
    d_num_exports = d_export_ids.size();

    // Create the communication plan.
    Teuchos::ArrayView<int> export_procs_view = export_procs();
    d_num_imports = d_distributor->createFromSends( export_procs_view );
    export_procs.clear();


    // Unroll the coordinates to handle cases where single source centers
    // may have multiple destinations.
    Teuchos::Array<unsigned>::const_iterator export_id_it;
    Teuchos::Array<double> src_coords( d_num_exports * DIM );
    for ( int n = 0; n < d_num_exports; ++n )
    {
        src_coords( DIM*n, DIM ).assign(
            source_centers(DIM*d_export_ids[n],DIM) );
    }

    // Move the source center coordinates to the target decomposition.
    Teuchos::ArrayView<const double> src_coords_view = src_coords();
    target_decomp_source_centers.resize( d_num_imports * DIM );
    d_distributor->doPostsAndWaits(
        src_coords_view, DIM, target_decomp_source_centers() );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Given a set of scalar values at the given source centers in the
 * source decomposition, distribute them to the target decomposition.
 */
template<int DIM>
template<class T>
void CenterDistributor<DIM>::distribute(
    const Teuchos::ArrayView<const T>& source_decomp_data,
    const Teuchos::ArrayView<T>& target_decomp_data ) const
{
    DTK_REQUIRE( d_num_imports == target_decomp_data.size() );

    // Unroll the source data to handle cases where single data points may
    // have multiple destinations.
    Teuchos::Array<unsigned>::const_iterator export_id_it;
    Teuchos::Array<T> src_data( d_num_exports );
    typename Teuchos::Array<T>::iterator src_it;
    for ( export_id_it = d_export_ids.begin(),
                src_it = src_data.begin();
          export_id_it != d_export_ids.end();
          ++export_id_it, ++src_it )
    {
        DTK_CHECK( *export_id_it < source_decomp_data.size() );
        *src_it = source_decomp_data[ *export_id_it ];
    }

    // Distribute.
    Teuchos::ArrayView<const T> src_view = src_data();
    d_distributor->doPostsAndWaits( src_view, 1, target_decomp_data );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compute the bounding domain of the local set of centers.
 */
template<int DIM>
CloudDomain<DIM> CenterDistributor<DIM>::localCloudDomain(
    const Teuchos::ArrayView<const double>& centers ) const
{
    Teuchos::Array<double> bounds( 2*DIM, 0.0 );

    if ( centers.size() > 0 )
    {
        for ( int d = 0; d < DIM; ++d )
        {
            bounds[2*d] = std::numeric_limits<double>::max();
            bounds[2*d+1] = std::numeric_limits<double>::min();
        }
    }

    Teuchos::ArrayView<const double>::const_iterator center_it;
    for ( center_it = centers.begin(); center_it != centers.end(); )
    {
        for ( int d = 0; d < DIM; ++d )
        {
            bounds[2*d] = std::min( bounds[2*d], *(center_it+d) );
            bounds[2*d+1] = std::max( bounds[2*d+1], *(center_it+d) );
        }
        center_it += DIM;
    }

    return CloudDomain<DIM>( bounds.getRawPtr() );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_CENTERDISTRIBUTOR_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_CenterDistributor_impl.hpp
//---------------------------------------------------------------------------//

