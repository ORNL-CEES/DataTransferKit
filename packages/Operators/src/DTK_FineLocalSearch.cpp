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
 * \file DTK_FineLocalSearch.cpp
 * \author Stuart R. Slattery
 * \brief FineLocalSearch definition.
 */
//---------------------------------------------------------------------------//

#include <cmath>
#include <algorithm>

#include "DTK_FineLocalSearch.hpp"
#include "DTK_DBC.hpp"
#include "DTK_SearchTreeFactory.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
FineLocalSearch::FineLocalSearch( 
    const Teuchos::RCP<EntityLocalMap>& local_map )
    : d_local_map( local_map )
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
FineLocalSearch::~FineLocalSearch()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Find the set of entities to which a point maps.
 */
void FineLocalSearch::search( 
    const Teuchos::ArrayView<const Entity>& neighbors,
    const Teuchos::ArrayView<const double>& point,
    const Teuchos::ParameterList& parameters,
    Teuchos::Array<Entity>& parents,
    Teuchos::Array<double>& reference_coordinates ) const
{
    bool return_all = false;
    if ( parameters.isParameter("Fine Local Search Return All") )
    {
	return_all = parameters.get<bool>("Fine Local Search Return All");
    }
    
    parents.clear();
    reference_coordinates.clear();
    int physical_dim = point.size();
    Teuchos::Array<double> ref_point( physical_dim );
    Teuchos::ArrayView<const Entity>::const_iterator neighbor_it;
    for ( neighbor_it = neighbors.begin();
	  neighbor_it != neighbors.end();
	  ++neighbor_it )
    {
	DTK_ENSURE( neighbor_it->physicalDimension() == point.size() );

	if ( d_local_map->isSafeToMapToReferenceFrame(*neighbor_it,point) )
	{
	    if ( d_local_map->mapToReferenceFrame(
		     *neighbor_it,point,ref_point()) )
	    {
		if ( d_local_map->checkPointInclusion(
			 *neighbor_it,ref_point()) )
		{
		    parents.push_back( *neighbor_it );
		    for ( int d = 0; d < physical_dim; ++d )
		    {
			reference_coordinates.push_back( ref_point[d] );
		    }
		    if ( !return_all )
		    {
			return;
		    }
		}
	    }
	}
    }
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_FineLocalSearch.cpp
//---------------------------------------------------------------------------//

