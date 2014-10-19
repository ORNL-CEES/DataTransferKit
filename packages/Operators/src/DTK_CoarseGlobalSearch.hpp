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
 * \file DTK_CoarseGlobalSearch.hpp
 * \author Stuart R. Slattery
 * \brief CoarseGlobalSearch declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_COARSEGLOBALSEARCH_HPP
#define DTK_COARSEGLOBALSEARCH_HPP

#include <unordered_map>

#include "DTK_Types.hpp"
#include "DTK_EntityIterator.hpp"
#include "DTK_EntityGlobalMap.hpp"
#include "DTK_StaticSearchTree.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_Tuple.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \class CoarseGlobalSearch 
 * \brief A CoarseGlobalSearch data structure for global entity coarse search.
 */
//---------------------------------------------------------------------------//
class CoarseGlobalSearch
{
  public:

    // Constructor.
    CoarseGlobalSearch( const EntityIterator& domain_iterator,
			const Teuchos::ParameterList& parameters ); 

    // Destructor.
    ~CoarseGlobalSearch();

    // Redistribute a set of range entity centroid coordinates to the correct
    // process.
    void search( const EntityIterator& range_iterator,
		 const Teuchos::RCP<EntityLocalMap>& range_local_map,
		 const Teuchos::ParameterList& parameters,
		 Teuchos::Array<double>& range_centroids ) const;

  private:

    // Assemble the local bounding box around an iterator.
    void assembleBoundingBox( const EntityIterator& entity_iterator,
			      Teuchos::Tuple<double,6>& bounding_box ) const;
    
  private:

    // Domain bounding boxes.
    Teuchos::Array<Teuchos::Tuple<double,6> > d_domain_boxes;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // DTK_COARSEGLOBALSEARCH_HPP

//---------------------------------------------------------------------------//
// end CoarseGlobalSearch.hpp
//---------------------------------------------------------------------------//

