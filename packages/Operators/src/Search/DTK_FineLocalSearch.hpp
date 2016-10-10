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
 * \file DTK_FineLocalSearch.hpp
 * \author Stuart R. Slattery
 * \brief FineLocalSearch declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_FINELOCALSEARCH_HPP
#define DTK_FINELOCALSEARCH_HPP

#include <unordered_map>

#include "DTK_Types.hpp"
#include "DTK_Entity.hpp"
#include "DTK_EntityLocalMap.hpp"

#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_ParameterList.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class FineLocalSearch
  \brief A FineLocalSearch data structure for local entity fine search.

  Find the entites in a subset into which a point parametrically maps.
 */
//---------------------------------------------------------------------------//
class FineLocalSearch
{
  public:

    // Constructor.
    FineLocalSearch( const Teuchos::RCP<EntityLocalMap>& local_map );

    // Find the set of entities to which a point maps.
    void search( const Teuchos::ArrayView<const Entity>& neighbors,
                 const Teuchos::ArrayView<const double>& point,
                 const Teuchos::ParameterList& parameters,
                 Teuchos::Array<Entity>& parents,
                 Teuchos::Array<double>& reference_coordinates ) const;

  private:

    // Local map for the fine search.
    Teuchos::RCP<EntityLocalMap> d_local_map;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // DTK_FINELOCALSEARCH_HPP

//---------------------------------------------------------------------------//
// end FineLocalSearch.hpp
//---------------------------------------------------------------------------//

