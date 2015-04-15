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
 * \brief DTK_IntegrationPointSet.hpp
 * \author Stuart R. Slattery
 * \brief Integration point set.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_INTEGRATIONPOINTSET_HPP
#define DTK_INTEGRATIONPOINTSET_HPP

#include "DTK_Types.hpp"
#include "DTK_EntitySet.hpp"
#include "DTK_IntegrationPoint.hpp"

#include <Teuchos_ArrayView.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class IntegrationPointSet
  \brief EntitySet of integration points.

  IntegrationPointSet is a special entity set of integration points.
*/
//---------------------------------------------------------------------------//
class IntegrationPointSet : public EntitySet
{
  public:

    /*!
     * \brief Constructor.
     */
    IntegrationPointSet( const Teuchos::RCP<const Teuchos::Comm<int> >& comm );

    // Add an integration point to the set.
    void addPoint( const IntegrationPoint& ip );

    // Finalize the point set to construct global ids.
    void finalize();
    
  private:

    // Communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > d_comm;

    // Local integration points.
    Teuchos::Array<IntegrationPoint> d_points;
    
    // Global ids of the local integration points.
    Teuchos::Array<EntityId> d_ip_global_ids;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_INTEGRATIONPOINTSET_HPP

//---------------------------------------------------------------------------//
// end DTK_IntegrationPointSet.hpp
//---------------------------------------------------------------------------//
