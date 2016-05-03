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
 * \file DTK_POD_PointCloudEntity.hpp
 * \author Stuart R. Slattery
 * \brief POD_PointCloudEntity declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_POD_POINTCLOUDENTITY_HPP
#define DTK_POD_POINTCLOUDENTITY_HPP

#include <iostream>

#include "DTK_Entity.hpp"
#include "DTK_POD_PointCloudEntityImpl.hpp"
#include "DTK_POD_Types.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class POD_PointCloudEntity
  \brief POD_PointCloudEntity interface.
  
  POD_PointCloudEntity gives an interface for entities in POD point clouds.
*/
//---------------------------------------------------------------------------//
class POD_PointCloudEntity : public Entity
{
  public:

    // Default constructor.
    POD_PointCloudEntity( const double* cloud_coords,
                          const unsigned num_points,
                          const int space_dim,
                          const DataLayout layout,
                          const EntityId global_id,
                          const int local_id,
                          const int owner_rank );
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_POD_POINTCLOUDENTITY_HPP

//---------------------------------------------------------------------------//
// end DTK_POD_PointCloudEntity.hpp
//---------------------------------------------------------------------------//

