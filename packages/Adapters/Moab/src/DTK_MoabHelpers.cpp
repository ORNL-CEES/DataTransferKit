//---------------------------------------------------------------------------//
/*
  Copyright (c) 2014, Stuart R. Slattery
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

  *: Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

  *: Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

  *: Neither the name of the Oak Ridge National Laboratory nor the
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
 * \brief DTK_MoabHelpers.cpp
 * \author Stuart R. Slattery
 * \brief Moab helper functions.
 */
//---------------------------------------------------------------------------//

#include "DTK_MoabHelpers.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Given a Moab EntityType, get the DTK EntityType.
EntityType MoabHelpers::getEntityTypeFromMoabType( 
    const moab::EntityType moab_type )
{
    EntityType dtk_type = ENTITY_TYPE_INVALID;
    switch( moab_type )
    {
	case moab::MBVERTEX:
	    dtk_type = ENTITY_TYPE_NODE;
	    break;
	case moab::MBEDGE:
	    dtk_type = ENTITY_TYPE_EDGE;
	    break;
	case moab::MBTRI:
	    dtk_type = ENTITY_TYPE_FACE;
	    break;
	case moab::MBQUAD:
	    dtk_type = ENTITY_TYPE_FACE;
	    break;
	case moab::MBPOLYGON:
	    dtk_type = ENTITY_TYPE_FACE;
	    break;
	case moab::MBTET:
	    dtk_type = ENTITY_TYPE_VOLUME;
	    break;
	case moab::MBPYRAMID:
	    dtk_type = ENTITY_TYPE_VOLUME;
	    break;
	case moab::MBPRISM:
	    dtk_type = ENTITY_TYPE_VOLUME;
	    break;
	case moab::MBKNIFE:
	    dtk_type = ENTITY_TYPE_VOLUME;
	    break;
	case moab::MBHEX:
	    dtk_type = ENTITY_TYPE_VOLUME;
	    break;
	case moab::MBPOLYHEDRON:
	    dtk_type = ENTITY_TYPE_VOLUME;
	    break;
	case moab::MBVERTEX:
	    dtk_type = ENTITY_TYPE_INVALID;
	    break;
	case moab::MBVERTEX:
	    dtk_type = ENTITY_TYPE_INVALID;
	    break;
	default:
	    dtk_type = ENTITY_TYPE_INVALID;
	    break;
    }
    return dtk_type;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_MoabHelpers.cpp
//---------------------------------------------------------------------------//
