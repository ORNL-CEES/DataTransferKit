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
 * \file DTK_MeshAdapter.hpp
 * \author Stuart R. Slattery
 * \brief FieldTraits adapter for mesh coordinates.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MESHADAPTER_HPP
#define DTK_MESHADAPTER_HPP

#include "DTK_FieldTraits.hpp"
#include <DTK_MeshTraits.hpp>
#include <DTK_MeshTools.hpp>

namespace DataTransferKit
{

//---------------------------------------------------------------------------//
/*! 
 * \class FieldTraits<Mesh>
 * \brief FieldTraits adapter for mesh coordinates. 
 *
 * This allows meshes to be viewed as a field of coordinates. Other mesh
 * information is not available through this interface.
 */
//---------------------------------------------------------------------------//
template<>
template<class MeshType>
class FieldTraits<MeshType>
{
    typedef MeshTraits<MeshType>                      MT;
    typedef MeshTools<MeshType>                       Tools;
    typedef MeshType                                  field_type;
    typedef double                                    value_type;
    typedef typename MT::global_ordinal_type          size_type;
    typedef typename MT::const_coordinate_iterator    iterator;
    typedef typename MT::const_coordinate_iterator    const_iterator;

    static inline std::size_t dim( const MeshType& mesh )
    { return MT::nodeDim( mesh ); }

    static inline std::size_t size( const MeshType& mesh )
    { return Tools::numNodes( mesh ) * MT::nodeDim( mesh ); }

    static inline bool empty( const MeshType& mesh )
    {
	if ( Tools::numNodes( mesh ) < 1 )
	{ 
	    return true;
	}
	else 
	{
	    return false;
	}
    }

    static inline iterator begin( MeshType& mesh )
    { return MT::coordsBegin( mesh ); }

    static inline const_iterator begin( const MeshType& mesh )
    { return MT::coordsBegin( mesh ); }

    static inline iterator end( MeshType& mesh )
    { return MT::coordsEnd( mesh ); }

    static inline const_iterator end( const MeshType& mesh )
    { return MT::coordsEnd( mesh ); }
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_MESHADAPTER

//---------------------------------------------------------------------------//
// end DTK_MeshAdapater.hpp
//---------------------------------------------------------------------------//

