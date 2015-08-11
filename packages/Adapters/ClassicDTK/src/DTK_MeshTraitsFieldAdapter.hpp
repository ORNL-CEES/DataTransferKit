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
 * \file DTK_MeshTraitsFieldAdapter.hpp
 * \author Stuart R. Slattery
 * \brief FieldTraits adapter for objects that have mesh traits.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MESHTRAITSFIELDADAPTER_HPP
#define DTK_MESHTRAITSFIELDADAPTER_HPP

#include "DTK_FieldTraits.hpp"
#include "DTK_MeshTraits.hpp"
#include "DTK_MeshTools.hpp"
#include "DTK_MeshContainer.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*! 
 * \brief FieldTraits adapter for mesh coordinates in MeshContainers.

 This allows the coordinate information stored in a block of mesh to be
 accessed as a field with field traits. Other mesh information is not
 available through this interface. The mesh object must have mesh traits.

 */
//---------------------------------------------------------------------------//
template<class GlobalOrdinal>
class FieldTraits< MeshContainer<GlobalOrdinal> >
{
  public:

    typedef MeshContainer<GlobalOrdinal>              MeshType;
    typedef MeshTraits<MeshType>                      MT;
    typedef MeshTools<MeshType>                       Tools;
    typedef MeshType                                  field_type;
    typedef double                                    value_type;
    typedef typename MT::global_ordinal_type          size_type;
    typedef typename MT::const_coordinate_iterator    iterator;
    typedef typename MT::const_coordinate_iterator    const_iterator;

    static inline int dim( const MeshType& mesh )
    { return MT::vertexDim( mesh ); }

    static inline size_type size( const MeshType& mesh )
    { return Tools::numVertices( mesh ) * MT::vertexDim( mesh ); }

    static inline bool empty( const MeshType& mesh )
    {
	if ( Tools::numVertices( mesh ) < 1 )
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

#endif // end DTK_MESHTRAITSFIELDADAPTER

//---------------------------------------------------------------------------//
// end DTK_MeshTraitsFieldAdapater.hpp
//---------------------------------------------------------------------------//

