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
 * \brief DTK_GeometryTraits.hpp
 * \author Stuart R. Slattery
 * \brief Declaration of geometry traits.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_GEOMETRYTRAITS_HPP
#define DTK_GEOMETRYTRAITS_HPP

#include "DTK_BoundingBox.hpp"

#include <Teuchos_Array.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Dummy struct. If a type does not create a specialization this will
 * not compile.
 */
template<typename UndefinedGeometryType>
struct UndefinedGeometryTraits
{
    static inline UndefinedGeometryType notDefined() 
    { return UndefinedGeometryType::this_type_is_missing_a_specialization(); }
};

//---------------------------------------------------------------------------//
/*!
  \class GeometryTraits
  \brief Geometry traits definitions.

  GeometryTraits provide access to basic properites of geometric objects.
*/
//---------------------------------------------------------------------------//
template<typename GeometryType>
class GeometryTraits
{
  public:

    //! Typedef for geometry type.
    typedef GeometryType geometry_type;

    /*!
     * \brief Return the dimension of the geometry.
     */
    static inline int dim( const GeometryType& geometry )
    { UndefinedGeometryTraits<GeometryType>::notDefined(); return 0; }

    /*!
     * \brief Return the geometry measure (volume for a 3D geometry, area for
     * 2D, and length for 1D).
     */
    static inline double measure( const GeometryType& geometry )
    { UndefinedGeometryTraits<GeometryType>::notDefined(); return 0; }

    /*!
     * \brief Return whether or not a point is in the geometry within a
     * specified tolerance.
     */
    static inline bool pointInGeometry( const GeometryType& geometry,
					const Teuchos::Array<double>& coords,
					const double tolerance )
    { UndefinedGeometryTraits<GeometryType>::notDefined(); return 0; }

    /*!
     * \brief Return the axis-aligned bounding box around the geometry.
     */
    static inline BoundingBox boundingBox( const GeometryType& geometry )
    { UndefinedGeometryTraits<GeometryType>::notDefined(); return 0; }
};

} // end namespace DataTransferKit

#endif // end DTK_GEOMETRYTRAITS_HPP

//---------------------------------------------------------------------------//
// end DTK_GeometryTraits.hpp
//---------------------------------------------------------------------------//
