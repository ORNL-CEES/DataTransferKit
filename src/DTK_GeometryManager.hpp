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
 * \file DTK_GeometryManager.hpp
 * \author Stuart R. Slattery
 * \brief Geometry manager declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_GEOMETRYMANAGER_HPP
#define DTK_GEOMETRYMANAGER_HPP

#include "DTK_GeometryTraits.hpp"
#include "DTK_BoundingBox.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_ArrayRCP.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 \class GeometryManager
 \brief Manager object for geometry.

 The geometry manager manages a collection of geometric objects and their
 parallel decomposition. A geometry has a dimension of arbitrary size. For
 example, a square has a dimension of 2 while a cylinder has a dimension of
 3. A collection of geometric objects need not know their parallel
 decomposition, but they exist on a shared parallel communicator. Individual
 geometric objects in the collection may not span a parallel process in
 pieces. If a geometric object exists on the domain owned by many processors,
 it must exist on those processors. All geometric objects managed by the
 manager may not intersect. This is true for both local and global objects.

 */
//---------------------------------------------------------------------------//
template<class Geometry,class GlobalOrdinal>
class GeometryManager
{
  public:

    //@{
    //! Typedefs.
    typedef Geometry                                      geometry_type;
    typedef GeometryTraits<Geometry>                      GT;
    typedef GlobalOrdinal                                 global_ordinal_type;
    typedef Teuchos::Comm<int>                            CommType;
    typedef Teuchos::RCP<const CommType>                  RCP_Comm;
    //@}

    // Constructor.
    GeometryManager( const Teuchos::ArrayRCP<Geometry>& geometry,
		     const Teuchos::ArrayRCP<GlobalOrdinal>& geom_gids,
		     const RCP_Comm& comm, const int dim );

    // Destructor.
    ~GeometryManager();

    //! Get the geometric objects managed by this manager.
    const Teuchos::ArrayRCP<Geometry>& geometry() const 
    { return d_geometry; }

    //! Get the global ids of the geometric objects managed by this manager.
    const Teuchos::ArrayRCP<GlobalOrdinal>& gids() const 
    { return d_geom_gids; }

    //! Get the communicator for the geometry.
    const RCP_Comm& comm() const
    { return d_comm; }

    //! Get the dimension of the geometry.
    const int dim() const
    { return d_dim; }

    //! Get the local number of objects owned by this manager.
    const typename Teuchos::ArrayRCP<Geometry>::size_type
    localNumGeometry() const 
    { return d_geometry.size(); }

    // Get the global number of objects owned by this manager.
    const typename Teuchos::ArrayRCP<Geometry>::size_type
    globalNumGeometry() const;

    // Get the bounding boxes for the objects owned by this manager.
    Teuchos::Array<BoundingBox> boundingBoxes() const;

    // Get the local bounding box for the objects owned by this manager.
    BoundingBox localBoundingBox() const;

    // Get the global bounding box for the objects owned by this manager.
    BoundingBox globalBoundingBox() const;

  private:

    // Validate the geometric objects to the domain model.
    void validate();

  private:

    // Geometric objects.
    Teuchos::ArrayRCP<Geometry> d_geometry;

    // Geometric objects global ids.
    Teuchos::ArrayRCP<GlobalOrdinal> d_geom_gids;
    
    // Communicator over which the geometry is defined.
    RCP_Comm d_comm;

    // The dimension of the geometry.
    int d_dim;
};

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_GeometryManager_def.hpp"

//---------------------------------------------------------------------------//

#endif // DTK_GEOMETRYMANAGER_HPP

//---------------------------------------------------------------------------//
// end DTK_GeometryManager.hpp
//---------------------------------------------------------------------------//

