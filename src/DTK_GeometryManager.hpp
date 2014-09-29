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

#include "DTK_GeometricEntity.hpp"
#include "DTK_BoundingBox.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayView.hpp>

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
class GeometryManager
{
  public:

    // Constructor.
    GeometryManager( const Teuchos::ArrayRCP<GeometricEntity>& geometry,
		     const Teuchos::ArrayRCP<MeshId>& geom_gids,
		     const Teuchos::RCP<const Teuchos::Comm<int> >& comm );

    // Destructor.
    ~GeometryManager();

    //! Get the local geometric objects managed by this manager.
    const Teuchos::ArrayRCP<GeometricEntity>& geometry() const 
    { return d_geometry; }

    //! Get the global ids of the geometric objects managed by this manager.
    const Teuchos::ArrayRCP<MeshId>& gids() const 
    { return d_geom_gids; }

    //! Get the communicator for the geometry.
    const Teuchos::RCP<const Teuchos::Comm<int> >& comm() const
    { return d_comm; }

    //! Get the local number of objects owned by this manager.
    const typename Teuchos::ArrayRCP<GeometricEntity>::size_type
    localNumGeometry() const 
    { return d_geometry.size(); }

    // Get the global number of objects owned by this manager.
    const typename Teuchos::ArrayRCP<GeometricEntity>::size_type
    globalNumGeometry() const;

    // Get the bounding boxes for the objects owned by this manager.
    Teuchos::Array<BoundingBox> boundingBoxes() const;

    // Get the local bounding box for the objects owned by this manager.
    BoundingBox localBoundingBox() const;

    // Get the global bounding box for the objects owned by this manager.
    BoundingBox globalBoundingBox() const;

    //! Set the active geometries.
    void setActiveGeometry( const Teuchos::Array<short int>& active_geometry )
    { d_active_geometry = active_geometry; }

    //! Get the active geometry.
    Teuchos::ArrayView<short int> getActiveGeometry()
    { return d_active_geometry(); }

  private:

    // Geometric objects.
    Teuchos::ArrayRCP<GeometricEntity> d_geometry;

    // Geometric objects global ids.
    Teuchos::ArrayRCP<MeshId> d_geom_gids;
    
    // Communicator over which the geometry is defined.
    Teuchos::RCP<const Teuchos::Comm<int> > d_comm;

    // Active geometry.
    Teuchos::Array<short int> d_active_geometry;
};

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // DTK_GEOMETRYMANAGER_HPP

//---------------------------------------------------------------------------//
// end DTK_GeometryManager.hpp
//---------------------------------------------------------------------------//

