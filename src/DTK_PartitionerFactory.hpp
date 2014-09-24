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
 * \file DTK_PartitionerFactory.hpp
 * \author Stuart R. Slattery
 * \brief PartitionerFactory interface definition.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_PARTITIONERFACTORY_HPP
#define DTK_PARTITIONERFACTORY_HPP

#include "DTK_Partitioner.hpp"
#include "DTK_MeshManager.hpp"
#include "DTK_GeometryManager.hpp"
#include "DTK_SerialPartitioner.hpp"

#ifdef HAVE_DTK_MPI
#include "DTK_RCB.hpp"
#include "DTK_GeometryRCB.hpp"
#endif

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>

namespace DataTransferKit
{

//---------------------------------------------------------------------------//
/*!
 * \class PartitionerFactory
 * \brief Factory for generating partitioners.
 */
//---------------------------------------------------------------------------//
class PartitionerFactory
{
  public:
    
    // Constructor.
    PartitionerFactory()
    { /* ... */ }

    // Destructor.
    ~PartitionerFactory()
    { /* ... */ }

    // Mesh factory method.
    static inline Teuchos::RCP<Partitioner> 
    createMeshPartitioner( 
	const Teuchos::RCP<const Teuchos::Comm<int> > comm,
	const Teuchos::RCP<MeshManager> mesh_manager,
	const int dimension );

    // Geometry factory method.
    template<class Geometry, class GlobalOrdinal>
    static inline Teuchos::RCP<Partitioner> 
    createGeometryPartitioner( 
	const Teuchos::RCP<const Teuchos::Comm<int> > comm,
	const Teuchos::RCP<GeometryManager<Geometry,GlobalOrdinal> > geometry_manager,
	const int dimension );
};

//---------------------------------------------------------------------------//
// Inline functions.
//---------------------------------------------------------------------------//
Teuchos::RCP<Partitioner> PartitionerFactory::createMeshPartitioner(
    const Teuchos::RCP<const Teuchos::Comm<int> > comm,
    const Teuchos::RCP<MeshManager> mesh_manager,
    const int dimension )
{
#ifdef HAVE_DTK_MPI
    return Teuchos::rcp( new RCB( comm, mesh_manager, dimension ) );
#else
    return Teuchos::rcp( new SerialPartitioner() );
#endif
}

//---------------------------------------------------------------------------//
template<class Geometry, class GlobalOrdinal>
Teuchos::RCP<Partitioner> PartitionerFactory::createGeometryPartitioner(
    const Teuchos::RCP<const Teuchos::Comm<int> > comm,
    const Teuchos::RCP<GeometryManager<Geometry,GlobalOrdinal> > geometry_manager,
    const int dimension )
{
#ifdef HAVE_DTK_MPI
    return Teuchos::rcp( new GeometryRCB<Geometry,GlobalOrdinal>( 
			     comm, geometry_manager, dimension ) );
#else
    return Teuchos::rcp( new SerialPartitioner() );
#endif
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_PARTITIONERFACTORY_HPP

//---------------------------------------------------------------------------//
// end DTK_PartitionerFactory.hpp
//---------------------------------------------------------------------------//
