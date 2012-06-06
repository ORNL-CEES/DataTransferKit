//---------------------------------------------------------------------------//
/*! 
 * \file tstMeshContainer.cpp
 * \author Stuart R. Slattery
 * \brief MeshContainer unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <set>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <DTK_MeshContainer.hpp>
#include <DTK_CoreTypes.hpp>
#include <DTK_MeshTraits.hpp>

#include <mpi.h>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_TypeTraits.hpp>

//---------------------------------------------------------------------------//
// MPI Setup
//---------------------------------------------------------------------------//

template<class Ordinal>
Teuchos::RCP<const Teuchos::Comm<Ordinal> > getDefaultComm()
{
#ifdef HAVE_MPI
    return Teuchos::DefaultComm<Ordinal>::getComm();
#else
    return Teuchos::rcp(new Teuchos::SerialComm<Ordinal>() );
#endif
}

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//

TEUCHOS_UNIT_TEST( MeshContainer, mesh_container_test )
{ 
    using namespace DataTransferKit;
 
    std::set<int> nodes;
    std::vector<double> coords;
    std::set<int> elements;
    std::vector<int> connectivity;

    MeshContainer<int> mesh_container( nodes, coords,
				       DTK_REGION, DTK_TETRAHEDRON, 4,
				       elements, connectivity );
}

TEUCHOS_UNIT_TEST( MeshContainer, serialization_test )
{
    using namespace DataTransferKit;

    int my_rank = getDefaultComm<int>()->getRank();

    std::set<int> nodes;
    std::vector<double> coords;
    std::set<int> elements;
    std::vector<int> connectivity;

    if ( my_rank == 0 )
    {
	nodes.insert( 1 );
	coords.push_back( 1.0 );
	elements.insert( 1 );
	connectivity.push_back( 1 );
    }

    MeshContainer<int> mesh_container( nodes, coords,
				       DTK_REGION, DTK_TETRAHEDRON, 4,
				       elements, connectivity );
    
    Teuchos::broadcast( *getDefaultComm<int>(), 0, &mesh_container );

    TEST_ASSERT( *mesh_container.nodesBegin() == 1 );
    TEST_ASSERT( *mesh_container.coordsBegin() == 1.0 );
    TEST_ASSERT( *mesh_container.elementsBegin() == 1 );
    TEST_ASSERT( *mesh_container.connectivityBegin() == 1 );
}

//---------------------------------------------------------------------------//
// end tstMeshContainer.cpp
//---------------------------------------------------------------------------//

