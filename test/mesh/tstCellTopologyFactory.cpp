//---------------------------------------------------------------------------//
/*!
 * \file tstCellTopologyFactory.cpp
 * \author Stuart R. Slattery
 * \brief Cell topology factory unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <DTK_CellTopologyFactory.hpp>

#include <mpi.h>

#include <MBInterface.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_TypeTraits.hpp>

#include <Shards_CellTopology.hpp>

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
// Unit tests.
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( CellTopologyFactory, cell_topology_factory_test )
{
    using namespace DataTransferKit;
    typedef Teuchos::RCP<shards::CellTopology>   RCP_CellTopology;

    RCP_CellTopology test_topology;

    test_topology = CellTopologyFactory::create( DTK_LINE_SEGMENT, 2 );
    TEST_ASSERT( test_topology->getDimension() == 1 );
    TEST_ASSERT( test_topology->getNodeCount() == 2 );
    TEST_ASSERT( test_topology->getVertexCount() == 2 );

    test_topology = CellTopologyFactory::create( DTK_LINE_SEGMENT, 3 );
    TEST_ASSERT( test_topology->getDimension() == 1 );
    TEST_ASSERT( test_topology->getNodeCount() == 3 );
    TEST_ASSERT( test_topology->getVertexCount() == 2 );

    test_topology = CellTopologyFactory::create( DTK_LINE_SEGMENT, 100 );
    TEST_ASSERT( test_topology->getDimension() == 1 );
    TEST_ASSERT( test_topology->getNodeCount() == 2 );
    TEST_ASSERT( test_topology->getVertexCount() == 2 );

    test_topology = CellTopologyFactory::create( DTK_TRIANGLE, 3 );
    TEST_ASSERT( test_topology->getDimension() == 2 );
    TEST_ASSERT( test_topology->getNodeCount() == 3 );
    TEST_ASSERT( test_topology->getVertexCount() == 3 );

    test_topology = CellTopologyFactory::create( DTK_TRIANGLE, 4 );
    TEST_ASSERT( test_topology->getDimension() == 2 );
    TEST_ASSERT( test_topology->getNodeCount() == 4 );
    TEST_ASSERT( test_topology->getVertexCount() == 3 );

    test_topology = CellTopologyFactory::create( DTK_TRIANGLE, 6 );
    TEST_ASSERT( test_topology->getDimension() == 2 );
    TEST_ASSERT( test_topology->getNodeCount() == 6 );
    TEST_ASSERT( test_topology->getVertexCount() == 3 );

    test_topology = CellTopologyFactory::create( DTK_TRIANGLE, 100 );
    TEST_ASSERT( test_topology->getDimension() == 2 );
    TEST_ASSERT( test_topology->getNodeCount() == 3 );
    TEST_ASSERT( test_topology->getVertexCount() == 3 );

    test_topology = CellTopologyFactory::create( DTK_QUADRILATERAL, 4 );
    TEST_ASSERT( test_topology->getDimension() == 2 );
    TEST_ASSERT( test_topology->getNodeCount() == 4 );
    TEST_ASSERT( test_topology->getVertexCount() == 4 );

    test_topology = CellTopologyFactory::create( DTK_QUADRILATERAL, 8 );
    TEST_ASSERT( test_topology->getDimension() == 2 );
    TEST_ASSERT( test_topology->getNodeCount() == 8 );
    TEST_ASSERT( test_topology->getVertexCount() == 4 );

    test_topology = CellTopologyFactory::create( DTK_QUADRILATERAL, 9 );
    TEST_ASSERT( test_topology->getDimension() == 2 );
    TEST_ASSERT( test_topology->getNodeCount() == 9 );
    TEST_ASSERT( test_topology->getVertexCount() == 4 );

    test_topology = CellTopologyFactory::create( DTK_QUADRILATERAL, 100 );
    TEST_ASSERT( test_topology->getDimension() == 2 );
    TEST_ASSERT( test_topology->getNodeCount() == 4 );
    TEST_ASSERT( test_topology->getVertexCount() == 4 );

    test_topology = CellTopologyFactory::create( DTK_TETRAHEDRON, 4 );
    TEST_ASSERT( test_topology->getDimension() == 3 );
    TEST_ASSERT( test_topology->getNodeCount() == 4 );
    TEST_ASSERT( test_topology->getVertexCount() == 4 );

    test_topology = CellTopologyFactory::create( DTK_TETRAHEDRON, 8 );
    TEST_ASSERT( test_topology->getDimension() == 3 );
    TEST_ASSERT( test_topology->getNodeCount() == 8 );
    TEST_ASSERT( test_topology->getVertexCount() == 4 );

    test_topology = CellTopologyFactory::create( DTK_TETRAHEDRON, 10 );
    TEST_ASSERT( test_topology->getDimension() == 3 );
    TEST_ASSERT( test_topology->getNodeCount() == 10 );
    TEST_ASSERT( test_topology->getVertexCount() == 4 );
    
    test_topology = CellTopologyFactory::create( DTK_TETRAHEDRON, 11 );
    TEST_ASSERT( test_topology->getDimension() == 3 );
    TEST_ASSERT( test_topology->getNodeCount() == 11 );
    TEST_ASSERT( test_topology->getVertexCount() == 4 );

    test_topology = CellTopologyFactory::create( DTK_TETRAHEDRON, 100 );
    TEST_ASSERT( test_topology->getDimension() == 3 );
    TEST_ASSERT( test_topology->getNodeCount() == 4 );
    TEST_ASSERT( test_topology->getVertexCount() == 4 );

    test_topology = CellTopologyFactory::create( DTK_PYRAMID, 5 );
    TEST_ASSERT( test_topology->getDimension() == 3 );
    TEST_ASSERT( test_topology->getNodeCount() == 5 );
    TEST_ASSERT( test_topology->getVertexCount() == 5 );

    test_topology = CellTopologyFactory::create( DTK_PYRAMID, 13 );
    TEST_ASSERT( test_topology->getDimension() == 3 );
    TEST_ASSERT( test_topology->getNodeCount() == 5 );
    TEST_ASSERT( test_topology->getVertexCount() == 5 );

    test_topology = CellTopologyFactory::create( DTK_PYRAMID, 14 );
    TEST_ASSERT( test_topology->getDimension() == 3 );
    TEST_ASSERT( test_topology->getNodeCount() == 14 );
    TEST_ASSERT( test_topology->getVertexCount() == 5 );

    test_topology = CellTopologyFactory::create( DTK_PYRAMID, 100 );
    TEST_ASSERT( test_topology->getDimension() == 3 );
    TEST_ASSERT( test_topology->getNodeCount() == 5 );
    TEST_ASSERT( test_topology->getVertexCount() == 5 );

    test_topology = CellTopologyFactory::create( DTK_WEDGE, 6 );
    TEST_ASSERT( test_topology->getDimension() == 3 );
    TEST_ASSERT( test_topology->getNodeCount() == 6 );
    TEST_ASSERT( test_topology->getVertexCount() == 6 );

    test_topology = CellTopologyFactory::create( DTK_WEDGE, 15 );
    TEST_ASSERT( test_topology->getDimension() == 3 );
    TEST_ASSERT( test_topology->getNodeCount() == 6 );
    TEST_ASSERT( test_topology->getVertexCount() == 6 );

    test_topology = CellTopologyFactory::create( DTK_WEDGE, 18 );
    TEST_ASSERT( test_topology->getDimension() == 3 );
    TEST_ASSERT( test_topology->getNodeCount() == 18 );
    TEST_ASSERT( test_topology->getVertexCount() == 6 );

    test_topology = CellTopologyFactory::create( DTK_WEDGE, 100 );
    TEST_ASSERT( test_topology->getDimension() == 3 );
    TEST_ASSERT( test_topology->getNodeCount() == 6 );
    TEST_ASSERT( test_topology->getVertexCount() == 6 );

    test_topology = CellTopologyFactory::create( DTK_HEXAHEDRON, 8 );
    TEST_ASSERT( test_topology->getDimension() == 3 );
    TEST_ASSERT( test_topology->getNodeCount() == 8 );
    TEST_ASSERT( test_topology->getVertexCount() == 8 );

    test_topology = CellTopologyFactory::create( DTK_HEXAHEDRON, 20 );
    TEST_ASSERT( test_topology->getDimension() == 3 );
    TEST_ASSERT( test_topology->getNodeCount() == 8 );
    TEST_ASSERT( test_topology->getVertexCount() == 8 );

    test_topology = CellTopologyFactory::create( DTK_HEXAHEDRON, 27 );
    TEST_ASSERT( test_topology->getDimension() == 3 );
    TEST_ASSERT( test_topology->getNodeCount() == 27 );
    TEST_ASSERT( test_topology->getVertexCount() == 8 );

    test_topology = CellTopologyFactory::create( DTK_HEXAHEDRON, 100 );
    TEST_ASSERT( test_topology->getDimension() == 3 );
    TEST_ASSERT( test_topology->getNodeCount() == 8 );
    TEST_ASSERT( test_topology->getVertexCount() == 8 );
}

//---------------------------------------------------------------------------//
// end tstCellTopologyFactory.cpp
//---------------------------------------------------------------------------//

