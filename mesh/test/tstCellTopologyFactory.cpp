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
#include <Teuchos_DefaultMpiComm.hpp>
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

    test_topology = CellTopologyFactory::create( moab::MBEDGE, 2 );
    TEST_ASSERT( test_topology->getDimension() == 1 );
    TEST_ASSERT( test_topology->getNodeCount() == 2 );
    TEST_ASSERT( test_topology->getVertexCount() == 2 );

    test_topology = CellTopologyFactory::create( moab::MBEDGE, 3 );
    TEST_ASSERT( test_topology->getDimension() == 1 );
    TEST_ASSERT( test_topology->getNodeCount() == 3 );
    TEST_ASSERT( test_topology->getVertexCount() == 2 );

    test_topology = CellTopologyFactory::create( moab::MBEDGE, 100 );
    TEST_ASSERT( test_topology->getDimension() == 1 );
    TEST_ASSERT( test_topology->getNodeCount() == 2 );
    TEST_ASSERT( test_topology->getVertexCount() == 2 );

    test_topology = CellTopologyFactory::create( moab::MBTRI, 3 );
    TEST_ASSERT( test_topology->getDimension() == 2 );
    TEST_ASSERT( test_topology->getNodeCount() == 3 );
    TEST_ASSERT( test_topology->getVertexCount() == 3 );

    test_topology = CellTopologyFactory::create( moab::MBTRI, 4 );
    TEST_ASSERT( test_topology->getDimension() == 2 );
    TEST_ASSERT( test_topology->getNodeCount() == 4 );
    TEST_ASSERT( test_topology->getVertexCount() == 3 );

    test_topology = CellTopologyFactory::create( moab::MBTRI, 6 );
    TEST_ASSERT( test_topology->getDimension() == 2 );
    TEST_ASSERT( test_topology->getNodeCount() == 6 );
    TEST_ASSERT( test_topology->getVertexCount() == 3 );

    test_topology = CellTopologyFactory::create( moab::MBTRI, 100 );
    TEST_ASSERT( test_topology->getDimension() == 2 );
    TEST_ASSERT( test_topology->getNodeCount() == 3 );
    TEST_ASSERT( test_topology->getVertexCount() == 3 );

    test_topology = CellTopologyFactory::create( moab::MBQUAD, 4 );
    TEST_ASSERT( test_topology->getDimension() == 2 );
    TEST_ASSERT( test_topology->getNodeCount() == 4 );
    TEST_ASSERT( test_topology->getVertexCount() == 4 );

    test_topology = CellTopologyFactory::create( moab::MBQUAD, 8 );
    TEST_ASSERT( test_topology->getDimension() == 2 );
    TEST_ASSERT( test_topology->getNodeCount() == 8 );
    TEST_ASSERT( test_topology->getVertexCount() == 4 );

    test_topology = CellTopologyFactory::create( moab::MBQUAD, 9 );
    TEST_ASSERT( test_topology->getDimension() == 2 );
    TEST_ASSERT( test_topology->getNodeCount() == 9 );
    TEST_ASSERT( test_topology->getVertexCount() == 4 );

    test_topology = CellTopologyFactory::create( moab::MBQUAD, 100 );
    TEST_ASSERT( test_topology->getDimension() == 2 );
    TEST_ASSERT( test_topology->getNodeCount() == 4 );
    TEST_ASSERT( test_topology->getVertexCount() == 4 );

    test_topology = CellTopologyFactory::create( moab::MBTET, 4 );
    TEST_ASSERT( test_topology->getDimension() == 3 );
    TEST_ASSERT( test_topology->getNodeCount() == 4 );
    TEST_ASSERT( test_topology->getVertexCount() == 4 );

    test_topology = CellTopologyFactory::create( moab::MBTET, 8 );
    TEST_ASSERT( test_topology->getDimension() == 3 );
    TEST_ASSERT( test_topology->getNodeCount() == 8 );
    TEST_ASSERT( test_topology->getVertexCount() == 4 );

    test_topology = CellTopologyFactory::create( moab::MBTET, 10 );
    TEST_ASSERT( test_topology->getDimension() == 3 );
    TEST_ASSERT( test_topology->getNodeCount() == 10 );
    TEST_ASSERT( test_topology->getVertexCount() == 4 );
    
    test_topology = CellTopologyFactory::create( moab::MBTET, 11 );
    TEST_ASSERT( test_topology->getDimension() == 3 );
    TEST_ASSERT( test_topology->getNodeCount() == 11 );
    TEST_ASSERT( test_topology->getVertexCount() == 4 );

    test_topology = CellTopologyFactory::create( moab::MBTET, 100 );
    TEST_ASSERT( test_topology->getDimension() == 3 );
    TEST_ASSERT( test_topology->getNodeCount() == 4 );
    TEST_ASSERT( test_topology->getVertexCount() == 4 );

    test_topology = CellTopologyFactory::create( moab::MBPYRAMID, 5 );
    TEST_ASSERT( test_topology->getDimension() == 3 );
    TEST_ASSERT( test_topology->getNodeCount() == 5 );
    TEST_ASSERT( test_topology->getVertexCount() == 5 );

    test_topology = CellTopologyFactory::create( moab::MBPYRAMID, 13 );
    TEST_ASSERT( test_topology->getDimension() == 3 );
    TEST_ASSERT( test_topology->getNodeCount() == 5 );
    TEST_ASSERT( test_topology->getVertexCount() == 5 );

    test_topology = CellTopologyFactory::create( moab::MBPYRAMID, 14 );
    TEST_ASSERT( test_topology->getDimension() == 3 );
    TEST_ASSERT( test_topology->getNodeCount() == 14 );
    TEST_ASSERT( test_topology->getVertexCount() == 5 );

    test_topology = CellTopologyFactory::create( moab::MBPYRAMID, 100 );
    TEST_ASSERT( test_topology->getDimension() == 3 );
    TEST_ASSERT( test_topology->getNodeCount() == 5 );
    TEST_ASSERT( test_topology->getVertexCount() == 5 );

    test_topology = CellTopologyFactory::create( moab::MBPRISM, 6 );
    TEST_ASSERT( test_topology->getDimension() == 3 );
    TEST_ASSERT( test_topology->getNodeCount() == 6 );
    TEST_ASSERT( test_topology->getVertexCount() == 6 );

    test_topology = CellTopologyFactory::create( moab::MBPRISM, 15 );
    TEST_ASSERT( test_topology->getDimension() == 3 );
    TEST_ASSERT( test_topology->getNodeCount() == 6 );
    TEST_ASSERT( test_topology->getVertexCount() == 6 );

    test_topology = CellTopologyFactory::create( moab::MBPRISM, 18 );
    TEST_ASSERT( test_topology->getDimension() == 3 );
    TEST_ASSERT( test_topology->getNodeCount() == 18 );
    TEST_ASSERT( test_topology->getVertexCount() == 6 );

    test_topology = CellTopologyFactory::create( moab::MBPRISM, 100 );
    TEST_ASSERT( test_topology->getDimension() == 3 );
    TEST_ASSERT( test_topology->getNodeCount() == 6 );
    TEST_ASSERT( test_topology->getVertexCount() == 6 );

    test_topology = CellTopologyFactory::create( moab::MBHEX, 8 );
    TEST_ASSERT( test_topology->getDimension() == 3 );
    TEST_ASSERT( test_topology->getNodeCount() == 8 );
    TEST_ASSERT( test_topology->getVertexCount() == 8 );

    test_topology = CellTopologyFactory::create( moab::MBHEX, 20 );
    TEST_ASSERT( test_topology->getDimension() == 3 );
    TEST_ASSERT( test_topology->getNodeCount() == 8 );
    TEST_ASSERT( test_topology->getVertexCount() == 8 );

    test_topology = CellTopologyFactory::create( moab::MBHEX, 27 );
    TEST_ASSERT( test_topology->getDimension() == 3 );
    TEST_ASSERT( test_topology->getNodeCount() == 27 );
    TEST_ASSERT( test_topology->getVertexCount() == 8 );

    test_topology = CellTopologyFactory::create( moab::MBHEX, 100 );
    TEST_ASSERT( test_topology->getDimension() == 3 );
    TEST_ASSERT( test_topology->getNodeCount() == 8 );
    TEST_ASSERT( test_topology->getVertexCount() == 8 );
}

//---------------------------------------------------------------------------//
// end tstCellTopologyFactory.cpp
//---------------------------------------------------------------------------//

