/****************************************************************************
 * Copyright (c) 2012-2017 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 ****************************************************************************/



#include <Kokkos_View.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_DefaultComm.hpp>

#include <Tpetra_Vector.hpp>
#include <Tpetra_Import.hpp>
#include <Tpetra_Export.hpp>

//---------------------------------------------------------------------------//
// Helper function to create uniquely owned maps.
template<class LO,class GO, class NO>
Teuchos::RCP<const Tpetra::Map<LO,GO,NO> >
createUniqueMap( const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                 const GO global_num_elements )
{
    // Create the map.
    return Tpetra::createUniformContigMapWithNode<LO,GO,NO>(
        global_num_elements, comm );
}

//---------------------------------------------------------------------------//
// Helper function to create ghosted maps.
template<class LO,class GO, class NO>
Teuchos::RCP<const Tpetra::Map<LO,GO,NO> >
createGhostedMap( const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                  const GO global_num_elements,
                  const int num_overlap )
{
    // First create a unique map to build from.
    auto unique_map = createUniqueMap<LO,GO,NO>( comm, global_num_elements );

    // Extract the contiguous global ids owned by this node.
    Teuchos::Array<GO> node_global_ids( unique_map->getNodeElementList() );

    // Check to see if this node has the first or last elements.
    GO front = node_global_ids.front();
    GO back = node_global_ids.back();
    bool has_first = ( 0 == front );
    bool has_last = ( global_num_elements -1 == back );

    // Add values for the overlap at the front.
    if ( !has_first )
        for ( int n = 0; n < num_overlap; ++n )
            node_global_ids.push_back( front-n-1 );

    // Add values for the overlap at the back.
    if ( !has_last )
        for ( int n = 0; n < num_overlap; ++n )
            node_global_ids.push_back( back+n+1 );

    // Create the overlapping map.
    return Tpetra::createNonContigMapWithNode<LO,GO,NO>(
        node_global_ids, comm );
}

//---------------------------------------------------------------------------//
// Unique distribution to unique distribution test.
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TpetraImportExport, unique_to_unique, Node )
{
    // Template aliases.
    using SC = double;
    using LO = int;
    using GO = unsigned long int;
    using NO = Node;

    // Get the communicator.
    auto comm = Teuchos::DefaultComm<int>::getComm();

    // Set the number of global elements.
    LO num_local = 100;
    GO num_global = num_local * comm->getSize();

    // Create maps.
    auto map_1 = createUniqueMap<LO,GO,NO>( comm, num_global );
    auto map_2 = createUniqueMap<LO,GO,NO>( comm, num_global );

    // Create vectors.
    auto vec_1 = Tpetra::createVector<SC,LO,GO,NO>( map_1 );
    auto vec_2 = Tpetra::createVector<SC,LO,GO,NO>( map_2 );

    // Create an export object.
    auto tpetra_export = Tpetra::createExport( map_1, map_2 );

    // Put some data in the vectors.
    SC value = 1.0;
    vec_1->putScalar( value );
    vec_2->putScalar( 0.0 );

    // Do the export.
    vec_2->doExport( *vec_1, *tpetra_export, Tpetra::ADD );

    // Check the 1-norm of the result.
    SC vec_norm_1 = vec_2->norm1();
    SC test_norm_1 = num_global * value;
    TEST_EQUALITY( vec_norm_1, test_norm_1 );

    // Check the infinity norm of the result.
    SC vec_norm_inf = vec_2->normInf();
    SC test_norm_inf = 1.0;
    TEST_EQUALITY( vec_norm_inf, test_norm_inf );
}

//---------------------------------------------------------------------------//
// Ghosted distribution to unique distribution test.
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TpetraImportExport, ghosted_to_unique, Node )
{
    // Template aliases.
    using SC = double;
    using LO = int;
    using GO = unsigned long int;
    using NO = Node;

    // Get the communicator.
    auto comm = Teuchos::DefaultComm<int>::getComm();

    // Set the number of global elements.
    LO num_local = 100;
    GO num_global = num_local * comm->getSize();

    // Create maps.
    int num_overlap = 4;
    auto map_1 = createGhostedMap<LO,GO,NO>( comm, num_global, num_overlap );
    auto map_2 = createUniqueMap<LO,GO,NO>( comm, num_global );

    // Create vectors.
    auto vec_1 = Tpetra::createVector<SC,LO,GO,NO>( map_1 );
    auto vec_2 = Tpetra::createVector<SC,LO,GO,NO>( map_2 );

    // Create an export object.
    auto tpetra_export = Tpetra::createExport( map_1, map_2 );

    // Put some data in the vectors.
    SC value = 1.0;
    vec_1->putScalar( value );
    vec_2->putScalar( 0.0 );

    // Do the export.
    vec_2->doExport( *vec_1, *tpetra_export, Tpetra::ADD );

    // Check the 1-norm of the result.
    SC vec_norm_1 = vec_2->norm1();
    SC test_norm_1 = num_global * value + (comm->getSize()-1) * 2 * num_overlap;
    TEST_EQUALITY( vec_norm_1, test_norm_1 );

    // Check the infinity norm of the result.
    SC vec_norm_inf = vec_2->normInf();
    SC test_norm_inf = (comm->getSize() > 1 ) ? 2.0 : 1.0;
    TEST_EQUALITY( vec_norm_inf, test_norm_inf );
}

//---------------------------------------------------------------------------//
// Unique distribution to ghosted distribution test.
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TpetraImportExport, unique_to_ghosted, Node )
{
    // Template aliases.
    using SC = double;
    using LO = int;
    using GO = unsigned long int;
    using NO = Node;

    // Get the communicator.
    auto comm = Teuchos::DefaultComm<int>::getComm();

    // Set the number of global elements.
    LO num_local = 100;
    GO num_global = num_local * comm->getSize();

    // Create maps.
    int num_overlap = 4;
    auto map_1 = createUniqueMap<LO,GO,NO>( comm, num_global );
    auto map_2 = createGhostedMap<LO,GO,NO>( comm, num_global, num_overlap );

    // Create vectors.
    auto vec_1 = Tpetra::createVector<SC,LO,GO,NO>( map_1 );
    auto vec_2 = Tpetra::createVector<SC,LO,GO,NO>( map_2 );

    // Create an import object.
    auto tpetra_import = Tpetra::createImport( map_1, map_2 );

    // Put some data in the vectors.
    SC value = 1.0;
    vec_1->putScalar( value );
    vec_2->putScalar( 0.0 );

    // Do the import.
    vec_2->doImport( *vec_1, *tpetra_import, Tpetra::ADD );

    // Check the 1-norm of the result.
    SC vec_norm_1 = vec_2->norm1();
    SC test_norm_1 = num_global * value + (comm->getSize()-1) * 2 * num_overlap;
    TEST_EQUALITY( vec_norm_1, test_norm_1 );

    // Check the infinity norm of the result.
    SC vec_norm_inf = vec_2->normInf();
    SC test_norm_inf = 1.0;
    TEST_EQUALITY( vec_norm_inf, test_norm_inf );
}

//---------------------------------------------------------------------------//
// Ghosted distribution to ghosted distribution test.
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TpetraImportExport, ghosted_to_ghosted, Node )
{
    // Template aliases.
    using SC = double;
    using LO = int;
    using GO = unsigned long int;
    using NO = Node;

    // Get the communicator.
    auto comm = Teuchos::DefaultComm<int>::getComm();

    // Set the number of global elements.
    LO num_local = 100;
    GO num_global = num_local * comm->getSize();

    // Create maps.
    int num_overlap_1 = 4;
    int num_overlap_2 = 8;
    auto map_1 = createGhostedMap<LO,GO,NO>( comm, num_global, num_overlap_1 );
    auto map_2 = createGhostedMap<LO,GO,NO>( comm, num_global, num_overlap_2 );
    auto unique_map_2 = Tpetra::createOneToOne( map_2 );

    // Create vectors.
    auto vec_1 = Tpetra::createVector<SC,LO,GO,NO>( map_1 );
    auto vec_2 = Tpetra::createVector<SC,LO,GO,NO>( map_2 );
    auto unique_vec_2 = Tpetra::createVector<SC,LO,GO,NO>( unique_map_2 );

    // Create an export object from map 1 to unique.
    auto tpetra_export_1_to_unique = Tpetra::createExport( map_1, unique_map_2 );

    // Create an import object from unique to map 2.
    auto tpetra_import_unique_to_2 = Tpetra::createImport( unique_map_2, map_2 );

    // Put some data in the vectors.
    SC value = 1.0;
    vec_1->putScalar( value );
    vec_2->putScalar( 0.0 );
    unique_vec_2->putScalar( 0.0 );

    // Export 1 to unique.
    unique_vec_2->doExport( *vec_1, *tpetra_export_1_to_unique, Tpetra::ADD );

    // Import unique to 2.
    vec_2->doImport( *unique_vec_2, *tpetra_import_unique_to_2, Tpetra::REPLACE );

    // Check the 1-norm of the result.
    SC vec_norm_1 = vec_2->norm1();
    SC test_norm_1 = num_global * value +
                     (comm->getSize()-1) * 4 * num_overlap_1 +
                     (comm->getSize()-1) * 2 * num_overlap_2;
    TEST_EQUALITY( vec_norm_1, test_norm_1 );

    // Check the infinity norm of the result.
    SC vec_norm_inf = vec_2->normInf();
    SC test_norm_inf = (comm->getSize() > 1 ) ? 2.0 : 1.0;
    TEST_EQUALITY( vec_norm_inf, test_norm_inf );
}

//---------------------------------------------------------------------------//

// Include the test macros.
#include "DataTransferKitDiscretization_ETIHelperMacros.h"

// Create the test group
#define UNIT_TEST_GROUP( NODE )                                                 \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TpetraImportExport, unique_to_unique, \
                                          NODE )                                \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TpetraImportExport, ghosted_to_unique, \
                                          NODE )                                \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TpetraImportExport, unique_to_ghosted, \
                                          NODE )                                \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TpetraImportExport, ghosted_to_ghosted, \
                                          NODE )

// Demangle the types
DTK_ETI_MANGLING_TYPEDEFS()

// Instantiate the tests
DTK_INSTANTIATE_N( UNIT_TEST_GROUP )
