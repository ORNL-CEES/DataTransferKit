//---------------------------------------------------------------------------//
/*!
 * \file cxx_main.cpp
 * \author Stuart R. Slattery
 * \brief Integral assembly transfer example.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <algorithm>
#include <cassert>

#include <DTK_IntegralAssemblyMap.hpp>
#include <DTK_MeshManager.hpp>
#include <DTK_MeshContainer.hpp>
#include <DTK_MeshTypes.hpp>
#include <DTK_MeshTools.hpp>
#include <DTK_ElementMeasure.hpp>
#include <DTK_FieldIntegrator.hpp>
#include <DTK_FieldManager.hpp>
#include <DTK_FieldContainer.hpp>
#include <DTK_FieldTools.hpp>
#include <DTK_CommTools.hpp>
#include <DTK_GeometryManager.hpp>
#include <DTK_Box.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_Ptr.hpp>

//---------------------------------------------------------------------------//
// Source mesh creation.
//---------------------------------------------------------------------------//
Teuchos::RCP<DataTransferKit::MeshContainer<int> >
createSourceMesh( const Teuchos::RCP<const Teuchos::Comm<int> >& comm )
{
    int src_rank = comm->getRank();

    // Create the vertex ids.
    int num_verts = 42;
    Teuchos::ArrayRCP<int> vertex_ids( num_verts );
    for ( int i = 0; i < num_verts; ++i )
    {
        vertex_ids[i] = src_rank*21 + i;
    }

    // Create the vertex coordinates. Blocked and split in parallel in Z
    // direction.
    int num_coords = num_verts*3;
    Teuchos::ArrayRCP<double> vertex_coords( num_coords );
    for ( int k = src_rank; k < src_rank+2; ++k ) {
        for ( int j = 0; j < 3; ++j ) {
            for ( int i = 0; i < 7; ++i ) {
                int index = i + 7*j + 21*(k-src_rank);
                vertex_coords[index] = 0.5*i;
                vertex_coords[index+num_verts] = 0.5*j;
                vertex_coords[index+2*num_verts] = 0.5*k;
            }
        }
    }

    // Create the element ids.
    int num_elements = 12;
    Teuchos::ArrayRCP<int> element_ids( num_elements );
    for ( int k = 0; k < 1; ++k ) {
        for ( int j = 0; j < 2; ++j ) {
            for ( int i = 0; i < 6; ++i ) {
                int elem_index = i + 6*j + 12*k;
                element_ids[elem_index] = src_rank*12 + elem_index;
            }
        }
    }

    // Create the element connectivity from the vertex ids.
    int conn_size = 8;
    Teuchos::ArrayRCP<int> element_conn( num_elements*conn_size );
    for ( int k = 0; k < 1; ++k ) {
        for ( int j = 0; j < 2; ++j ) {
            for ( int i = 0; i < 6; ++i ) {
                int vert_index = i + 7*j + 13*k;
                int elem_index = i + 6*j + 12*k;

                element_conn[elem_index] = vertex_ids[vert_index];
                element_conn[elem_index+num_elements] = vertex_ids[vert_index+1];
                element_conn[elem_index+2*num_elements] = vertex_ids[vert_index+8];
                element_conn[elem_index+3*num_elements] = vertex_ids[vert_index+7];
                element_conn[elem_index+4*num_elements] = vertex_ids[vert_index+21];
                element_conn[elem_index+5*num_elements] = vertex_ids[vert_index+22];
                element_conn[elem_index+6*num_elements] = vertex_ids[vert_index+29];
                element_conn[elem_index+7*num_elements] = vertex_ids[vert_index+28];
            }
        }
    }

    // Create the permutation list.
    Teuchos::ArrayRCP<int> permutation( conn_size );
    for ( int i = 0; i < conn_size; ++i ) permutation[i] = i;

    // Create the mesh container.
    return Teuchos::rcp( new DataTransferKit::MeshContainer<int>(
                             3,
                             vertex_ids,
                             vertex_coords,
                             DataTransferKit::DTK_HEXAHEDRON,
                             conn_size,
                             element_ids,
                             element_conn,
                             permutation ) );
}

//---------------------------------------------------------------------------//
// Target geometry creation.
//---------------------------------------------------------------------------//
Teuchos::RCP<DataTransferKit::GeometryManager<DataTransferKit::Box,int> >
createTargetGeometry( const Teuchos::RCP<const Teuchos::Comm<int> >& comm )
{
    int tgt_rank = comm->getRank();

    // Allocate space for 2 boxes.
    Teuchos::ArrayRCP<DataTransferKit::Box> boxes( 2 );
    Teuchos::ArrayRCP<int> box_ids( 2 );

    // Create 2 boxes. 1 box will be shared with the other proc.
    boxes[0] = DataTransferKit::Box( tgt_rank, 0, 0,
                                     tgt_rank+1, 1, 1 );
    box_ids[0] = tgt_rank;

    boxes[1] = DataTransferKit::Box( tgt_rank+1, 0, 0,
                                     tgt_rank+2, 1, 1 );
    box_ids[1] = tgt_rank + 1;

    return Teuchos::rcp(
        new DataTransferKit::GeometryManager<DataTransferKit::Box,int>(
            boxes, box_ids, comm, 3 ) );
}

//---------------------------------------------------------------------------//
// Source function integrator.
//---------------------------------------------------------------------------//
class SourceIntegrator:
    public DataTransferKit::FieldIntegrator<DataTransferKit::MeshContainer<int>,
                                            DataTransferKit::FieldContainer<double> >
{
  public:

    SourceIntegrator( const Teuchos::ArrayRCP<int>& element_gids )
        : d_element_gids( element_gids )
    {
        std::sort( element_gids.begin(), element_gids.end() );
    }

    ~SourceIntegrator()
    { /* ... */ }

    DataTransferKit::FieldContainer<double> integrate(
        const Teuchos::ArrayRCP<int>& gids )
    {
        std::cout << gids() << std::endl;
        Teuchos::ArrayRCP<double> integrated_data( gids.size() );
        for ( int n = 0; n < gids.size(); ++n )
        {
            if ( gids[n] == 0  || gids[n] == 1  || gids[n] == 6  || gids[n] == 7 ||
                 gids[n] == 12 || gids[n] == 13 || gids[n] == 18 || gids[n] == 19 )
            {
                integrated_data[n] = 1.0;
            }
            else if ( gids[n] == 2  || gids[n] == 3  || gids[n] == 8  || gids[n] == 9 ||
                      gids[n] == 14 || gids[n] == 15 || gids[n] == 20 || gids[n] == 21 )
            {
                integrated_data[n] = 2.0;
            }
            else if ( gids[n] == 4  || gids[n] == 5  || gids[n] == 10  || gids[n] == 11 ||
                      gids[n] == 16 || gids[n] == 17 || gids[n] == 22  || gids[n] == 23 )
            {
                integrated_data[n] = 3.0;
            }
            else
            {
                integrated_data[n] = 0.0;
            }
        }

        return DataTransferKit::FieldContainer<double>( integrated_data, 1 );
    }

  private:

    Teuchos::ArrayRCP<int> d_element_gids;
};

//---------------------------------------------------------------------------//
// Source mesh element measure.
//---------------------------------------------------------------------------//
class SourceMeasure:
    public DataTransferKit::ElementMeasure<DataTransferKit::MeshContainer<int> >
{
  public:

    SourceMeasure( const Teuchos::ArrayRCP<int>& element_gids )
        : d_element_gids( element_gids )
    {
        std::sort( element_gids.begin(), element_gids.end() );
    }

    ~SourceMeasure()
    { /* ... */ }

    Teuchos::Array<double> measure(
        const Teuchos::ArrayRCP<int>& gids )
    {
        Teuchos::Array<double> measure_data( gids.size() );
        for ( int n = 0; n < gids.size(); ++n )
        {
            if ( std::binary_search( d_element_gids.begin(), d_element_gids.end(),
                                     gids[n] ) )
            {
                measure_data[n] = 1.0;
            }
            else
            {
                measure_data[n] = 0.0;
            }
        }

        return measure_data;
    }

  private:

    Teuchos::ArrayRCP<int> d_element_gids;
};

//---------------------------------------------------------------------------//
// Main function driver for the integral assembly transfer problem.
int main(int argc, char* argv[])
{
    // ---------------//
    // PARALLEL SETUP
    // ---------------//

    // Setup communication.
    Teuchos::GlobalMPISession mpiSession(&argc,&argv);
    Teuchos::RCP<const Teuchos::Comm<int> > comm_default =
        Teuchos::DefaultComm<int>::getComm();
    int num_procs = comm_default->getSize();

    // Only 4 procs for this example.
    if( num_procs == 4 )
    {
        // Split the main communicator into 2 separate groups, one for the source
        // geometry and one for the target mesh.
        Teuchos::Array<int> sub_ranks_src(2), sub_ranks_tgt(2);
        for ( int i = 0; i < 2; ++i )
        {
            sub_ranks_src[i] = i;
            sub_ranks_tgt[i] = i+2;
        }

        // Generate the source and target communicators from the sub ranks.
        Teuchos::RCP<const Teuchos::Comm<int> > src_comm =
            comm_default->createSubcommunicator( sub_ranks_src() );

        Teuchos::RCP<const Teuchos::Comm<int> > tgt_comm =
            comm_default->createSubcommunicator( sub_ranks_tgt() );

        // Build the union communicator for the source and target. This is the
        // communicator over which we will operate the coupling.
        Teuchos::RCP<const Teuchos::Comm<int> > comm_union;
        DataTransferKit::CommTools::unite( src_comm, tgt_comm, comm_union );

        // Set a boolean for source/target existence.
        bool src_exists = false;
        if ( !src_comm.is_null() ) src_exists = true;
        bool tgt_exists = false;
        if ( !tgt_comm.is_null() ) tgt_exists = true;
        comm_union->barrier();



        // ---------------//
        // SOURCE SETUP
        // ---------------//

        // Set required variables in the scope of the global communicator.
        Teuchos::RCP<DataTransferKit::MeshManager<DataTransferKit::MeshContainer<int> > >
            src_mesh;
        Teuchos::RCP<DataTransferKit::ElementMeasure<DataTransferKit::MeshContainer<int> > >
            src_measure;
        Teuchos::RCP<DataTransferKit::FieldIntegrator<DataTransferKit::MeshContainer<int>,
                                                      DataTransferKit::FieldContainer<double> > >
            src_integrator;

        // If the source code exists on this process, build its data structures.
        if ( src_exists )
        {

            // Get the source mesh.
            Teuchos::ArrayRCP<Teuchos::RCP<DataTransferKit::MeshContainer<int> > > mesh_blocks(1);
            mesh_blocks[0] = createSourceMesh( src_comm );
            src_mesh = Teuchos::rcp(
                new DataTransferKit::MeshManager<DataTransferKit::MeshContainer<int> >(
                    mesh_blocks, src_comm, 3 ) );

            // Get the element measure.
            Teuchos::ArrayRCP<int> mesh_element_ids =
                DataTransferKit::MeshTools<DataTransferKit::MeshContainer<int> >::elementsNonConstView( *mesh_blocks[0] );
            src_measure = Teuchos::rcp( new SourceMeasure( mesh_element_ids ) );

            // Get the source field integrator.
            src_integrator = Teuchos::rcp( new SourceIntegrator( mesh_element_ids ) );
        }
        comm_union->barrier();



        // ---------------//
        // TARGET SETUP
        // ---------------//

        // Set required variables in the scope of the global communicator.
        Teuchos::RCP<DataTransferKit::GeometryManager<DataTransferKit::Box,int> >
            tgt_geometry;
        Teuchos::ArrayRCP<double> tgt_array;
        Teuchos::RCP<
            DataTransferKit::FieldManager<DataTransferKit::FieldContainer<double> > >
            tgt_space;

        // If the target code exists on this process, build its data structures.
        if ( tgt_exists )
        {
            // Get target geometry.
            tgt_geometry = createTargetGeometry( tgt_comm );

            // Get the target data space.
            tgt_array = Teuchos::ArrayRCP<double>( 2 );
            Teuchos::RCP<DataTransferKit::FieldContainer<double> > tgt_container =
                Teuchos::rcp( new DataTransferKit::FieldContainer<double>( tgt_array, 1 ) );
            tgt_space = Teuchos::rcp(
                new DataTransferKit::FieldManager<DataTransferKit::FieldContainer<double> >(
                    tgt_container, tgt_comm ) );
        }
        comm_union->barrier();



        // ---------------//
        // MAPPING SETUP
        // ---------------//

        // Create the mapping for the source-to-target transfer. The mapping
        // will occur over the union communicator in 3 dimensions. Keep track
        // of missed points as we expect to miss some.
        DataTransferKit::IntegralAssemblyMap<DataTransferKit::MeshContainer<int>,
                                             DataTransferKit::Box>
            src_to_tgt_map( comm_union, 3 );

        // Setup the source-to-target map with the source geometry as the source
        // and the target coordinates as the target.
        src_to_tgt_map.setup( src_mesh, src_measure, tgt_geometry );



        // ---------------//
        // DATA TRANSFER
        // ---------------//

        // Apply the mapping to transfer the data.
        src_to_tgt_map.apply( src_integrator, tgt_space );

        // Check the resulting data transfer.
        int fail_count = 0;
        if ( tgt_exists )
        {
            int target_rank = tgt_comm->getRank();
            if ( target_rank == 0 )
            {
                if ( tgt_array[0] != 1.0 ) ++fail_count;
                if ( tgt_array[1] != 2.0 ) ++fail_count;
            }
            else if ( target_rank == 1 )
            {
                if ( tgt_array[0] != 2.0 ) ++fail_count;
                if ( tgt_array[1] != 3.0 ) ++fail_count;
            }

            std::cout << std::endl;
            if ( fail_count == 0 )
            {
                std::cout << "TEST PASSED: target proc "
                          << tgt_comm->getRank() << std::endl;
            }
            else
            {
                std::cout << "TEST FAILED " << fail_count << ": target proc "
                          << tgt_comm->getRank() << std::endl;
            }
            std::cout << std::endl;
        }
        comm_union->barrier();

    }
    return 0;
}

//---------------------------------------------------------------------------//
// end cxx_main.hpp
//---------------------------------------------------------------------------//

