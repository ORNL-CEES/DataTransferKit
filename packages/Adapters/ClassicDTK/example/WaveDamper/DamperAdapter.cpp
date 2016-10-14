//---------------------------------------------------------------------------//
/*!
 * \file DamperAdapter.cpp
 * \author Stuart R. Slattery
 * \brief Damper code adapters for DTK.
 */
//---------------------------------------------------------------------------//

#include <vector>

#include "DamperAdapter.hpp"

#include <DTK_MeshTypes.hpp>

//---------------------------------------------------------------------------//
// DAMPER ADAPTERS
//---------------------------------------------------------------------------//
/*!
 * \brief Get the damper mesh.
 */
Teuchos::RCP<DataTransferKit::MeshManager<DamperAdapter::MeshType>>
DamperAdapter::getMesh( const RCP_Damper &damper )
{
    // Get the process rank.
    int my_rank = damper->get_comm()->getRank();

    // Set the vertex dimension - this is a 1D problem.
    int vertex_dimension = 1;

    // Compute globally unique vertex ids.
    Teuchos::RCP<std::vector<double>> grid = damper->get_grid();
    Teuchos::ArrayRCP<int> vertices( grid->size() );
    for ( int i = 0; i < (int)vertices.size(); ++i )
    {
        vertices[i] = i + my_rank * vertices.size();
    }

    // Get the grid vertex coordinates.
    Teuchos::ArrayRCP<double> coordinates( &( *grid )[0], 0, grid->size(),
                                           false );

    // Set the grid topology - this is 1D so we are using line segments.
    DataTransferKit::DTK_ElementTopology element_topology =
        DataTransferKit::DTK_LINE_SEGMENT;

    // Each line segment will be constructed by 2 vertices.
    int vertices_per_element = 2;

    // Compute globally unique element ids.
    Teuchos::ArrayRCP<int> elements( grid->size() - 1 );
    for ( int i = 0; i < (int)elements.size(); ++i )
    {
        elements[i] = i + my_rank * elements.size();
    }

    // Generate element connectivity. The global vertex ids are used to
    // describe the construction of each line segment.
    Teuchos::ArrayRCP<int> connectivity( vertices_per_element *
                                         elements.size() );
    for ( int i = 0; i < (int)elements.size(); ++i )
    {
        connectivity[i] = vertices[i];
        connectivity[elements.size() + i] = vertices[i + 1];
    }

    // Define the permutation list. Here our line segments are ordered the
    // same as DTK canonical ordering so the list is an monotonically
    // increasing set of integers.
    Teuchos::ArrayRCP<int> permutation_list( vertices_per_element );
    for ( int i = 0; i < vertices_per_element; ++i )
    {
        permutation_list[i] = i;
    }

    // Build a DTK Mesh container with the data. The MeshType typedef is set
    // in the header file.
    Teuchos::RCP<MeshType> mesh_container = Teuchos::rcp( new MeshType(
        vertex_dimension, vertices, coordinates, element_topology,
        vertices_per_element, elements, connectivity, permutation_list ) );

    // We only have 1 element topology in this grid so we make just one mesh
    // block.
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType>> mesh_blocks( 1 );
    mesh_blocks[0] = mesh_container;

    // Return a mesh manager.
    return Teuchos::rcp( new DataTransferKit::MeshManager<MeshType>(
        mesh_blocks, damper->get_comm(), 1 ) );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the damper field evaluator.
 */
DamperAdapter::RCP_Evaluator
DamperAdapter::getFieldEvaluator( const RCP_Damper &damper )
{
    return Teuchos::rcp( new DamperEvaluator( damper ) );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the damper target coordinates directly from the mesh. We can do
 * this because we used the mesh container (See
 * DTK_MeshTraitsFieldAdapter.hpp).
 */
Teuchos::RCP<DataTransferKit::FieldManager<DamperAdapter::MeshType>>
DamperAdapter::getTargetCoords( const RCP_Damper &damper )
{
    return Teuchos::rcp( new DataTransferKit::FieldManager<MeshType>(
        getMesh( damper )->getBlock( 0 ), damper->get_comm() ) );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the damper target space.
 */
Teuchos::RCP<DataTransferKit::FieldManager<DamperAdapter::FieldType>>
DamperAdapter::getTargetSpace( const RCP_Damper &damper )
{
    // Get the damper data vector we will write into.
    Teuchos::RCP<std::vector<double>> target_space =
        damper->get_external_data();

    // The data we are transferring has 1 dimension.
    int field_dim = 1;

    // Build an ArrayRCP from the data vector.
    Teuchos::ArrayRCP<double> data_space( &( *target_space )[0], 0,
                                          target_space->size(), false );

    // Build a field container from the data space.
    Teuchos::RCP<FieldType> field_container =
        Teuchos::rcp( new FieldType( data_space, field_dim ) );

    // Return a field manager for the target data space.
    return Teuchos::rcp( new DataTransferKit::FieldManager<FieldType>(
        field_container, damper->get_comm() ) );
}

//---------------------------------------------------------------------------//
// end DamperAdapter.cpp
//---------------------------------------------------------------------------//
