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
Teuchos::RCP<DataTransferKit::MeshManager<DamperAdapter::MeshType> >
DamperAdapter::getMesh( const RCP_Damper& damper )
{
    int my_rank = damper->get_comm()->getRank();

    int vertex_dimension = 1;

    std::vector<double> grid = damper->get_grid();
    Teuchos::ArrayRCP<int> vertices( grid.size() );
    for ( int i = 0; i < (int) vertices.size(); ++i )
    {
	vertices[i] = i + my_rank*vertices.size();
    }
    Teuchos::ArrayRCP<double> coordinates( &grid[0], 0, grid.size(), false );

    DataTransferKit::DTK_ElementTopology element_topology = 
	DataTransferKit::DTK_LINE_SEGMENT;
    int vertices_per_element = 2;

    Teuchos::ArrayRCP<int> elements( grid.size() - 1 );
    for ( int i = 0; i < (int) elements.size(); ++i )
    {
	elements[i] = i + my_rank*elements.size();
    }

    Teuchos::ArrayRCP<int> connectivity( vertices_per_element*elements.size() );
    for ( int i = 0; i < (int) elements.size(); ++i )
    {
	connectivity[i] = i;
	connectivity[ elements.size() + i ] = i+1;
    }

    Teuchos::ArrayRCP<int> permutation_list( vertices_per_element );
    for ( int i = 0; i < vertices_per_element; ++i )
    {
	permutation_list[i] = i;
    }

    MeshType mesh_container( vertex_dimension, vertices, coordinates, 
			     element_topology, vertices_per_element, 
			     elements, connectivity, permutation_list );

    Teuchos::ArrayRCP<MeshType> mesh_blocks( 1 );
    mesh_blocks[0] = mesh_container;

    return Teuchos::rcp( new DataTransferKit::MeshManager<MeshType>( 
			     mesh_blocks, damper->get_comm(), 1 ) );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the damper field evaluator.
 */
DamperAdapter::RCP_Evaluator 
DamperAdapter::getFieldEvaluator( const RCP_Damper& damper )
{
    return Teuchos::rcp( new DamperEvaluator( damper ) );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the damper target coordinates directly from the mesh. We can do
 * this because we used the mesh container (See
 * DTK_MeshTraitsFieldAdapter.hpp).
 */
Teuchos::RCP<DataTransferKit::FieldManager<DamperAdapter::MeshType> >
DamperAdapter::getTargetCoords( const RCP_Damper& damper )
{
    return Teuchos::rcp( 
	new DataTransferKit::FieldManager<MeshType>( 
	    getMesh( damper )->getBlock(0), damper->get_comm() ) );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the damper target space.
 */
Teuchos::RCP<DataTransferKit::FieldManager<DamperAdapter::FieldType> >
DamperAdapter::getTargetSpace( const RCP_Damper& damper )
{
    std::vector<double> external_space = damper->get_external_data();
    int field_dim = 1;
    Teuchos::ArrayRCP<double> data_space( &external_space[0], 0, 
					  external_space.size(), false );

    FieldType field_container( data_space, field_dim );
    
    return Teuchos::rcp(
	new DataTransferKit::FieldManager<FieldType>( field_container, 
						      damper->get_comm() ) );
}

//---------------------------------------------------------------------------//
// end DamperAdapter.cpp
//---------------------------------------------------------------------------//
