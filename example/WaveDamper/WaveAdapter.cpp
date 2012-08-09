//---------------------------------------------------------------------------//
/*!
 * \file WaveAdapter.cpp
 * \author Stuart R. Slattery
 * \brief Wave code adapters for DTK.
 */
//---------------------------------------------------------------------------//

#include <vector>

#include "WaveAdapter.hpp"

#include <DTK_MeshTypes.hpp>

//---------------------------------------------------------------------------//
/*!
 * \brief Get the wave mesh.
 */
DataTransferKit::MeshManager<Wave::MeshType>
Wave::getMesh( const RCP_Wave& wave )
{
    int vertex_dimension = 1;

    std::vector<double> grid = wave->get_grid();
    Teuchos::ArrayRCP<int> vertices( grid.size() );
    for ( int i = 0; i < (int) vertices.size(); ++i )
    {
	vertices[i] = i;
    }
    Teuchos::ArrayRCP<double> coordinates( &grid[0], 0, grid.size(), false );

    DataTransferKit::DTK_ELementTopology element_topology = DTK_LINE_SEGMENT;
    int vertices_per_element = 2;

    Teuchos::ArrayRCP<int> elements( grid.size() - 1 );
    for ( int i = 0; i < (int) elements.size(); ++i )
    {
	elements[i] = i;
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

    MeshType mesh_container( vertex_dimension, coordinates, element_topology,
			     vertices_per_element, elements, connectivity,
			     permutation_list );

    Teuchos::ArrayRCP<MeshType> mesh_blocks( 1 );
    mesh_blocks[0] = mesh_container;

    return DataTransferKit::MeshManager<MeshType>( 
	mesh_container, wave->get_comm(), 1 );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the wave field evaluator.
 */
Wave::RCP_Evaluator 
Wave::getFieldEvaluator( const RCP_Wave& wave )
{
    return Teuchos::rcp( new WaveEvaluator( wave ) );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the wave target coordinates.
 */
DataTransferKit::FieldManager<Wave::MT> 
Wave::getTargetCoords( const RCP_Wave& wave )
{
    return DataTransferKit::FieldManager<MT>( getMesh( wave ), wave->get_comm );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the wave target space.
 */
DataTransferKit::FieldManager<Wave::FieldType> 
Wave::getTargetSpace( const RCP_Wave& wave )
{
    std::vector<double> damping_space = wave->get_damping();
    int field_dim = 1;
    Teuchos::ArrayRCP<double> data_space( &damping_space, 0, 
					  damping_space.size(), false );

    FieldType field_container( data_space, field_dim );
    
    return DataTransferKit::FieldManager<FieldType>( field_container, 
						     wave->get_comm() );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Function evaluator.
 */
DataTransferKit::FieldContainer<double> 
WaveEvaluator::evaluate( const Teuchos::ArrayRCP<int>& elements,
			 const Teuchos::ArrayRCP<double>& coords )
{

}

//---------------------------------------------------------------------------//
// end WaveAdapter.cpp
//---------------------------------------------------------------------------//


