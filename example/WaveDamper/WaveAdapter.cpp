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
// WAVE ADAPTERS
//---------------------------------------------------------------------------//
/*!
 * \brief Get the wave mesh.
 */
Teuchos::RCP<DataTransferKit::MeshManager<WaveAdapter::MeshType> >
WaveAdapter::getMesh( const RCP_Wave& wave )
{
    int my_rank = wave->get_comm()->getRank();
    int vertex_dimension = 1;

    Teuchos::RCP<std::vector<double> > grid = wave->get_grid();
    Teuchos::ArrayRCP<int> vertices( grid->size() );
    for ( int i = 0; i < (int) vertices.size(); ++i )
    {
	vertices[i] = i + my_rank*vertices.size();
    }
    Teuchos::ArrayRCP<double> coordinates( &(*grid)[0], 0, grid->size(), false );

    DataTransferKit::DTK_ElementTopology element_topology = 
	DataTransferKit::DTK_LINE_SEGMENT;
    int vertices_per_element = 2;

    Teuchos::ArrayRCP<int> elements( grid->size() - 1 );
    for ( int i = 0; i < (int) elements.size(); ++i )
    {
	elements[i] = i + my_rank*elements.size();
    }

    Teuchos::ArrayRCP<int> connectivity( vertices_per_element*elements.size() );
    for ( int i = 0; i < (int) elements.size(); ++i )
    {
	connectivity[i] = vertices[i];
	connectivity[ elements.size() + i ] = vertices[i+1];
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
			     mesh_blocks, wave->get_comm(), 1 ) );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the wave field evaluator.
 */
WaveAdapter::RCP_Evaluator 
WaveAdapter::getFieldEvaluator( const RCP_Wave& wave )
{
    return Teuchos::rcp( new WaveEvaluator( wave ) );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the wave target coordinates directly from the mesh. We can do
 * this because we used the mesh container (See
 * DTK_MeshTraitsFieldAdapter.hpp).
 */
Teuchos::RCP<DataTransferKit::FieldManager<WaveAdapter::MeshType> >
WaveAdapter::getTargetCoords( const RCP_Wave& wave )
{
    return Teuchos::rcp( 
	new DataTransferKit::FieldManager<MeshType>( 
	    getMesh( wave )->getBlock(0), wave->get_comm() ) );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the wave target space.
 */
Teuchos::RCP<DataTransferKit::FieldManager<WaveAdapter::FieldType> >
WaveAdapter::getTargetSpace( const RCP_Wave& wave )
{
    Teuchos::RCP<std::vector<double> > damping_space = wave->get_damping();
    int field_dim = 1;
    Teuchos::ArrayRCP<double> data_space( &(*damping_space)[0], 0, 
					  damping_space->size(), false );

    FieldType field_container( data_space, field_dim );
    
    return Teuchos::rcp( 
	new DataTransferKit::FieldManager<FieldType>( field_container, 
						      wave->get_comm() ) );
}

//---------------------------------------------------------------------------//
// end WaveAdapter.cpp
//---------------------------------------------------------------------------//
