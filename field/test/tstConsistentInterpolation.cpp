//---------------------------------------------------------------------------//
/*!
 * \file tstConsistentInterpolation.cpp
 * \author Stuart R. Slattery
 * \brief Consistent interpolation unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <DTK_ConsistentInterpolation.hpp>
#include <DTK_FieldTraits.hpp>
#include <DTK_FieldEvaluator.hpp>
#include <DTK_MeshTypes.hpp>
#include <DTK_MeshTraits.hpp>

#include <mpi.h>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_Array.hpp>
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
// Mesh Implementation
//---------------------------------------------------------------------------//

class MyMesh
{
  public:

    typedef int    global_ordinal_type;
    
    MyMesh() 
    { /* ... */ }

    MyMesh( const Teuchos::Array<int>& node_handles,
	    const Teuchos::Array<double>& coords,
	    const Teuchos::Array<int>& quad_handles,
	    const Teuchos::Array<int>& quad_connectivity )
	: d_node_handles( node_handles )
	, d_coords( coords )
	, d_quad_handles( quad_handles )
	, d_quad_connectivity( quad_connectivity )
    { /* ... */ }

    ~MyMesh()
    { /* ... */ }

    Teuchos::Array<int>::const_iterator nodesBegin() const
    { return d_node_handles.begin(); }

    Teuchos::Array<int>::const_iterator nodesEnd() const
    { return d_node_handles.end(); }

    Teuchos::Array<double>::const_iterator coordsBegin() const
    { return d_coords.begin(); }

    Teuchos::Array<double>::const_iterator coordsEnd() const
    { return d_coords.end(); }

    Teuchos::Array<int>::const_iterator quadsBegin() const
    { return d_quad_handles.begin(); }

    Teuchos::Array<int>::const_iterator quadsEnd() const
    { return d_quad_handles.end(); }

    Teuchos::Array<int>::const_iterator connectivityBegin() const
    { return d_quad_connectivity.begin(); }

    Teuchos::Array<int>::const_iterator connectivityEnd() const
    { return d_quad_connectivity.end(); }
    

  private:

    Teuchos::Array<int> d_node_handles;
    Teuchos::Array<double> d_coords;
    Teuchos::Array<int> d_quad_handles;
    Teuchos::Array<int> d_quad_connectivity;
};

//---------------------------------------------------------------------------//
// DTK Traits Specializations
//---------------------------------------------------------------------------//
namespace DataTransferKit
{

//---------------------------------------------------------------------------//
// Mesh traits specialization for MyMesh
template<>
class MeshTraits<MyMesh>
{
  public:

    typedef MyMesh::global_ordinal_type global_ordinal_type;
    typedef Teuchos::Array<int>::const_iterator const_node_iterator;
    typedef Teuchos::Array<double>::const_iterator const_coordinate_iterator;
    typedef Teuchos::Array<int>::const_iterator const_element_iterator;
    typedef Teuchos::Array<int>::const_iterator const_connectivity_iterator;

    static inline std::size_t nodeDim( const MyMesh& mesh )
    { return 2; }

    static inline const_node_iterator nodesBegin( const MyMesh& mesh )
    { return mesh.nodesBegin(); }

    static inline const_node_iterator nodesEnd( const MyMesh& mesh )
    { return mesh.nodesEnd(); }

    static inline const_coordinate_iterator coordsBegin( const MyMesh& mesh )
    { return mesh.coordsBegin(); }

    static inline const_coordinate_iterator coordsEnd( const MyMesh& mesh )
    { return mesh.coordsEnd(); }


    static inline std::size_t elementType( const MyMesh& mesh )
    { return DTK_FACE; }

    static inline std::size_t elementTopology( const MyMesh& mesh )
    { return DTK_QUADRILATERAL; }

    static inline std::size_t nodesPerElement( const MyMesh& mesh )
    { return 4; }


    static inline const_element_iterator elementsBegin( const MyMesh& mesh )
    { return mesh.quadsBegin(); }

    static inline const_element_iterator elementsEnd( const MyMesh& mesh )
    { return mesh.quadsEnd(); }

    static inline const_connectivity_iterator connectivityBegin( const MyMesh& mesh )
    { return mesh.connectivityBegin(); }

    static inline const_connectivity_iterator connectivityEnd( const MyMesh& mesh )
    { return mesh.connectivityEnd(); }
};

//---------------------------------------------------------------------------//
// Field Traits specification for Teuchos::Array
template<>
class FieldTraits< Teuchos::Array<double> >
{
  public:

    typedef Teuchos::Array<double>                    field_type;
    typedef double                                    value_type;
    typedef Teuchos::Array<double>::size_type         size_type;
    typedef Teuchos::Array<double>::iterator          iterator;
    typedef Teuchos::Array<double>::const_iterator    const_iterator;

    static inline size_type dim( const Teuchos::Array<double>& field )
    { return 1; }

    static inline size_type size( const Teuchos::Array<double>& field )
    { return field.size(); }

    static inline bool empty( const Teuchos::Array<double>& field )
    { return field.empty(); }

    static inline iterator begin( Teuchos::Array<double>& field )
    { return field.begin(); }

    static inline const_iterator begin( Teuchos::Array<double>& field ) const
    { return field.begin(); }

    static inline iterator end( Teuchos::Array<double>& field )
    { return field.end(); }

    static inline const_iterator end( Teuchos::Array<double>& field ) const
    { return field.end(); }
};

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Mesh create funciton.
//---------------------------------------------------------------------------//
MyMesh buildMyMesh()
{
    int my_rank = getDefaultComm<int>()->getRank();

    // Make some nodes.
    int num_nodes = 10;
    int node_dim = 2;
    Teuchos::Array<int> node_handles( num_nodes );
    Teuchos::Array<double> coords( node_dim*num_nodes );

    for ( int i = 0; i < num_nodes; ++i )
    {
	node_handles[i] = (num_nodes / 2)*my_rank + i;
    }
    for ( int i = 0; i < num_nodes / 2; ++i )
    {
	coords[ i ] = my_rank;
	coords[ num_nodes + i ] = i;
    }
    for ( int i = num_nodes / 2; i < num_nodes; ++i )
    {
	coords[ i ] = my_rank + 1;
	coords[ num_nodes + i ] = i - num_nodes/2;
    }
    
    // Make the quads.
    int num_quads = 4;
    Teuchos::Array<int> quad_handles( num_quads );
    Teuchos::Array<int> quad_connectivity( 4*num_quads );
    
    for ( int i = 0; i < num_quads; ++i )
    {
	quad_handles[ i ] = num_quads*my_rank + i;
	quad_connectivity[ i ] = node_handles[i];
	quad_connectivity[ num_quads + i ] = node_handles[num_nodes/2 + i];
	quad_connectivity[ 2*num_quads + i ] = node_handles[num_nodes/2 + i + 1];
	quad_connectivity[ 3*num_quads + i ] = node_handles[i+1];
    }

    return MyMesh( node_handles, coords, quad_handles, quad_connectivity );
}

//---------------------------------------------------------------------------//
// Unit tests
//---------------------------------------------------------------------------//

TEUCHOS_UNIT_TEST( ConsistentInterpolation, consistent_interpolation_test )
{
    // Setup communication.
    Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm<int>();
    int my_rank = comm->getRank();
    int my_size = comm->getSize();

    // This is a 4 processor test.
    if ( my_size == 4 )
    {
	// Setup source mesh.
	MyMesh source_mesh = buildMyMesh();

	// Setup target coordinate field.
	Teuchos::Array<double> target_coords = buildCoordinateField();

	// Create field evaluator.
	Teuchos::RCP< FieldEvaluator<MyMesh,Teuchos::Array<double> > 
		      my_evaluator = Teuchos::rcp( new MyEvaluator() );

	// Create data target.
	Teuchos::Array<double> my_target;

	// Create a map for a consistent interpolation scheme.
	ConsistentInterpolation map( comm );

	// Setup and apply the interpolation to the field.
	ConsistentInterpolation<MyMesh,Teuchos::Array<double> > 
		      transfer_operator( map );
	transfer_operator.setup( source_mesh, target_coords );
	transfer_operator.apply( my_evaluator, my_target );
    }
}

//---------------------------------------------------------------------------//
// end tstTransferOperator.cpp
//---------------------------------------------------------------------------//
