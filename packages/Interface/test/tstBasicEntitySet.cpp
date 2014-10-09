//---------------------------------------------------------------------------//
/*!
 * \file tstBasicEntitySet.cpp
 * \author Stuart R. Slattery
 * \brief BasicEntitySet unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <DTK_BasicEntitySet.hpp>
#include <DTK_Entity.hpp>
#include <DTK_AbstractObjectRegistry.hpp>
#include <DTK_DataSerializer.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_TypeTraits.hpp>
#include <Teuchos_Tuple.hpp>

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
// NodeImpl implementation for testing.
//---------------------------------------------------------------------------//
class NodeImpl : public DataTransferKit::EntityImpl
{
  public:

    // Default constructor.
    NodeImpl() 
	: d_global_id( DataTransferKit::dtk_invalid_entity_id )
	, d_owner_rank( -1 )
	, d_coordinates( 0 )
    { /* ... */ }

    // Array constructor.
    NodeImpl( const DataTransferKit::EntityId global_id, 
	      const int owner_rank,
	      const Teuchos::Array<double>& coordinates )
	: d_global_id( global_id )
	, d_owner_rank( owner_rank )
	, d_coordinates( coordinates )
    { /* ... */ }

    // Destructor.
    ~NodeImpl()
    { /* ... */ }

    // Return a string indicating the derived entity type.
    std::string name() const
    {
	return std::string("Unit Test Node");
    }

    // Get the entity type.
    DataTransferKit::EntityType entityType() const
    {
	return DataTransferKit::NODE;
    }

    // Get the unique global identifier for the entity.
    DataTransferKit::EntityId id() const
    {
	return d_global_id;
    }
    
    // Get the parallel rank that owns the entity.
    int ownerRank() const
    {
	return d_owner_rank;
    }

    // Return the physical dimension of the entity.
    int physicalDimension() const
    {
	return d_coordinates.size();
    }

    // Return the parametric dimension of the entity.
    int parametricDimension() const
    {
	return 0;
    }

    // Return the entity measure with respect to the parameteric
    double measure() const
    {
	return 0.0;
    }

    // Return the centroid of the entity.
    void centroid( Teuchos::ArrayView<const double>& centroid ) const
    {
	centroid = d_coordinates();
    }

    // Return the axis-aligned bounding box around the entity.
    void boundingBox( Teuchos::Tuple<double,6>& bounds ) const
    {
	bounds = Teuchos::tuple( d_coordinates[0], d_coordinates[1], d_coordinates[2], 
				 d_coordinates[0], d_coordinates[1], d_coordinates[2] );
    }

    // Perform a safeguard check for mapping a point to the reference
    void safeguardMapToReferenceFrame(
	const Teuchos::ParameterList& parameters,
	const Teuchos::ArrayView<const double>& point,
	DataTransferKit::MappingStatus& status ) const
    { /* ... */ }

    // Map a point to the reference space of an entity. Return the
    void mapToReferenceFrame( 
	const Teuchos::ParameterList& parameters,
	const Teuchos::ArrayView<const double>& point,
	const Teuchos::ArrayView<double>& reference_point,
	DataTransferKit::MappingStatus& status ) const
    { /* ... */ }

    // Determine if a reference point is in the parameterized space of
    bool checkPointInclusion( 
	const Teuchos::ParameterList& parameters,
	const Teuchos::ArrayView<const double>& reference_point ) const
    { return false; }

    // Map a reference point to the physical space of an entity.
    void mapToPhysicalFrame( 
	const Teuchos::ArrayView<const double>& reference_point,
	const Teuchos::ArrayView<double>& point ) const
    { /* ... */ }
     
    // Serialize the entity into a buffer.
    void serialize( const Teuchos::ArrayView<char>& buffer ) const
    {
	DataTransferKit::DataSerializer serializer;
	serializer.setBuffer( buffer );
	serializer << d_global_id << d_owner_rank;

	Teuchos::Array<double>::const_iterator coord_it;
	for ( coord_it = d_coordinates.begin(); 
	      coord_it != d_coordinates.end();
	      ++coord_it )
	{
	    serializer << *coord_it;
	}
    }

    // Deserialize an entity from a buffer.
    void deserialize( const Teuchos::ArrayView<const char>& buffer )
    {
	Teuchos::ArrayView<char> buffer_nonconst(
	    const_cast<char*>(buffer.getRawPtr()), buffer.size() );

	DataTransferKit::DataDeserializer deserializer;
	deserializer.setBuffer( buffer_nonconst );
	deserializer >> d_global_id >> d_owner_rank;

	d_coordinates.resize( 3 );
	Teuchos::Array<double>::iterator coord_it;
	for ( coord_it = d_coordinates.begin(); 
	      coord_it != d_coordinates.end();
	      ++coord_it )
	{
	    deserializer >> *coord_it;
	}
    }

    // Get the byte size for the box.
    static std::size_t byteSize();

  private:

    // Global id.
    DataTransferKit::EntityId d_global_id;

    // Owning parallel rank.
    int d_owner_rank;

    // Coordinates.
    Teuchos::Array<double> d_coordinates;

    // Packed size in bytes.
    static std::size_t d_byte_size;
};

// Byte size of the point.
std::size_t 
NodeImpl::d_byte_size = 
    sizeof(DataTransferKit::EntityId) + sizeof(int) + 3*sizeof(double);

//---------------------------------------------------------------------------//
// Get the byte size of the point.
std::size_t NodeImpl::byteSize()
{
    return d_byte_size;
}

//---------------------------------------------------------------------------//
// Node implementation for testing.
//---------------------------------------------------------------------------//
class Node : public DataTransferKit::Entity
{
  public:

    // Default constructor.
    Node()
    {
	this->b_entity_impl = Teuchos::rcp( new NodeImpl() );
    }


    // Array constructor.
    Node( const DataTransferKit::EntityId global_id, 
	  const int owner_rank,
	  const Teuchos::Array<double>& coordinates )
    {
	this->b_entity_impl =
	    Teuchos::rcp( new NodeImpl(global_id,owner_rank,coordinates) );
    }

    // Destructor.
    ~Node() { /* ... */ }

    // Get the byte size for the node.
    static std::size_t byteSize();
};

std::size_t Node::byteSize()
{
    return NodeImpl::byteSize();
}

//---------------------------------------------------------------------------//
namespace DataTransferKit
{
template<>
class AbstractObjectRegistrationPolicy<Node>
{
  public:

    //! Base class type.
    typedef Node object_type;

    /*!
     * \brief Register a derived class with a base class.
     */
    static void registerDerivedClassWithBaseClass()
    {
	// Register the constructor with the base class
	// AbstractBuildableObject interface.
	Entity::setDerivedClassFactory<Node>();

	// Register the byte size with the base class
	// AbstractSerializableObject interface.
	Entity::setDerivedClassByteSize( Node::byteSize() );
    }
};
} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
// Local point set test.
TEUCHOS_UNIT_TEST( Point, set_test )
{
    using namespace DataTransferKit;

    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > comm = 
	Teuchos::DefaultComm<int>::getComm();
    int comm_rank = comm->getRank();
    int comm_size = comm->getSize();

    // Make point.
    double x_1 = 3.2 + comm_rank;
    double y_1 = -9.233 + comm_rank;
    double z_1 = 1.3 + comm_rank;
    Teuchos::Array<double> p1(3);
    p1[0] = x_1;
    p1[1] = y_1;
    p1[2] = z_1;
    Entity point_1 = Node(0, comm_rank, p1);

    // Make a second point.
    double x_2 = 3.2 - comm_rank;
    double y_2 = -9.233 - comm_rank;
    double z_2 = 1.3 - comm_rank;
    Teuchos::Array<double> p2(3);
    p2[0] = x_2;
    p2[1] = y_2;
    p2[2] = z_2;
    Entity point_2 = Node(1, comm_rank, p2);

    // Make an entity set.
    Teuchos::RCP<EntitySet> entity_set = Teuchos::rcp(
	new BasicEntitySet(comm, 3) );

    // Add the points to the set.
    entity_set->addEntity( point_1 );
    entity_set->addEntity( point_2 );

    // Get an iterator to the entity set objects.
    AbstractIterator<Entity> node_it =
	entity_set->entityIterator( NODE );
    AbstractIterator<Entity> edge_it =
	entity_set->entityIterator( EDGE );
    AbstractIterator<Entity> face_it =
	entity_set->entityIterator( FACE );
    AbstractIterator<Entity> volume_it =
	entity_set->entityIterator( VOLUME );

    // Check the entity set.
    TEST_EQUALITY( entity_set->name(), "DTK Basic Entity Set" );
    TEST_EQUALITY( entity_set->physicalDimension(), 3 );
    TEST_EQUALITY( node_it.size(), 2 );
    TEST_EQUALITY( edge_it.size(), 0 );
    TEST_EQUALITY( face_it.size(), 0 );
    TEST_EQUALITY( volume_it.size(), 0 );

    // Check the nodes.
    node_it = node_it.begin();
    Entity ge0 = *node_it;
    ++node_it;
    Entity ge1 = *node_it;
    if ( ge1.id() == 1 )
    {
	TEST_EQUALITY( ge0.id(), 0 );
    }
    else if ( ge1.id() == 0 )
    {
	TEST_EQUALITY( ge0.id(), 1 );
    }
    else
    {
	TEST_ASSERT( false );
    }
    Entity entity;
    entity_set->getEntity( 0, entity );
    TEST_EQUALITY( 0, entity.id() );
    entity_set->getEntity( 1, entity );
    TEST_EQUALITY( 1, entity.id() );

    // Check the bounding boxes.
    Teuchos::Tuple<double,6> local_bounds;
    entity_set->localBoundingBox( local_bounds );
    TEST_EQUALITY( local_bounds[0], x_2 );
    TEST_EQUALITY( local_bounds[1], y_2 );
    TEST_EQUALITY( local_bounds[2], z_2 );
    TEST_EQUALITY( local_bounds[3], x_1 );
    TEST_EQUALITY( local_bounds[4], y_1 );
    TEST_EQUALITY( local_bounds[5], z_1 );

    Teuchos::Tuple<double,6> global_bounds;
    entity_set->globalBoundingBox( global_bounds );
    TEST_FLOATING_EQUALITY( global_bounds[0], 3.2 - comm_size + 1.0, 1.0e-12 );
    TEST_FLOATING_EQUALITY( global_bounds[1], -9.233 - comm_size + 1.0, 1.0e-12 );
    TEST_FLOATING_EQUALITY( global_bounds[2], 1.3 - comm_size + 1.0, 1.0e-12 );
    TEST_FLOATING_EQUALITY( global_bounds[3], 3.2 + comm_size - 1.0, 1.0e-12 );
    TEST_FLOATING_EQUALITY( global_bounds[4], -9.233 + comm_size - 1.0, 1.0e-12 );
    TEST_FLOATING_EQUALITY( global_bounds[5], 1.3 + comm_size - 1.0, 1.0e-12 );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( Point, modification_test )
{
    using namespace DataTransferKit;

    // Register the Entity classes.
    AbstractObjectRegistry<Entity,Node >::registerDerivedClasses();

    // Register the EntitySet classes.
    AbstractObjectRegistry<EntitySet,BasicEntitySet>::registerDerivedClasses();

    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > comm = 
	Teuchos::DefaultComm<int>::getComm();
    int comm_size = comm->getSize();

    // Create a builder for the entity sets.
    Teuchos::RCP<AbstractBuilder<EntitySet> > builder = EntitySet::getBuilder();

    // Make an entity set on process 0.
    Teuchos::RCP<EntitySet> entity_set;
    Teuchos::Array<Entity> points(2);
    int entity_set_key = -1;
    double x_1 = 3.2 + comm_size;
    double y_1 = -9.233 + comm_size;
    double z_1 = 1.3 + comm_size;
    double x_2 = 3.2 - comm_size;
    double y_2 = -9.233 - comm_size;
    double z_2 = 1.3 - comm_size;
    if ( 0 == comm->getRank() )
    {
	Teuchos::Array<double> p1(3);
	p1[0] = x_1;
	p1[1] = y_1;
	p1[2] = z_1;
	points[0] = Node(0, 0, p1);
	Teuchos::Array<double> p2(3);
	p2[0] = x_2;
	p2[1] = y_2;
	p2[2] = z_2;
	points[1] = Node(1, 0, p2);

	entity_set = Teuchos::rcp(new BasicEntitySet(comm,3) );
	entity_set_key = builder->getIntegralKey( entity_set->name() );
    }

    // Create an entity set.
    Teuchos::broadcast( *comm, 0, Teuchos::Ptr<int>(&entity_set_key) );
    entity_set = builder->create( entity_set_key );
    entity_set->assignCommunicator( comm );
    TEST_EQUALITY( entity_set->physicalDimension(), 0 );

    // Broadcast the points with indirect serialization through the geometric
    // entity api.
    Teuchos::broadcast( *comm, 0, points() );

    // Add the points to the entity set.
    entity_set->addEntity( points[0] );
    entity_set->addEntity( points[1] );
    TEST_EQUALITY( entity_set->physicalDimension(), 3 );

    // Check the bounding boxes.
    Teuchos::Tuple<double,6> local_bounds;
    entity_set->localBoundingBox( local_bounds );
    TEST_EQUALITY( local_bounds[0], x_2 );
    TEST_EQUALITY( local_bounds[1], y_2 );
    TEST_EQUALITY( local_bounds[2], z_2 );
    TEST_EQUALITY( local_bounds[3], x_1 );
    TEST_EQUALITY( local_bounds[4], y_1 );
    TEST_EQUALITY( local_bounds[5], z_1 );

    Teuchos::Tuple<double,6> global_bounds;
    entity_set->globalBoundingBox( global_bounds );
    TEST_FLOATING_EQUALITY( global_bounds[0], x_2, 1.0e-12 );
    TEST_FLOATING_EQUALITY( global_bounds[1], y_2, 1.0e-12 );
    TEST_FLOATING_EQUALITY( global_bounds[2], z_2, 1.0e-12 );
    TEST_FLOATING_EQUALITY( global_bounds[3], x_1, 1.0e-12 );
    TEST_FLOATING_EQUALITY( global_bounds[4], y_1, 1.0e-12 );
    TEST_FLOATING_EQUALITY( global_bounds[5], z_1, 1.0e-12 );
}

//---------------------------------------------------------------------------//
// end tstPoint.cpp
//---------------------------------------------------------------------------//

