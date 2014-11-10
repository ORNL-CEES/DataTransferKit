//---------------------------------------------------------------------------//
/*! 
 * \file tstConsistentInterpolationOperatorFEM.cpp
 * \author Stuart R. Slattery
 * \brief ConsistentInterpolationOperator unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <DTK_ConsistentInterpolationOperator.hpp>
#include <DTK_EntitySelector.hpp>
#include <DTK_FunctionSpace.hpp>
#include <DTK_BasicEntitySet.hpp>
#include <DTK_BasicGeometryLocalMap.hpp>
#include <DTK_EntityShapeFunction.hpp>
#include <DTK_EntityCenteredShapeFunction.hpp>
#include <DTK_EntityCenteredDOFVector.hpp>
#include <DTK_Box.hpp>
#include <DTK_Point.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_TypeTraits.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>

#include <Thyra_TpetraThyraWrappers.hpp>


//---------------------------------------------------------------------------//
// Shape function implementation.
//---------------------------------------------------------------------------//
class TestShapeFunction : public DataTransferKit::EntityShapeFunction
{
  public:

    TestShapeFunction( const Teuchos::ArrayView<const std::size_t>& entity_ids,
		       const int dofs_per_entity ) 
	: d_entity_ids( entity_ids )
	, d_dofs_per_entity( dofs_per_entity )
    { /* ... */ }
    ~TestShapeFunction() { /* ... */ }

    void allDOFIds( Teuchos::Array<std::size_t>& dof_ids ) const
    {
	dof_ids = d_entity_ids;
    }

    void entityDOFIds( const DataTransferKit::Entity& entity,
		       Teuchos::Array<std::size_t>& dof_ids ) const
    {
	dof_ids.resize( d_dofs_per_entity );
	for ( int i = 0; i < d_dofs_per_entity; ++i )
	{
	    dof_ids[i] = entity.id()*d_dofs_per_entity + i;
	}
    }

    void evaluateValue( 
	const DataTransferKit::Entity& entity,
	const Teuchos::ArrayView<const double>& reference_point,
	Teuchos::Array<double> & values ) const
    {
	values.assign( d_dofs_per_entity, 1.0 / d_dofs_per_entity );
    }

  private:
    Teuchos::Array<std::size_t> d_entity_ids;
    int d_dofs_per_entity;
};

//---------------------------------------------------------------------------//
// DOF vector.
//---------------------------------------------------------------------------//
template<class Scalar>
Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > 
createTestDOFVector( 
    const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
    const Teuchos::ArrayView<const std::size_t>& entity_ids,
    const Teuchos::ArrayRCP<Scalar>& dof_data,
    const std::size_t lda,
    const std::size_t num_vectors )
{
    // Construct a map.
    int num_entity = entity_ids.size();
    int dofs_per_entity = lda / num_entity;
    Teuchos::Array<std::size_t> dof_ids( lda );
    for ( int i = 0; i < num_entity; ++i )
    {
	for ( int n = 0; n < dofs_per_entity; ++n )
	{
	    dof_ids[ i*dofs_per_entity + n ] = entity_ids[i]*dofs_per_entity + n;
	}
    }
    Teuchos::RCP<const Tpetra::Map<int,std::size_t> > map =
	Tpetra::createNonContigMap<int,std::size_t>( dof_ids, comm );

    // Build a tpetra multivector.
    Teuchos::RCP<Tpetra::MultiVector<Scalar,int,std::size_t> > tpetra_mv =
	Tpetra::createMultiVectorFromView<Scalar,int,std::size_t>( 
	    map, dof_data, lda, num_vectors );

    // Return a thyra multivector.
    return Thyra::createMultiVector( tpetra_mv );
}

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( ConsistentInterpolationOperator, one_to_one_test )
{
    using namespace DataTransferKit;

    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
	Teuchos::DefaultComm<int>::getComm();
    int comm_rank = comm->getRank();
    int comm_size = comm->getSize();

    // DOMAIN SETUP
    // Make a domain entity set.
    Teuchos::RCP<EntitySet> domain_set = 
	Teuchos::rcp( new BasicEntitySet(comm,3) );
    int num_boxes = 5;
    int dofs_per_entity = 4;
    Teuchos::Array<std::size_t> box_ids( num_boxes );
    Teuchos::ArrayRCP<double> box_dofs( dofs_per_entity * num_boxes );
    Teuchos::Array<std::size_t> box_dof_ids( dofs_per_entity * num_boxes );
    for ( int i = 0; i < num_boxes; ++i )
    {
	box_ids[i] = num_boxes*(comm_size-comm_rank-1) + i;
	for ( int n = 0; n < dofs_per_entity; ++n )
	{
	    box_dofs[i*dofs_per_entity+n] = 2.0*box_ids[i] + n;
	    box_dof_ids[i*dofs_per_entity+n] = box_ids[i]*dofs_per_entity + n;
	}
	Teuchos::rcp_dynamic_cast<BasicEntitySet>(domain_set)->addEntity(
	    Box(box_ids[i],comm_rank,box_ids[i],
		0.0,0.0,box_ids[i],1.0,1.0,box_ids[i]+1.0) );
    }

    // Construct a local map for the boxes.
    Teuchos::RCP<EntityLocalMap> domain_map = 
	Teuchos::rcp( new BasicGeometryLocalMap() );

    // Construct a shape function for the boxes.
    Teuchos::RCP<EntityShapeFunction> domain_shape =
	Teuchos::rcp( new TestShapeFunction(box_dof_ids,dofs_per_entity) );

    // Construct a function space for the boxes.
    Teuchos::RCP<FunctionSpace> domain_space =
	Teuchos::rcp( new FunctionSpace(domain_set,domain_map,domain_shape) );

    // Construct a selector for the boxes.
    Teuchos::RCP<EntitySelector> domain_selector = 
	Teuchos::rcp( new EntitySelector(ENTITY_TYPE_VOLUME) );

    // Construct a DOF vector for the boxes.
    Teuchos::RCP<Thyra::MultiVectorBase<double> > domain_dofs =
	createTestDOFVector(
	    comm, box_ids, box_dofs, dofs_per_entity*num_boxes, 1 );

    // RANGE SETUP
    // Make a range entity set.
    Teuchos::RCP<EntitySet> range_set =
	Teuchos::rcp( new BasicEntitySet(comm,3) );
    int num_points = 5;
    Teuchos::Array<double> point(3);
    Teuchos::Array<std::size_t> point_ids( num_points );
    Teuchos::ArrayRCP<double> point_dofs( num_points );
    for ( int i = 0; i < num_points; ++i )
    {
	bool on_surface = (i%2==0);
	point_ids[i] = num_points*comm_rank + i;
	point_dofs[i] = 0.0;
	point[0] = 0.5;
	point[1] = 0.5;
	point[2] = point_ids[i] + 0.5;
	Teuchos::rcp_dynamic_cast<BasicEntitySet>(range_set)->addEntity(
	    Point(point_ids[i],comm_rank,point,on_surface) );
    }

    // Construct a local map for the points.
    Teuchos::RCP<EntityLocalMap> range_map = 
	Teuchos::rcp( new BasicGeometryLocalMap() );

    // Construct a shape function for the points.
    Teuchos::RCP<EntityShapeFunction> range_shape =
	Teuchos::rcp( new EntityCenteredShapeFunction(point_ids) );

    // Construct a function space for the points.
    Teuchos::RCP<FunctionSpace> range_space =
	Teuchos::rcp( new FunctionSpace(range_set,range_map,range_shape) );

    // Construct a selector for the points.
    Teuchos::RCP<EntitySelector> range_selector = 
	Teuchos::rcp( new EntitySelector(ENTITY_TYPE_NODE) );

    // Construct a DOF vector for the points.
    Teuchos::RCP<Thyra::MultiVectorBase<double> > range_dofs =
	EntityCenteredDOFVector::createThyraMultiVector(
	    comm, point_ids, point_dofs, num_points, 1 );

    // MAPPING
    // Create a map.
    Teuchos::RCP<MapOperator<double> > map_op = Teuchos::rcp(
	new ConsistentInterpolationOperator<double>(
	    comm,domain_selector,range_selector) );

    // Setup the map.
    Teuchos::RCP<Teuchos::ParameterList> parameters = Teuchos::parameterList();
    map_op->setup( domain_space, range_space, parameters );

    // Apply the map.
    map_op->apply( Thyra::NOTRANS, *domain_dofs, range_dofs.ptr(), 1.0, 0.0 );

    // Check the results of the mapping.
    for ( int i = 0; i < num_points; ++i )
    {
	double test_val = 0.0;
	for ( int n = 0; n < dofs_per_entity; ++n )
	{
	    test_val += 2.0*point_ids[i] + n;
	}
	test_val /= dofs_per_entity;
	TEST_EQUALITY( test_val, point_dofs[i] );
    }
}

//---------------------------------------------------------------------------//
// end tstConsistentInterpolationOperatorFEM.cpp
//---------------------------------------------------------------------------//
