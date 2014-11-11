//---------------------------------------------------------------------------//
/*! 
 * \file tstMapOperator.cpp
 * \author Stuart R. Slattery
 * \brief MapOperator unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <DTK_MapOperator.hpp>
#include <DTK_FunctionSpace.hpp>
#include <DTK_BasicEntitySet.hpp>
#include <DTK_BasicGeometryLocalMap.hpp>
#include <DTK_EntityCenteredShapeFunction.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_TypeTraits.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Tpetra_Map.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_MultiVector.hpp>

#include <Thyra_VectorSpaceBase.hpp>
#include <Thyra_MultiVectorBase.hpp>
#include <Thyra_LinearOpBase.hpp>
#include <Thyra_MultiVectorStdOps.hpp>
#include <Thyra_TpetraThyraWrappers.hpp>
#include <Thyra_TpetraLinearOp.hpp>

//---------------------------------------------------------------------------//
// Create a diagonally scaled linear operator.
//---------------------------------------------------------------------------//
Teuchos::RCP<Thyra::LinearOpBase<double> > 
createScaledIdentity( 
    const Teuchos::RCP<const Teuchos::Comm<int> >& comm, 
    int local_size,
    double val )
{
    Teuchos::RCP<const Tpetra::Map<int,std::size_t> > map = 
	Tpetra::createUniformContigMap<int,std::size_t>( local_size*comm->getSize(),
						 comm );
    Teuchos::RCP<Tpetra::CrsMatrix<double,int,std::size_t> > mat = Teuchos::rcp( 
	new Tpetra::CrsMatrix<double,int,std::size_t>(map, 1, Tpetra::StaticProfile) );
    Teuchos::Array<double> values( 1, val );
    Teuchos::Array<std::size_t> indices( 1 );
    int global_row = 0;
    for ( int i = 0; i < local_size; ++i )
    {
	global_row = local_size*comm->getRank() + i;
	indices[0] = global_row;
	mat->insertGlobalValues( global_row, indices(), values() );
    }
    mat->fillComplete();

    Teuchos::RCP<const Thyra::VectorSpaceBase<double> > space =
	Thyra::createVectorSpace<double>( map );

    Teuchos::RCP<Thyra::TpetraLinearOp<double,int,std::size_t> > tpetra_op =
	Teuchos::rcp( new Thyra::TpetraLinearOp<double,int,std::size_t>() );
    tpetra_op->initialize( space, space, mat );

    return tpetra_op;
}

//---------------------------------------------------------------------------//
// Create a vector filled with a scalar.
//---------------------------------------------------------------------------//
Teuchos::RCP<Thyra::MultiVectorBase<double> > 
createDOFVector( 
    const Teuchos::RCP<const Teuchos::Comm<int> >& comm, 
    int local_size,
    double val,
    int num_vec )
{
    Teuchos::RCP<const Tpetra::Map<int,std::size_t> > map = 
	Tpetra::createUniformContigMap<int,std::size_t>( local_size*comm->getSize(),
						 comm );
    Teuchos::RCP<Tpetra::MultiVector<double,int,std::size_t> > vec =
	Tpetra::createMultiVector<double>( map, num_vec );
    vec->putScalar( val );

    return Thyra::createMultiVector( vec );
}

//---------------------------------------------------------------------------//
// Map operator implementation.
//---------------------------------------------------------------------------//
class TestOperator : public DataTransferKit::MapOperator<double>
{
  public:

    TestOperator( const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
		  int local_size,
		  double scalar_a,
		  double scalar_b,
		  double scalar_c,
		  int num_vec )
	: d_comm( comm )
	, d_local_size( local_size )
	, d_a( scalar_a )
	, d_b( scalar_b )
	, d_c( scalar_c )
	, d_num_vec( num_vec )
    { /* ... */ }

    ~TestOperator() { /* ... */ }

    void setup( const Teuchos::RCP<DataTransferKit::FunctionSpace>& domain_space,
		const Teuchos::RCP<DataTransferKit::FunctionSpace>& range_space,
		const Teuchos::RCP<Teuchos::ParameterList>& parameters )
    {
	this->b_mass_matrix_inv = createScaledIdentity( d_comm, d_local_size, d_a );
	this->b_coupling_matrix = createScaledIdentity( d_comm, d_local_size, d_b );
	this->b_forcing_vector = createDOFVector( d_comm, d_local_size, d_c, d_num_vec );
    }

    Teuchos::RCP<const Thyra::LinearOpBase<double> > clone() const
    {
	return Teuchos::rcp( 
	    new TestOperator(d_comm,d_local_size,d_a,d_b,d_c,d_num_vec) );
    }

  private:
    Teuchos::RCP<const Teuchos::Comm<int> > d_comm;
    int d_local_size;
    double d_a;
    double d_b;
    double d_c;
    int d_num_vec;
};

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( MapOperator, scaled_identity_test )
{
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
	Teuchos::DefaultComm<int>::getComm();

    int local_size = 10;
    double scalar = 2.0;
    int num_vec = 3;

    Teuchos::RCP<Thyra::LinearOpBase<double> > identity_op =
	createScaledIdentity( comm, local_size, scalar );

    // Test the apply on some vectors.
    Teuchos::RCP<const Tpetra::Map<int,std::size_t> > map = 
	Tpetra::createUniformContigMap<int,std::size_t>( local_size*comm->getSize(),
						 comm );
    Teuchos::RCP<Tpetra::MultiVector<double,int,std::size_t> > X =
	Tpetra::createMultiVector<double,int,std::size_t>( map, num_vec );
    X->randomize();
    Teuchos::RCP<Tpetra::MultiVector<double,int,std::size_t> > Y =
	Tpetra::createMultiVector<double,int,std::size_t>( map, num_vec );

    Teuchos::RCP<Thyra::MultiVectorBase<double> > thyra_X =
	Thyra::createMultiVector( X );	
    Teuchos::RCP<Thyra::MultiVectorBase<double> > thyra_Y =
	Thyra::createMultiVector( Y );

    identity_op->apply( Thyra::NOTRANS, *thyra_X, thyra_Y.ptr(), 1.0, 0.0 );

    for ( int n = 0; n < num_vec; ++n )
    {
	Teuchos::ArrayRCP<const double> x_view = X->getVector(n)->getData();
	Teuchos::ArrayRCP<const double> y_view = Y->getVector(n)->getData();
	for ( int i = 0; i < local_size; ++i )
	{
	    TEST_FLOATING_EQUALITY( y_view[i], scalar*x_view[i], 1.0e-14 );
	}
    }
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( MapOperator, apply_test )
{
    using namespace DataTransferKit;
    
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
	Teuchos::DefaultComm<int>::getComm();

    // Parameters.
    int local_size = 10;
    double a = 3.2;
    double b = 2.3;
    double c = -1.2;
    int num_vec = 3;

    // Make a function space.
    Teuchos::RCP<EntitySet> entity_set = Teuchos::rcp(
	new BasicEntitySet(comm, 1) );
    Teuchos::RCP<EntityLocalMap> local_map = 
	Teuchos::rcp( new BasicGeometryLocalMap() );
    Teuchos::RCP<EntityShapeFunction> shape_function =
	Teuchos::rcp( new EntityCenteredShapeFunction() );
    Teuchos::RCP<const Tpetra::Map<int,std::size_t> > dof_map = 
	Tpetra::createUniformContigMap<int,std::size_t>( 
	    local_size*comm->getSize(), comm );
    Teuchos::RCP<FunctionSpace> function_space = Teuchos::rcp( 
	new FunctionSpace(entity_set, local_map, shape_function, dof_map) );

    // Make a map.
    Teuchos::RCP<MapOperator<double> > map_op = Teuchos::rcp(
	new TestOperator(comm,local_size,a,b,c,num_vec) );
    map_op->setup( function_space, function_space, Teuchos::parameterList() );

    // Test the apply on some vectors.
    Teuchos::RCP<Tpetra::MultiVector<double,int,std::size_t> > X =
	Tpetra::createMultiVector<double,int,std::size_t>( dof_map, num_vec );
    X->randomize();
    Teuchos::RCP<Tpetra::MultiVector<double,int,std::size_t> > Y =
	Tpetra::createMultiVector<double,int,std::size_t>( dof_map, num_vec );
    Y->putScalar( 1.0 );

    Teuchos::RCP<Thyra::MultiVectorBase<double> > thyra_X =
	Thyra::createMultiVector( X );	
    Teuchos::RCP<Thyra::MultiVectorBase<double> > thyra_Y =
	Thyra::createMultiVector( Y );

    double alpha = 1.3;
    double beta = -2.1;
    map_op->apply( Thyra::NOTRANS, *thyra_X, thyra_Y.ptr(), alpha, beta );

    for ( int n = 0; n < num_vec; ++n )
    {
	Teuchos::ArrayRCP<const double> x_view = X->getVector(n)->getData();
	Teuchos::ArrayRCP<const double> y_view = Y->getVector(n)->getData();
	for ( int i = 0; i < local_size; ++i )
	{
	    TEST_FLOATING_EQUALITY( y_view[i], alpha*a*(c + x_view[i]*b) + beta,
				    1.0e-14 );
	}
    }

    // Clone the operator and apply again.
    Teuchos::RCP<const Thyra::LinearOpBase<double> > map_op_clone =
	map_op->clone();
    Teuchos::rcp_const_cast<MapOperator<double> >(
	Teuchos::rcp_dynamic_cast<const MapOperator<double> >(map_op_clone) )->setup( 
	function_space, function_space, Teuchos::parameterList() );

    Teuchos::RCP<Tpetra::MultiVector<double,int,std::size_t> > Z =
	Tpetra::createMultiVector<double,int,std::size_t>( dof_map, num_vec );
    Z->putScalar( 1.0 );

    Teuchos::RCP<Thyra::MultiVectorBase<double> > thyra_Z =
	Thyra::createMultiVector( Z );

    map_op_clone->apply( Thyra::NOTRANS, *thyra_X, thyra_Z.ptr(), alpha, beta );

    for ( int n = 0; n < num_vec; ++n )
    {
	Teuchos::ArrayRCP<const double> x_view = X->getVector(n)->getData();
	Teuchos::ArrayRCP<const double> z_view = Z->getVector(n)->getData();
	for ( int i = 0; i < local_size; ++i )
	{
	    TEST_FLOATING_EQUALITY( z_view[i], alpha*a*(c + x_view[i]*b) + beta,
				    1.0e-14 );
	}
    }

}

//---------------------------------------------------------------------------//
// end tstMapOperator.cpp
//---------------------------------------------------------------------------//
