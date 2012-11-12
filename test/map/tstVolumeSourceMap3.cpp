//---------------------------------------------------------------------------//
/*!
 * \file tstVolumeSourceMap3.cpp
 * \author Stuart R. Slattery
 * \brief Volume source map unit test 3 for repeated geometry transfer.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <cassert>
#include <ctime>
#include <cstdlib>

#include <DTK_VolumeSourceMap.hpp>
#include <DTK_FieldTraits.hpp>
#include <DTK_FieldEvaluator.hpp>
#include <DTK_FieldManager.hpp>
#include <DTK_FieldTools.hpp>
#include <DTK_FieldContainer.hpp>
#include <DTK_GeometryTraits.hpp>
#include <DTK_GeometryManager.hpp>
#include <DTK_Cylinder.hpp>
#include <DTK_Box.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Ptr.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayView.hpp>
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
// FieldEvaluator Implementation.
class MyEvaluator : 
    public DataTransferKit::FieldEvaluator<int,DataTransferKit::FieldContainer<double> >
{
  public:

    MyEvaluator( const Teuchos::ArrayRCP<int>& geom_gids, 
		 const Teuchos::RCP< const Teuchos::Comm<int> >& comm )
	: d_geom_gids( geom_gids )
	, d_comm( comm )
    { /* ... */ }

    ~MyEvaluator()
    { /* ... */ }

    DataTransferKit::FieldContainer<double> evaluate( 
	const Teuchos::ArrayRCP<int>& gids,
	const Teuchos::ArrayRCP<double>& coords )
    {
	Teuchos::ArrayRCP<double> evaluated_data( gids.size() );
	for ( int n = 0; n < gids.size(); ++n )
	{
	    if ( std::find( d_geom_gids.begin(),
			    d_geom_gids.end(),
			    gids[n] ) != d_geom_gids.end() )
	    {
		evaluated_data[n] = 1.0 + gids[n];
	    }
	    else
	    {
		evaluated_data[n] = 0.0;
	    }
	}
	return DataTransferKit::FieldContainer<double>( evaluated_data, 1 );
    }

  private:

    Teuchos::ArrayRCP<int> d_geom_gids;
    Teuchos::RCP< const Teuchos::Comm<int> > d_comm;
};

//---------------------------------------------------------------------------//
// Unit tests.
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( VolumeSourceMap, cylinder_test )
{
    using namespace DataTransferKit;
    typedef FieldContainer<double> FieldType;

    // Setup communication.
    Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm<int>();

    // Setup source geometry. Same on every proc.
    int geom_dim = 3;
    int num_geom = 4;
    double length = 2.5;
    double radius = 0.75;
    double center_z = 0.25;
    Teuchos::ArrayRCP<Cylinder> geometry(num_geom);
    geometry[0] = Cylinder( length, radius, -1.5, -1.5, center_z );
    geometry[1] = Cylinder( length, radius,  1.5, -1.5, center_z );
    geometry[2] = Cylinder( length, radius,  1.5,  1.5, center_z );
    geometry[3] = Cylinder( length, radius, -1.5,  1.5, center_z );

    Teuchos::ArrayRCP<int> geom_gids(num_geom);
    for ( int i = 0; i < num_geom; ++i )
    {
	geom_gids[i] = i;
    }

    Teuchos::RCP<GeometryManager<Cylinder,int> > source_geometry_manager =
	Teuchos::rcp( new GeometryManager<Cylinder,int>( 
			      geometry, geom_gids, comm, geom_dim ) );

    Teuchos::RCP<FieldEvaluator<int,FieldType> > source_evaluator = 
	Teuchos::rcp( new MyEvaluator( geom_gids, comm ) );

    // Setup target coords. Use the geometry centroids.
    Teuchos::ArrayRCP<double> target_coords( num_geom*geom_dim );
    for ( int i = 0; i < num_geom; ++i )
    {
	target_coords[i] = geometry[i].centroid()[0];
	target_coords[i + num_geom] = geometry[i].centroid()[1];
	target_coords[i + 2*num_geom] = geometry[i].centroid()[2];
    }
    Teuchos::RCP<FieldType > coord_field =
	Teuchos::rcp( new FieldType( target_coords, geom_dim ) );

    Teuchos::RCP<FieldManager<FieldType> > target_coord_manager = 
	Teuchos::rcp( new FieldManager<FieldType>( coord_field, comm ) );

    // Setup target field.
    int target_field_dim = 1;
    Teuchos::ArrayRCP<double> target_data( num_geom );
    Teuchos::RCP<FieldType> target_field =
	Teuchos::rcp( new FieldType( target_data, target_field_dim ) );

    Teuchos::RCP<FieldManager<FieldType> > target_space_manager = 
	Teuchos::rcp( new FieldManager<FieldType>( target_field, comm ) );

    // Setup and apply the volume source mapping.
    VolumeSourceMap<Cylinder,int,FieldType> volume_source_map( 
	comm, geom_dim, true, 1.0e-6 );
    volume_source_map.setup( source_geometry_manager, target_coord_manager );
    volume_source_map.apply( source_evaluator, target_space_manager );

    // Check the evaluation.
    std::cout << comm->getRank() << ": " << target_data() << std::endl;
    for ( int i = 0; i < num_geom; ++i )
    {
	TEST_ASSERT( target_data[i] == 1.0 + i );
    }

    // Make sure all points were found.
    TEST_ASSERT( volume_source_map.getMissedTargetPoints().size() == 0 );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( VolumeSourceMap, box_test )
{
}

//---------------------------------------------------------------------------//
// end tstVolumeSourceMap3.cpp
//---------------------------------------------------------------------------//
