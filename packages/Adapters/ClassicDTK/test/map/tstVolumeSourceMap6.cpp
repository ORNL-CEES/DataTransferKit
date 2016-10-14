//---------------------------------------------------------------------------//
/*!
 * \file tstVolumeSourceMap6.cpp
 * \author Stuart R. Slattery
 * \brief Volume source map unit test 6 for repeated geometry transfer.
 */
//---------------------------------------------------------------------------//

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <vector>

#include <DTK_Box.hpp>
#include <DTK_Cylinder.hpp>
#include <DTK_FieldContainer.hpp>
#include <DTK_FieldEvaluator.hpp>
#include <DTK_FieldManager.hpp>
#include <DTK_FieldTools.hpp>
#include <DTK_FieldTraits.hpp>
#include <DTK_GeometryManager.hpp>
#include <DTK_GeometryTraits.hpp>
#include <DTK_VolumeSourceMap.hpp>

#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_Ptr.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TypeTraits.hpp>
#include <Teuchos_UnitTestHarness.hpp>

//---------------------------------------------------------------------------//
// MPI Setup
//---------------------------------------------------------------------------//

template <class Ordinal>
Teuchos::RCP<const Teuchos::Comm<Ordinal>> getDefaultComm()
{
#ifdef HAVE_MPI
    return Teuchos::DefaultComm<Ordinal>::getComm();
#else
    return Teuchos::rcp( new Teuchos::SerialComm<Ordinal>() );
#endif
}

//---------------------------------------------------------------------------//
// FieldEvaluator Implementation.
class MyEvaluator
    : public DataTransferKit::FieldEvaluator<
          unsigned long int, DataTransferKit::FieldContainer<double>>
{
  public:
    MyEvaluator( const Teuchos::ArrayRCP<unsigned long int> &geom_gids,
                 const Teuchos::RCP<const Teuchos::Comm<int>> &comm )
        : d_geom_gids( geom_gids )
        , d_comm( comm )
    { /* ... */
    }

    ~MyEvaluator() { /* ... */}

    DataTransferKit::FieldContainer<double>
    evaluate( const Teuchos::ArrayRCP<unsigned long int> &gids,
              const Teuchos::ArrayRCP<double> &coords )
    {
        Teuchos::ArrayRCP<double> evaluated_data( gids.size() );
        for ( int n = 0; n < gids.size(); ++n )
        {
            if ( std::find( d_geom_gids.begin(), d_geom_gids.end(), gids[n] ) !=
                 d_geom_gids.end() )
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
    Teuchos::ArrayRCP<unsigned long int> d_geom_gids;
    Teuchos::RCP<const Teuchos::Comm<int>> d_comm;
};

//---------------------------------------------------------------------------//
// Unit tests. This is a many-to-one transfer.
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( VolumeSourceMap, cylinder_test )
{
    using namespace DataTransferKit;
    typedef FieldContainer<double> FieldType;

    // Setup communication.
    Teuchos::RCP<const Teuchos::Comm<int>> comm = getDefaultComm<int>();

    // Setup source geometry. Same on every proc.
    int geom_dim = 3;
    int num_geom = 4;
    double length = 2.5;
    double radius = 0.75;
    double center_z = 0.25;
    Teuchos::ArrayRCP<Cylinder> geometry( num_geom );
    geometry[0] = Cylinder( length, radius, -1.5, -1.5, center_z );
    geometry[1] = Cylinder( length, radius, 1.5, -1.5, center_z );
    geometry[2] = Cylinder( length, radius, 1.5, 1.5, center_z );
    geometry[3] = Cylinder( length, radius, -1.5, 1.5, center_z );

    Teuchos::ArrayRCP<unsigned long int> geom_gids( num_geom );
    for ( int i = 0; i < num_geom; ++i )
    {
        geom_gids[i] = i;
    }

    Teuchos::RCP<GeometryManager<Cylinder, unsigned long int>>
        source_geometry_manager =
            Teuchos::rcp( new GeometryManager<Cylinder, unsigned long int>(
                geometry, geom_gids, comm, geom_dim ) );

    Teuchos::RCP<FieldEvaluator<unsigned long int, FieldType>>
        source_evaluator = Teuchos::rcp( new MyEvaluator( geom_gids, comm ) );

    // Setup target coords on proc 0 only. Use the geometry centroids plus bogus
    // point.
    Teuchos::ArrayRCP<double> target_coords( 0 );
    if ( comm->getRank() == 0 )
    {
        target_coords =
            Teuchos::ArrayRCP<double>( ( num_geom + 1 ) * geom_dim );
        for ( int i = 0; i < num_geom; ++i )
        {
            target_coords[i] = geometry[i].centroid()[0];
            target_coords[i + num_geom + 1] = geometry[i].centroid()[1];
            target_coords[i + 2 * ( num_geom + 1 )] = geometry[i].centroid()[2];
        }
        target_coords[num_geom] = std::numeric_limits<int>::max();
        target_coords[2 * num_geom + 1] = std::numeric_limits<int>::max();
        target_coords[2 * ( num_geom + 1 ) + num_geom] =
            std::numeric_limits<int>::max();
    }
    Teuchos::RCP<FieldType> coord_field =
        Teuchos::rcp( new FieldType( target_coords, geom_dim ) );

    Teuchos::RCP<FieldManager<FieldType>> target_coord_manager =
        Teuchos::rcp( new FieldManager<FieldType>( coord_field, comm ) );

    // Setup target field.
    int target_field_dim = 1;
    Teuchos::ArrayRCP<double> target_data( 0 );
    if ( comm->getRank() == 0 )
    {
        target_data = Teuchos::ArrayRCP<double>( num_geom + 1 );
    }
    Teuchos::RCP<FieldType> target_field =
        Teuchos::rcp( new FieldType( target_data, target_field_dim ) );

    Teuchos::RCP<FieldManager<FieldType>> target_space_manager =
        Teuchos::rcp( new FieldManager<FieldType>( target_field, comm ) );

    // Setup and apply the volume source mapping.
    VolumeSourceMap<Cylinder, unsigned long int, FieldType> volume_source_map(
        comm, geom_dim, true, 1.0e-6 );
    volume_source_map.setup( source_geometry_manager, target_coord_manager );
    volume_source_map.apply( source_evaluator, target_space_manager );

    // Check the evaluation.
    if ( comm->getRank() == 0 )
    {
        for ( int i = 0; i < num_geom; ++i )
        {
            TEST_ASSERT( target_data[i] == 1.0 + i );
        }
        TEST_ASSERT( target_data[num_geom] == 0.0 );

        // Make sure all points were found except the bogus point.
        TEST_ASSERT( volume_source_map.getMissedTargetPoints().size() == 1 );
        TEST_EQUALITY(
            Teuchos::as<int>( volume_source_map.getMissedTargetPoints()[0] ),
            num_geom );
    }
}

//---------------------------------------------------------------------------//
// end tstVolumeSourceMap6.cpp
//---------------------------------------------------------------------------//
