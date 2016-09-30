//---------------------------------------------------------------------------//
/*!
 * \file tstVolumeSourceMap2.cpp
 * \author Stuart R. Slattery
 * \brief Volume source map unit test 1 for signed ordinals.
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
// Field implementation.
//---------------------------------------------------------------------------//
class MyField
{
  public:

    typedef double value_type;
    typedef Teuchos::Array<double>::size_type size_type;
    typedef Teuchos::Array<double>::iterator iterator;
    typedef Teuchos::Array<double>::const_iterator const_iterator;

    MyField( size_type size, int dim )
        : d_dim( dim )
        , d_data( dim*size, 0.0 )
    { /* ... */ }

    ~MyField()
    { /* ... */ }

    int dim() const
    { return d_dim; }

    size_type size() const
    { return d_data.size(); }

    bool empty() const
    { return d_data.empty(); }

    iterator begin()
    { return d_data.begin(); }

    const_iterator begin() const
    { return d_data.begin(); }

    iterator end()
    { return d_data.end(); }

    const_iterator end() const
    { return d_data.end(); }

    Teuchos::Array<double>& getData()
    { return d_data; }

    const Teuchos::Array<double>& getData() const
    { return d_data; }

  private:
    int d_dim;
    Teuchos::Array<double> d_data;
};

//---------------------------------------------------------------------------//
// DTK implementations.
//---------------------------------------------------------------------------//
namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Field Traits specification for MyField
template<>
class FieldTraits<MyField>
{
  public:

    typedef MyField                    field_type;
    typedef double                     value_type;
    typedef MyField::size_type         size_type;
    typedef MyField::iterator          iterator;
    typedef MyField::const_iterator    const_iterator;

    static inline size_type dim( const MyField& field )
    { return field.dim(); }

    static inline size_type size( const MyField& field )
    { return field.size(); }

    static inline bool empty( const MyField& field )
    { return field.empty(); }

    static inline iterator begin( MyField& field )
    { return field.begin(); }

    static inline const_iterator begin( const MyField& field )
    { return field.begin(); }

    static inline iterator end( MyField& field )
    { return field.end(); }

    static inline const_iterator end( const MyField& field )
    { return field.end(); }
};

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// FieldEvaluator Implementation.
class MyEvaluator :
    public DataTransferKit::FieldEvaluator<unsigned long int,MyField>
{
  public:

    MyEvaluator( const Teuchos::ArrayRCP<unsigned long int>& geom_gids,
                 const Teuchos::RCP< const Teuchos::Comm<int> >& comm )
        : d_geom_gids( geom_gids )
        , d_comm( comm )
    { /* ... */ }

    ~MyEvaluator()
    { /* ... */ }

    MyField evaluate(
        const Teuchos::ArrayRCP<unsigned long int>& gids,
        const Teuchos::ArrayRCP<double>& coords )
    {
        MyField evaluated_data( gids.size(), 1 );
        for ( int n = 0; n < gids.size(); ++n )
        {
            if ( std::find( d_geom_gids.begin(),
                            d_geom_gids.end(),
                            gids[n] ) != d_geom_gids.end() )
            {
                *(evaluated_data.begin() + n ) = 1.0;
            }
            else
            {
                *(evaluated_data.begin() + n ) = 0.0;
            }
        }
        return evaluated_data;
    }

  private:

    Teuchos::ArrayRCP<unsigned long int> d_geom_gids;
    Teuchos::RCP< const Teuchos::Comm<int> > d_comm;
};

//---------------------------------------------------------------------------//
// Geometry create functions. These geometries will span the entire domain,
// requiring them to be broadcast throughout the rendezvous.
//---------------------------------------------------------------------------//
void buildCylinderGeometry(
    int my_size, int edge_size,
    Teuchos::ArrayRCP<DataTransferKit::Cylinder>& cylinders,
    Teuchos::ArrayRCP<unsigned long int>& gids )
{
    Teuchos::ArrayRCP<DataTransferKit::Cylinder> new_cylinders(1);
    Teuchos::ArrayRCP<unsigned long int> new_gids(1,0);
    double length = (double) my_size;
    double radius = (double) (edge_size-1) / 2.0;
    double x_center = (double) (edge_size-1) / 2.0;
    double y_center = (double) (edge_size-1) / 2.0;
    double z_center = (double) my_size / 2.0;
    new_cylinders[0] = DataTransferKit::Cylinder( length, radius,
                                                  x_center, y_center, z_center );
    cylinders = new_cylinders;
    gids = new_gids;
}

//---------------------------------------------------------------------------//
void buildBoxGeometry( int my_size, int edge_size,
                       Teuchos::ArrayRCP<DataTransferKit::Box>& boxes,
                       Teuchos::ArrayRCP<unsigned long int>& gids )
{
    Teuchos::ArrayRCP<DataTransferKit::Box> new_boxes(1);
    Teuchos::ArrayRCP<unsigned long int> new_gids(1,0);
    new_boxes[0] = DataTransferKit::Box( 0.0, 0.0, 0.0, edge_size-1,
                                         edge_size-1, my_size );
    boxes = new_boxes;
    gids = new_gids;
}

//---------------------------------------------------------------------------//
// Coordinate field create function.
//---------------------------------------------------------------------------//
void buildCoordinateField( int my_rank, int my_size,
                           int num_points, int edge_size,
                           Teuchos::RCP<MyField>& coordinate_field )
{
    Teuchos::ArrayRCP<double> coords =
        DataTransferKit::FieldTools<MyField>::nonConstView( *coordinate_field );
    std::srand( my_rank*num_points*2 );
    for ( int i = 0; i < num_points; ++i )
    {
        coords[i] =
            my_size * (edge_size-1) * (double) std::rand() / RAND_MAX;
        coords[num_points + i] =
            (edge_size-1) * (double) std::rand() / RAND_MAX;
        coords[2*num_points + i] =
            (double) std::rand() / RAND_MAX;
    }
}

//---------------------------------------------------------------------------//
// Unit tests.
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( VolumeSourceMap, cylinder_test )
{
    using namespace DataTransferKit;

    // Setup communication.
    Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm<int>();
    int my_rank = comm->getRank();
    int my_size = comm->getSize();
    Teuchos::Array<int> geom_ranks(1,0);
    Teuchos::RCP<const Teuchos::Comm<int> > geom_comm =
        comm->createSubcommunicator( geom_ranks() );

    // Setup source geometry.
    double edge_size = 1.0;
    int geom_dim = 3;
    Teuchos::ArrayRCP<Cylinder> geometry(0);
    Teuchos::ArrayRCP<unsigned long int> geom_gids(0);
    Teuchos::RCP<FieldEvaluator<unsigned long int,MyField> > source_evaluator;
    Teuchos::RCP<GeometryManager<Cylinder,unsigned long int> > source_geometry_manager;
    if ( my_rank == 0 )
    {
        buildCylinderGeometry( my_size, edge_size, geometry, geom_gids );
            source_evaluator = Teuchos::rcp( new MyEvaluator( geom_gids, geom_comm ) );
        source_geometry_manager =
            Teuchos::rcp( new GeometryManager<Cylinder,unsigned long int>(
                              geometry, geom_gids, geom_comm, geom_dim ) );
    }
    comm->barrier();

    // Setup target.
    int target_dim = 1;
    int num_target_points = 100;
    Teuchos::RCP<MyField> target_coords =
        Teuchos::rcp( new MyField( num_target_points, geom_dim ) );
    Teuchos::RCP<MyField> target_field =
        Teuchos::rcp( new MyField( num_target_points, target_dim ) );

    buildCoordinateField( my_rank, my_size, num_target_points, edge_size,
                          target_coords );

    Teuchos::RCP<FieldManager<MyField> > target_coord_manager = Teuchos::rcp(
        new FieldManager<MyField>( target_coords, comm ) );

    Teuchos::RCP<FieldManager<MyField> > target_space_manager = Teuchos::rcp(
        new FieldManager<MyField>( target_field, comm ) );

    // Setup and apply the volume source mapping.
    VolumeSourceMap<Cylinder,unsigned long int,MyField> volume_source_map(
        comm, geom_dim, true, 1.0e-6 );
    volume_source_map.setup( source_geometry_manager, target_coord_manager );
    volume_source_map.apply( source_evaluator, target_space_manager );

    // Check the evaluation.
    Cylinder global_cylinder;
    if ( my_rank == 0 )
    {
        global_cylinder = geometry[0];
    }
    comm->barrier();
    Teuchos::broadcast( *comm, 0, Teuchos::Ptr<Cylinder>(&global_cylinder) );

    Teuchos::ArrayRCP<const double> coords =
        FieldTools<MyField>::view( *target_coords );

    Teuchos::ArrayRCP<const double> target_data =
        FieldTools<MyField>::view( *target_field );

    Teuchos::Array<double> vertex(3);
    double tol = 1.0e-6;
    int num_in_cylinder = 0;
    for ( int i = 0; i < num_target_points; ++i )
    {
        vertex[0] = coords[i];
        vertex[1] = coords[i + num_target_points];
        vertex[2] = coords[i + 2*num_target_points];

        if ( global_cylinder.pointInCylinder( vertex, tol ) )
        {
            ++num_in_cylinder;

            TEST_ASSERT( target_data[i] == 1.0 );
        }
        else
        {
            TEST_ASSERT( target_data[i] == 0.0 );
        }
    }
    comm->barrier();

    int global_num_in_cylinder = 0;
    Teuchos::reduceAll( *comm, Teuchos::REDUCE_SUM, 1,
                        &num_in_cylinder, &global_num_in_cylinder );

    int num_missed = volume_source_map.getMissedTargetPoints().size();
    int global_num_missed = 0;
    Teuchos::reduceAll( *comm, Teuchos::REDUCE_SUM, 1,
                        &num_missed, &global_num_missed );

    TEST_ASSERT( num_target_points*my_size ==
                 global_num_missed + global_num_in_cylinder );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( VolumeSourceMap, box_test )
{
    using namespace DataTransferKit;

    // Setup communication.
    Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm<int>();
    int my_rank = comm->getRank();
    int my_size = comm->getSize();
    Teuchos::Array<int> geom_ranks(1,0);
    Teuchos::RCP<const Teuchos::Comm<int> > geom_comm =
        comm->createSubcommunicator( geom_ranks() );

    // Setup source geometry.
    double edge_size = 1.0;
    int geom_dim = 3;
    Teuchos::ArrayRCP<Box> geometry(0);
    Teuchos::ArrayRCP<unsigned long int> geom_gids(0);
    Teuchos::RCP<FieldEvaluator<unsigned long int,MyField> > source_evaluator;
    Teuchos::RCP<GeometryManager<Box,unsigned long int> > source_geometry_manager;
    if ( my_rank == 0 )
    {
        buildBoxGeometry( my_size, edge_size, geometry, geom_gids );
            source_evaluator = Teuchos::rcp( new MyEvaluator( geom_gids, geom_comm ) );
        source_geometry_manager =
            Teuchos::rcp( new GeometryManager<Box,unsigned long int>(
                              geometry, geom_gids, geom_comm, geom_dim ) );
    }
    comm->barrier();

    // Setup target.
    int target_dim = 1;
    int num_target_points = 100;
    Teuchos::RCP<MyField> target_coords =
        Teuchos::rcp( new MyField( num_target_points, geom_dim ) );
    Teuchos::RCP<MyField> target_field =
        Teuchos::rcp( new MyField( num_target_points, target_dim ) );

    buildCoordinateField( my_rank, my_size, num_target_points, edge_size,
                          target_coords );

    Teuchos::RCP<FieldManager<MyField> > target_coord_manager = Teuchos::rcp(
        new FieldManager<MyField>( target_coords, comm ) );

    Teuchos::RCP<FieldManager<MyField> > target_space_manager = Teuchos::rcp(
        new FieldManager<MyField>( target_field, comm ) );

    // Setup and apply the volume source mapping.
    VolumeSourceMap<Box,unsigned long int,MyField> volume_source_map(
        comm, geom_dim, true, 1.0e-6 );
    volume_source_map.setup( source_geometry_manager, target_coord_manager );
    volume_source_map.apply( source_evaluator, target_space_manager );

    // Check the evaluation.
    Box global_box;
    if ( my_rank == 0 )
    {
        global_box = geometry[0];
    }
    comm->barrier();
    Teuchos::broadcast( *comm, 0, Teuchos::Ptr<Box>(&global_box) );

    Teuchos::ArrayRCP<const double> coords =
        FieldTools<MyField>::view( *target_coords );

    Teuchos::ArrayRCP<const double> target_data =
        FieldTools<MyField>::view( *target_field );

    Teuchos::Array<double> vertex(3);
    double tol = 1.0e-6;
    int num_in_box = 0;
    for ( int i = 0; i < num_target_points; ++i )
    {
        vertex[0] = coords[i];
        vertex[1] = coords[i + num_target_points];
        vertex[2] = coords[i + 2*num_target_points];

        if ( global_box.pointInBox( vertex, tol ) )
        {
            ++num_in_box;

            TEST_ASSERT( target_data[i] == 1.0 );
        }
        else
        {
            TEST_ASSERT( target_data[i] == 0.0 );
        }
    }
    comm->barrier();

    int global_num_in_box = 0;
    Teuchos::reduceAll( *comm, Teuchos::REDUCE_SUM, 1,
                        &num_in_box, &global_num_in_box );

    int num_missed = volume_source_map.getMissedTargetPoints().size();
    int global_num_missed = 0;
    Teuchos::reduceAll( *comm, Teuchos::REDUCE_SUM, 1,
                        &num_missed, &global_num_missed );

    TEST_ASSERT( num_target_points*my_size ==
                 global_num_missed + global_num_in_box );
}

//---------------------------------------------------------------------------//
// end tstVolumeSourceMap2.cpp
//---------------------------------------------------------------------------//
