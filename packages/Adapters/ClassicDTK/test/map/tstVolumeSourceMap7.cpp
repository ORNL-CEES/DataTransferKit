//---------------------------------------------------------------------------//
/*!
 * \file tstVolumeSourceMap7.cpp
 * \author Roger P. Pawlowski
 * \brief Volume source map unit test 5. One to many, all targets requested 
 *        on all processes with a single source geometry to test function
 * evaluation capability is maintained.
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
#include <Teuchos_FancyOStream.hpp>

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
    public DataTransferKit::FieldEvaluator<unsigned long int,DataTransferKit::FieldContainer<double> >
{
  public:

    MyEvaluator( const Teuchos::ArrayRCP<unsigned long int>& geom_gids, 
		 const Teuchos::RCP< const Teuchos::Comm<int> >& comm )
	: d_geom_gids( geom_gids )
	, d_comm( comm )
    { /* ... */ }

    ~MyEvaluator()
    { /* ... */ }

    DataTransferKit::FieldContainer<double> evaluate( 
	const Teuchos::ArrayRCP<unsigned long int>& gids,
	const Teuchos::ArrayRCP<double>& coords )
    {
	Teuchos::ArrayRCP<double> evaluated_data( gids.size() );
	for ( int n = 0; n < gids.size(); ++n )
	{
	  const Teuchos::ArrayRCP<int>::size_type num_vals = gids.size();

	  double z = coords[2*num_vals+n];
	  
	  std::cout << "n = " << n << ", z = " << z << std::endl;

	  if ((z >= 0.0) && (z < 1.0))
	    evaluated_data[n] = 1.0;
	  else if ((z >= 1.0) && (z < 2.0))
	    evaluated_data[n] = 2.0;
	  else if ((z >= 2.0) && (z < 3.0))
	    evaluated_data[n] = 3.0;
	  else if ((z >= 3.0) && (z < 4.0))
	    evaluated_data[n] = 4.0;

	  std::cout << "value = " << evaluated_data[n] << std::endl;

	}
	return DataTransferKit::FieldContainer<double>( evaluated_data, 1 );
    }

  private:

    Teuchos::ArrayRCP<unsigned long int> d_geom_gids;
    Teuchos::RCP< const Teuchos::Comm<int> > d_comm;
};

//---------------------------------------------------------------------------//
// Unit tests.
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( VolumeSourceMap, one_to_many_parallel)
{
  using namespace DataTransferKit;
  typedef FieldContainer<double> FieldType;

  Teuchos::FancyOStream pout(Teuchos::rcp(&std::cout,false));
  pout.setShowLinePrefix(true);
  pout.setShowProcRank(true);

  // Setup comms: rank 0 is source app, all ranks are target app 
  Teuchos::RCP<const Teuchos::Comm<int> > global_comm = getDefaultComm<int>();
  Teuchos::Array<int> subcomm_ranks;
  subcomm_ranks.push_back(0);
  Teuchos::RCP<const Teuchos::Comm<int> > source_comm = global_comm->createSubcommunicator(subcomm_ranks);
  Teuchos::RCP<const Teuchos::Comm<int> > target_comm = global_comm;

  // Setup source geometry. Only on proc zero
  const int geom_dim = 3;
  Teuchos::RCP<GeometryManager<Box,unsigned long int> > source_geometry_manager;
  Teuchos::RCP<FieldEvaluator<unsigned long int,FieldType> > source_evaluator;

  if (global_comm->getRank() == 0) {
      const int num_geom = 1;
      Teuchos::ArrayRCP<Box> geometry(num_geom);
      geometry[0] = Box(-1.0, -1.0, 0.0, 1.0, 1.0, 4.0);

      Teuchos::ArrayRCP<unsigned long int> geom_gids(num_geom);
      for ( int i = 0; i < num_geom; ++i )
      {
	  geom_gids[i] = 0;
      }
    
      source_geometry_manager = Teuchos::rcp(new GeometryManager<Box,unsigned long int>(geometry,geom_gids,source_comm,geom_dim));

      source_evaluator = Teuchos::rcp(new MyEvaluator(geom_gids,source_comm));
  }

  // Setup target coords.
  Teuchos::ArrayRCP<double> target_coords(4*geom_dim);
  if (global_comm->getRank() == 0) {
    target_coords[0] = 0.0; // x0
    target_coords[1] = 0.0; // x1
    target_coords[2] = 0.0; // x2
    target_coords[3] = 0.0; // x3
    target_coords[4] = 0.0; // y0
    target_coords[5] = 0.0; // y1
    target_coords[6] = 0.0; // y2
    target_coords[7] = 0.0; // y3
    target_coords[8] = 0.5; // z0
    target_coords[9] = 1.5; // z1
    target_coords[10] = 2.5; // z2
    target_coords[11] = 3.5; // z3
  }
  else if (global_comm->getRank() == 1) {
    target_coords[0] = 0.0; // x0
    target_coords[1] = 0.0; // x1
    target_coords[2] = 0.0; // x2
    target_coords[3] = 0.0; // x3
    target_coords[4] = 0.0; // y0
    target_coords[5] = 0.0; // y1
    target_coords[6] = 0.0; // y2
    target_coords[7] = 0.0; // y3
    target_coords[8] = 0.5; // z0
    target_coords[9] = 1.5; // z1
    target_coords[10] = 2.5; // z2
    target_coords[11] = 3.5; // z3
  }
  else if (global_comm->getRank() == 2) {
    target_coords[0] = 0.0; // x0
    target_coords[1] = 0.0; // x1
    target_coords[2] = 0.0; // x2
    target_coords[3] = 0.0; // x3
    target_coords[4] = 0.0; // y0
    target_coords[5] = 0.0; // y1
    target_coords[6] = 0.0; // y2
    target_coords[7] = 0.0; // y3
    target_coords[8] = 0.5; // z0
    target_coords[9] = 1.5; // z1
    target_coords[10] = 2.5; // z2
    target_coords[11] = 3.5; // z3
  }
  else if (global_comm->getRank() == 3) {
    target_coords[0] = 0.0; // x0
    target_coords[1] = 0.0; // x1
    target_coords[2] = 0.0; // x2
    target_coords[3] = 0.0; // x3
    target_coords[4] = 0.0; // y0
    target_coords[5] = 0.0; // y1
    target_coords[6] = 0.0; // y2
    target_coords[7] = 0.0; // y3
    target_coords[8] = 0.5; // z0
    target_coords[9] = 1.5; // z1
    target_coords[10] = 2.5; // z2
    target_coords[11] = 3.5; // z3
  }

  Teuchos::RCP<FieldType > coord_field =
    Teuchos::rcp( new FieldType( target_coords, geom_dim ) );
  
  Teuchos::RCP<FieldManager<FieldType> > target_coord_manager = 
    Teuchos::rcp( new FieldManager<FieldType>( coord_field, target_comm ) );

  // Setup target field.
  int target_field_dim = 1;
  Teuchos::ArrayRCP<double> target_data(4);

  Teuchos::RCP<FieldType> target_field =
    Teuchos::rcp( new FieldType( target_data, target_field_dim ) );

  Teuchos::RCP<FieldManager<FieldType> > target_space_manager = 
    Teuchos::rcp( new FieldManager<FieldType>( target_field, target_comm ) );

  // Setup and apply the volume source mapping.
  VolumeSourceMap<Box,unsigned long int,FieldType> volume_source_map( 
    global_comm, geom_dim, true, 1.0e-6 );
  volume_source_map.setup( source_geometry_manager, target_coord_manager );
  volume_source_map.apply( source_evaluator, target_space_manager );

  // Check the evaluation.
  pout << global_comm->getRank() << ": " << target_data() << std::endl;

  double tol = 100.0 * Teuchos::ScalarTraits<double>::eps();

  TEUCHOS_TEST_FLOATING_EQUALITY(target_data[0],1.0,tol,pout,success);
  TEUCHOS_TEST_FLOATING_EQUALITY(target_data[1],2.0,tol,pout,success);
  TEUCHOS_TEST_FLOATING_EQUALITY(target_data[2],3.0,tol,pout,success);
  TEUCHOS_TEST_FLOATING_EQUALITY(target_data[3],4.0,tol,pout,success);

  // Make sure all points were found.
  TEST_EQUALITY( volume_source_map.getMissedTargetPoints().size(), 0 );
}

//---------------------------------------------------------------------------//
// end tstVolumeSourceMap4.cpp
//---------------------------------------------------------------------------//
