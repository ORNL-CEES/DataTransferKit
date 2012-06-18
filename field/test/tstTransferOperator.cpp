//---------------------------------------------------------------------------//
/*!
 * \file tstTransferOperator.cpp
 * \author Stuart R. Slattery
 * \brief Transfer operator unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <DTK_DataSource.hpp>
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
// Unit tests
//---------------------------------------------------------------------------//

TEUCHOS_UNIT_TEST( TransferOperator, transfer_operator_test )
{
    // Setup communication.
    Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm<int>();

    // Setup source mesh.
    MyMesh source_mesh = buildMyMesh();

    // Setup target coordinate field.
    Teuchos::Array<double> target_coords = buildCoordinateField();

    // Create field evaluators.
    Teuchos::RCP< FieldEvaluator<MyMesh,Teuchos::Array<double> > 
		  temperature_evaluator = Teuchos::rcp(
		      new  TemperatureEvaluator() );

    Teuchos::RCP< FieldEvaluator<MyMesh,Teuchos::Array<float> > 
		  pressure_evaluator = Teuchos::rcp(
		      new  PressureEvaluator() );

    // Create a map for a consistent interpolation scheme.
    ConsistentInterpolation map( comm );

    // Setup and apply the transfer operator to the fields.
    TransferOperator<ConsistentInterpolation> transfer_operator( map );
    transfer_operator.setup( source_mesh, target_coords );
    transfer_operator.apply( temperature_evaluator, temperature_target );
    transfer_operator.apply( pressure_evaluator, pressure_target );
}

//---------------------------------------------------------------------------//
// end tstTransferOperator.cpp
//---------------------------------------------------------------------------//

