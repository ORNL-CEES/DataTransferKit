//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   tstEntityCenteredDOFVector.cpp
 * \author Stuart Slattery
 * \date   Wed May 25 12:36:14 2011
 * \brief  Entity-centered DOF vector test.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include "DTK_EntityCenteredDOFVector.hpp"

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Ptr.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_DefaultComm.hpp"

#include <Thyra_MultiVectorBase.hpp>
#include <Thyra_MultiVectorStdOps.hpp>

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( EntityCenteredDOFVector, vector_test )
{
    // Initialize parallel communication.
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
	Teuchos::DefaultComm<int>::getComm();

    // Vector parameters.
    int num_vec = 3;
    int vec_length = 10;

    // Create data.
    Teuchos::Array<std::size_t> ids( vec_length );
    Teuchos::ArrayRCP<double> in_data( num_vec*vec_length );
    Teuchos::ArrayRCP<double> out_data( num_vec*vec_length );
    for ( int i = 0; i < vec_length; ++i )
    {
	ids[i] = i+1;
	in_data[i] = 2.0*i + 1.0;
	out_data[i] = 0.0;
    }

    // Create an input vector.
    Teuchos::RCP<Thyra::MultiVectorBase<double> > in_vec =
	DataTransferKit::EntityCenteredDOFVector::createThyraMultiVector( 
	    comm, ids(), in_data, vec_length, num_vec );
        
    // Create an output vector.
    Teuchos::RCP<Thyra::MultiVectorBase<double> > out_vec =
	DataTransferKit::EntityCenteredDOFVector::createThyraMultiVector( 
	    comm, ids(), out_data, vec_length, num_vec );

    // Add the vectors together.
    Thyra::Vp_V( Teuchos::ptr(out_vec.getRawPtr()), *in_vec );

    // Check the results.
    for ( int i = 0; i < vec_length; ++i )
    {
	TEST_EQUALITY( in_data[i], out_data[i] );
    }
}

//---------------------------------------------------------------------------//
// end of tstEntityCenteredDOFVector.cpp
//---------------------------------------------------------------------------//
