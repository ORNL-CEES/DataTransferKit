//---------------------------------------------------------------------------//
/*
  Copyright (c) 2012, Stuart R. Slattery
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

  *: Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

  *: Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

  *: Neither the name of the University of Wisconsin - Madison nor the
  names of its contributors may be used to endorse or promote products
  derived from this software without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
//---------------------------------------------------------------------------//
/*! 
 * \file tstL2ProjectionOperator.cpp
 * \author Stuart R. Slattery
 * \brief L2ProjectionOperator unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <cassert>

#include "reference_implementation/DTK_ReferenceHexMesh.hpp"

#include <DTK_FieldMultiVector.hpp>

#include <DTK_L2ProjectionOperator.hpp>

#include <DTK_BasicEntityPredicates.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ParameterList.hpp>

//---------------------------------------------------------------------------//
// TEST EPSILON
//---------------------------------------------------------------------------//

const double epsilon = 1.0e-7;

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( L2ProjectionOperator, l2_projection )
{
    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
	Teuchos::DefaultComm<int>::getComm();
    DataTransferKit::LocalEntityPredicate local_pred( comm->getRank() );

    // Set the global problem bounds.
    double x_max = 3.1;
    double y_max = 5.2;
    double z_max = 8.3;
    
    // Create a target mesh and field.
    int num_sx = 8;
    int num_sy = 8;
    int num_sz = 8;
    double sx_width = x_max / (num_sx-1);
    double sy_width = y_max / (num_sy-1);
    double sz_width = z_max / (num_sz-1);    
    Teuchos::Array<double> sx( num_sx );
    Teuchos::Array<double> sy( num_sy );
    Teuchos::Array<double> sz( num_sz );
    for ( int i = 0; i < num_sx; ++i ) sx[i] = i*sx_width;
    for ( int i = 0; i < num_sy; ++i ) sy[i] = i*sy_width;
    for ( int i = 0; i < num_sz; ++i ) sz[i] = i*sz_width;

    DataTransferKit::UnitTest::ReferenceHexMesh source_mesh( comm, sx, sy, sz );
    auto source_field = source_mesh.nodalField( 1 );
    Teuchos::RCP<DataTransferKit::FieldMultiVector> source_vector =
        Teuchos::rcp( new DataTransferKit::FieldMultiVector(comm, source_field) );

    // Put some data on the source field.
    double data_x = 9.3;
    double data_y = 2.2;
    double data_z = 1.33;
    auto source_nodes = source_mesh.functionSpace()->entitySet()->entityIterator(
        0, local_pred.getFunction() );
    auto source_nodes_begin = source_nodes.begin();
    auto source_nodes_end = source_nodes.end();
    for ( source_nodes = source_nodes.begin();
          source_nodes != source_nodes.end();
          ++source_nodes )
    {
        unsigned k = source_nodes->id() / (num_sx*num_sy);
        unsigned j = (source_nodes->id() - k*num_sx*num_sy) / num_sx;
        unsigned i = source_nodes->id() - j*num_sx - k*num_sx*num_sy;
        TEST_EQUALITY( source_nodes->id(), i + j*num_sx + k*num_sx*num_sy );

        double data = sx[i]*data_x + sy[j]*data_y + sz[k]*data_z + 1.0;
        source_field->writeFieldData( source_nodes->id(), 0, data );
    }

    // Create a source mesh and field.
    int num_tx = 9;
    int num_ty = 7;
    int num_tz = 7;
    double tx_width = x_max / (num_tx-1);
    double ty_width = y_max / (num_ty-1);
    double tz_width = z_max / (num_tz-1);    
    Teuchos::Array<double> tx( num_tx );
    Teuchos::Array<double> ty( num_ty );
    Teuchos::Array<double> tz( num_tz );
    for ( int i = 0; i < num_tx; ++i ) tx[i] = i*tx_width;
    for ( int i = 0; i < num_ty; ++i ) ty[i] = i*ty_width;
    for ( int i = 0; i < num_tz; ++i ) tz[i] = i*tz_width;

    DataTransferKit::UnitTest::ReferenceHexMesh target_mesh( comm, tx, ty, tz );
    auto target_field = target_mesh.nodalField( 1 );
    Teuchos::RCP<DataTransferKit::FieldMultiVector> target_vector =
        Teuchos::rcp( new DataTransferKit::FieldMultiVector(comm, target_field) );
    
    // MAPPING
    // Create a map.
    Teuchos::RCP<Teuchos::ParameterList> parameters = Teuchos::parameterList();
    Teuchos::ParameterList& l2_list = parameters->sublist("L2 Projection");
    l2_list.set("Integration Order",3);
    Teuchos::ParameterList& search_list = parameters->sublist("Search");
    search_list.set("Point Inclusion Tolerance",1.0e-6);
    
    Teuchos::RCP<DataTransferKit::L2ProjectionOperator> map_op = Teuchos::rcp(
	new DataTransferKit::L2ProjectionOperator(
	    source_vector->getMap(),target_vector->getMap(),*parameters) );

    // Setup the map.
    map_op->setup(
	source_mesh.functionSpace(), target_mesh.functionSpace() );

    // Apply the map.
    map_op->apply( *source_vector, *target_vector );

    // Check the results of the mapping.
    auto target_nodes = target_mesh.functionSpace()->entitySet()->entityIterator(
        0, local_pred.getFunction() );
    auto target_nodes_begin = target_nodes.begin();
    auto target_nodes_end = target_nodes.end();
    for ( target_nodes = target_nodes.begin();
          target_nodes != target_nodes.end();
          ++target_nodes )
    {
        unsigned k = target_nodes->id() / (num_tx*num_ty);
        unsigned j = (target_nodes->id() - k*num_tx*num_ty) / num_tx;
        unsigned i = target_nodes->id() - j*num_tx - k*num_tx*num_ty;
        TEST_EQUALITY( target_nodes->id(), i + j*num_tx + k*num_tx*num_ty );

        double gold_data = tx[i]*data_x + ty[j]*data_y + tz[k]*data_z + 1.0;
        double target_data = target_field->readFieldData(target_nodes->id(), 0);
        TEST_FLOATING_EQUALITY( target_data, gold_data, epsilon );
    }
}

//---------------------------------------------------------------------------//
// end tstL2ProjectionOperator.cpp
//---------------------------------------------------------------------------//
