//---------------------------------------------------------------------------//
/*
  Copyright (c) 2014, Stuart R. Slattery
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

  *: Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

  *: Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

  *: Neither the name of the Oak Ridge National Laboratory nor the
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
 * \file   tstIntrepidBasisFactory.cpp
 * \author Stuart R. Slattery
 * \brief  Intrepid basis factory unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <cstdlib>

#include <DTK_IntrepidBasisFactory.hpp>

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Ptr.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"

#include <Shards_CellTopology.hpp>
#include <Shards_BasicTopologies.hpp>

#include <Intrepid_Basis.hpp>
#include <Intrepid_FieldContainer.hpp>

//---------------------------------------------------------------------------//
// Tests.
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( IntrepidBasisFactory, factory_test )
{
    typedef Intrepid::FieldContainer<double> MDArray;
    typedef double Scalar;
    typedef Intrepid::Basis<Scalar,MDArray> Basis;

    // Line 2.
    {
	shards::CellTopology topo(
	    shards::getCellTopologyData<shards::Line<2> >() );
	Teuchos::RCP<Basis> basis = DataTransferKit::IntrepidBasisFactory::create( topo );
	TEST_EQUALITY( basis->getCardinality(), 2 );
	TEST_EQUALITY( basis->getDegree(), 1 );
    }

    // Triangle 3.
    {
	shards::CellTopology topo(
	    shards::getCellTopologyData<shards::Triangle<3> >() );
	Teuchos::RCP<Basis> basis = DataTransferKit::IntrepidBasisFactory::create( topo );
	TEST_EQUALITY( basis->getCardinality(), 3 );
	TEST_EQUALITY( basis->getDegree(), 1 );
    }

    // Triangle 6.
    {
	shards::CellTopology topo(
	    shards::getCellTopologyData<shards::Triangle<6> >() );
	Teuchos::RCP<Basis> basis = DataTransferKit::IntrepidBasisFactory::create( topo );
	TEST_EQUALITY( basis->getCardinality(), 6 );
	TEST_EQUALITY( basis->getDegree(), 2 );
    }

    // Quadrilateral 4.
    {
	shards::CellTopology topo(
	    shards::getCellTopologyData<shards::Quadrilateral<4> >() );
	Teuchos::RCP<Basis> basis = DataTransferKit::IntrepidBasisFactory::create( topo );
	TEST_EQUALITY( basis->getCardinality(), 4 );
	TEST_EQUALITY( basis->getDegree(), 1 );
    }

    // Quadrilateral 9.
    {
	shards::CellTopology topo(
	    shards::getCellTopologyData<shards::Quadrilateral<9> >() );
	Teuchos::RCP<Basis> basis = DataTransferKit::IntrepidBasisFactory::create( topo );
	TEST_EQUALITY( basis->getCardinality(), 9 );
	TEST_EQUALITY( basis->getDegree(), 2 );
    }

    // Tetrahedron 4.
    {
	shards::CellTopology topo(
	    shards::getCellTopologyData<shards::Tetrahedron<4> >() );
	Teuchos::RCP<Basis> basis = DataTransferKit::IntrepidBasisFactory::create( topo );
	TEST_EQUALITY( basis->getCardinality(), 4 );
	TEST_EQUALITY( basis->getDegree(), 1 );
    }

    // Tetrahedron 10.
    {
	shards::CellTopology topo(
	    shards::getCellTopologyData<shards::Tetrahedron<10> >() );
	Teuchos::RCP<Basis> basis = DataTransferKit::IntrepidBasisFactory::create( topo );
	TEST_EQUALITY( basis->getCardinality(), 10 );
	TEST_EQUALITY( basis->getDegree(), 2 );
    }

    // Hexahedron 8.
    {
	shards::CellTopology topo(
	    shards::getCellTopologyData<shards::Hexahedron<8> >() );
	Teuchos::RCP<Basis> basis = DataTransferKit::IntrepidBasisFactory::create( topo );
	TEST_EQUALITY( basis->getCardinality(), 8 );
	TEST_EQUALITY( basis->getDegree(), 1 );
    }

    // Hexahedron 27.
    {
	shards::CellTopology topo(
	    shards::getCellTopologyData<shards::Hexahedron<27> >() );
	Teuchos::RCP<Basis> basis = DataTransferKit::IntrepidBasisFactory::create( topo );
	TEST_EQUALITY( basis->getCardinality(), 27 );
	TEST_EQUALITY( basis->getDegree(), 2 );
    }

    // Wedge 6.
    {
	shards::CellTopology topo(
	    shards::getCellTopologyData<shards::Wedge<6> >() );
	Teuchos::RCP<Basis> basis = DataTransferKit::IntrepidBasisFactory::create( topo );
	TEST_EQUALITY( basis->getCardinality(), 6 );
	TEST_EQUALITY( basis->getDegree(), 1 );
    }

    // Wedge 18.
    {
	shards::CellTopology topo(
	    shards::getCellTopologyData<shards::Wedge<18> >() );
	Teuchos::RCP<Basis> basis = DataTransferKit::IntrepidBasisFactory::create( topo );
	TEST_EQUALITY( basis->getCardinality(), 18 );
	TEST_EQUALITY( basis->getDegree(), 2 );
    }

    // Pyramid 5.
    {
	shards::CellTopology topo(
	    shards::getCellTopologyData<shards::Pyramid<5> >() );
	Teuchos::RCP<Basis> basis = DataTransferKit::IntrepidBasisFactory::create( topo );
	TEST_EQUALITY( basis->getCardinality(), 5 );
	TEST_EQUALITY( basis->getDegree(), 1 );
    }

    // Pyramid 13.
    {
	shards::CellTopology topo(
	    shards::getCellTopologyData<shards::Pyramid<13> >() );
	Teuchos::RCP<Basis> basis = DataTransferKit::IntrepidBasisFactory::create( topo );
	TEST_EQUALITY( basis->getCardinality(), 5 );
	TEST_EQUALITY( basis->getDegree(), 1 );
    }

    // Pyramid 14.
    {
	shards::CellTopology topo(
	    shards::getCellTopologyData<shards::Pyramid<14> >() );
	Teuchos::RCP<Basis> basis = DataTransferKit::IntrepidBasisFactory::create( topo );
	TEST_EQUALITY( basis->getCardinality(), 5 );
	TEST_EQUALITY( basis->getDegree(), 1 );
    }
}

//---------------------------------------------------------------------------//
// end tstIntrepidBasisFactory.cpp
//---------------------------------------------------------------------------//
