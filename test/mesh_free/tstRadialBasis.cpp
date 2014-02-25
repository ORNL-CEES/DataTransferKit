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
 * \file   tstRadialBasis.cpp
 * \author Stuart R. Slattery
 * \brief  Basis function unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <stdexcept>

#include <DTK_RadialBasisPolicy.hpp>
#include <DTK_WendlandBasis.hpp>
#include <DTK_EuclideanDistance.hpp>

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"

//---------------------------------------------------------------------------//
// HELPER FUNCTIONS
//---------------------------------------------------------------------------//

// Get the default communicator.
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
// Tests.
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( WendlandBasis, dim_1_order_0 )
{
    typedef DataTransferKit::WendlandBasis<0> BasisType;
    typedef DataTransferKit::RadialBasisPolicy<BasisType> BP;

    int dim = 1;

    Teuchos::Array<double> x1(dim, 0.5);
    Teuchos::Array<double> x2(dim, 0.75);

    double radius_1 = 1.0;
    Teuchos::RCP<BasisType> basis_1 = BP::create( radius_1 );

    double dist = 
	DataTransferKit::EuclideanDistance<1>::distance( x1.getRawPtr(), x2.getRawPtr() );
    double basis_value = BP::evaluateValue( *basis_1, dist );
    double basis_grad = BP::evaluateGradient( *basis_1, dist );

    double x = 0.0;
    for ( int i = 0; i < dim; ++i )
    {
	x += (x2[i]-x1[i])*(x2[i]-x1[i]);
    }
    x = std::sqrt(x);
    double test_value = (1.0-x)*(1.0-x);
    double test_grad = 2.0*x-2.0;

    TEST_EQUALITY( test_value, basis_value );
    TEST_EQUALITY( test_grad, basis_grad );

    double radius_2 = 0.1;
    BasisType basis_2( radius_2 );

    basis_value = BP::evaluateValue( basis_2, dist );
    basis_grad = BP::evaluateGradient( basis_2, dist );

    TEST_EQUALITY( 0.0, basis_value );
    TEST_EQUALITY( 0.0, basis_grad );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( WendlandBasis, dim_2_order_0 )
{
    typedef DataTransferKit::WendlandBasis<0> BasisType;
    typedef DataTransferKit::RadialBasisPolicy<BasisType> BP;

    int dim = 2;

    Teuchos::Array<double> x1(dim, 0.5);
    Teuchos::Array<double> x2(dim, 0.75);

    double radius_1 = 1.0;
    Teuchos::RCP<BasisType> basis_1 = BP::create( radius_1 );

    double dist = 
	DataTransferKit::EuclideanDistance<2>::distance( x1.getRawPtr(), x2.getRawPtr() );
    double basis_value = BP::evaluateValue( *basis_1, dist );
    double basis_grad = BP::evaluateGradient( *basis_1, dist );

    double x = 0.0;
    for ( int i = 0; i < dim; ++i )
    {
	x += (x2[i]-x1[i])*(x2[i]-x1[i]);
    }
    x = std::sqrt(x);
    double test_value = (1.0-x)*(1.0-x);
    double test_grad = 2.0*x-2.0;

    TEST_EQUALITY( test_value, basis_value );
    TEST_EQUALITY( test_grad, basis_grad );

    double radius_2 = 0.1;
    BasisType basis_2( radius_2 );

    basis_value = BP::evaluateValue( basis_2, dist );
    basis_grad = BP::evaluateGradient( basis_2, dist );

    TEST_EQUALITY( 0.0, basis_value );
    TEST_EQUALITY( 0.0, basis_grad );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( WendlandBasis, dim_3_order_0 )
{
    typedef DataTransferKit::WendlandBasis<0> BasisType;
    typedef DataTransferKit::RadialBasisPolicy<BasisType> BP;

    int dim = 3;

    Teuchos::Array<double> x1(dim, 0.5);
    Teuchos::Array<double> x2(dim, 0.75);

    double radius_1 = 1.0;
    Teuchos::RCP<BasisType> basis_1 = BP::create( radius_1 );

    double dist = 
	DataTransferKit::EuclideanDistance<3>::distance( x1.getRawPtr(), x2.getRawPtr() );
    double basis_value = BP::evaluateValue( *basis_1, dist );
    double basis_grad = BP::evaluateGradient( *basis_1, dist );

    double x = 0.0;
    for ( int i = 0; i < dim; ++i )
    {
	x += (x2[i]-x1[i])*(x2[i]-x1[i]);
    }
    x = std::sqrt(x);
    double test_value = (1.0-x)*(1.0-x);
    double test_grad = 2.0*x-2.0;

    TEST_EQUALITY( test_value, basis_value );
    TEST_EQUALITY( test_grad, basis_grad );

    double radius_2 = 0.1;
    BasisType basis_2( radius_2 );

    basis_value = BP::evaluateValue( basis_2, dist );
    basis_grad = BP::evaluateGradient( basis_2, dist );

    TEST_EQUALITY( 0.0, basis_value );
    TEST_EQUALITY( 0.0, basis_grad );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( WendlandBasis, dim_1_order_2 )
{
    typedef DataTransferKit::WendlandBasis<2> BasisType;
    typedef DataTransferKit::RadialBasisPolicy<BasisType> BP;

    int dim = 1;

    Teuchos::Array<double> x1(dim, 0.5);
    Teuchos::Array<double> x2(dim, 0.75);

    double radius_1 = 1.0;
    Teuchos::RCP<BasisType> basis_1 = BP::create( radius_1 );

    double dist = 
	DataTransferKit::EuclideanDistance<1>::distance( x1.getRawPtr(), x2.getRawPtr() );
    double basis_value = BP::evaluateValue( *basis_1, dist );
    double basis_grad = BP::evaluateGradient( *basis_1, dist );

    double x = 0.0;
    for ( int i = 0; i < dim; ++i )
    {
	x += (x2[i]-x1[i])*(x2[i]-x1[i]);
    }
    x = std::sqrt(x);
    double test_value = (1.0-x)*(1.0-x)*(1.0-x)*(1.0-x)*(4.0*x+1.0);
    double test_grad = 20.0*x*(x-1.0)*(x-1.0)*(x-1.0);

    TEST_EQUALITY( test_value, basis_value );
    TEST_EQUALITY( test_grad, basis_grad );

    double radius_2 = 0.1;
    BasisType basis_2( radius_2 );

    basis_value = BP::evaluateValue( basis_2, dist );
    basis_grad = BP::evaluateGradient( basis_2, dist );

    TEST_EQUALITY( 0.0, basis_value );
    TEST_EQUALITY( 0.0, basis_grad );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( WendlandBasis, dim_2_order_2 )
{
    typedef DataTransferKit::WendlandBasis<2> BasisType;
    typedef DataTransferKit::RadialBasisPolicy<BasisType> BP;

    int dim = 2;

    Teuchos::Array<double> x1(dim, 0.5);
    Teuchos::Array<double> x2(dim, 0.75);

    double radius_1 = 1.0;
    Teuchos::RCP<BasisType> basis_1 = BP::create( radius_1 );

    double dist = 
	DataTransferKit::EuclideanDistance<2>::distance( x1.getRawPtr(), x2.getRawPtr() );
    double basis_value = BP::evaluateValue( *basis_1, dist );
    double basis_grad = BP::evaluateGradient( *basis_1, dist );

    double x = 0.0;
    for ( int i = 0; i < dim; ++i )
    {
	x += (x2[i]-x1[i])*(x2[i]-x1[i]);
    }
    x = std::sqrt(x);
    double test_value = (1.0-x)*(1.0-x)*(1.0-x)*(1.0-x)*(4.0*x+1.0);
    double test_grad = 20.0*x*(x-1.0)*(x-1.0)*(x-1.0);

    TEST_EQUALITY( test_value, basis_value );
    TEST_EQUALITY( test_grad, basis_grad );

    double radius_2 = 0.1;
    BasisType basis_2( radius_2 );

    basis_value = BP::evaluateValue( basis_2, dist );
    basis_grad = BP::evaluateGradient( basis_2, dist );

    TEST_EQUALITY( 0.0, basis_value );
    TEST_EQUALITY( 0.0, basis_grad );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( WendlandBasis, dim_3_order_2 )
{
    typedef DataTransferKit::WendlandBasis<2> BasisType;
    typedef DataTransferKit::RadialBasisPolicy<BasisType> BP;

    int dim = 3;

    Teuchos::Array<double> x1(dim, 0.5);
    Teuchos::Array<double> x2(dim, 0.75);

    double radius_1 = 1.0;
    Teuchos::RCP<BasisType> basis_1 = BP::create( radius_1 );

    double dist = 
	DataTransferKit::EuclideanDistance<3>::distance( x1.getRawPtr(), x2.getRawPtr() );
    double basis_value = BP::evaluateValue( *basis_1, dist );
    double basis_grad = BP::evaluateGradient( *basis_1, dist );

    double x = 0.0;
    for ( int i = 0; i < dim; ++i )
    {
	x += (x2[i]-x1[i])*(x2[i]-x1[i]);
    }
    x = std::sqrt(x);
    double test_value = (1.0-x)*(1.0-x)*(1.0-x)*(1.0-x)*(4.0*x+1.0);
    double test_grad = 20.0*x*(x-1.0)*(x-1.0)*(x-1.0);

    TEST_EQUALITY( test_value, basis_value );
    TEST_EQUALITY( test_grad, basis_grad );

    double radius_2 = 0.1;
    BasisType basis_2( radius_2 );

    basis_value = BP::evaluateValue( basis_2, dist );
    basis_grad = BP::evaluateGradient( basis_2, dist );

    TEST_EQUALITY( 0.0, basis_value );
    TEST_EQUALITY( 0.0, basis_grad );
}

//---------------------------------------------------------------------------//
// end tstRadialBasis.cpp
//---------------------------------------------------------------------------//
