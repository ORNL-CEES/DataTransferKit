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
  LIMITED TO, PROCUREMENT OF SUBSTITUTE std::size_tODS OR SERVICES; LOSS OF USE,
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
//---------------------------------------------------------------------------//
/*!
 * \file   tstSplineCoefficientMatrix.cpp
 * \author Stuart R. Slattery
 * \brief  Spline coeffcient matrix.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <cstdlib>

#include <DTK_SplineCoefficientMatrix.hpp>
#include <DTK_SplineInterpolationPairing.hpp>
#include <DTK_WuBasis.hpp>

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Ptr.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"

#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>

#include "BelosPseudoBlockGmresSolMgr.hpp"
#include "BelosLinearProblem.hpp"

#include <Thyra_TpetraThyraWrappers.hpp>
#include <Thyra_DefaultMultipliedLinearOp.hpp>
#include <Thyra_DefaultScaledAdjointLinearOp.hpp>
#include <Thyra_DefaultAddedLinearOp.hpp>

#include <BelosThyraAdapter.hpp>

//---------------------------------------------------------------------------//
// Tests.
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( SplineCoefficientMatrix, coeff_mtx_test )
{
    // Initialize.
    Teuchos::RCP<const Teuchos::Comm<int> > comm = 
	Teuchos::DefaultComm<int>::getComm();
    const int dim = 3;
    double radius = 1.1;
    int num_src_points = 5;
    int num_src_coords = dim*num_src_points;
    int comm_rank = comm->getRank();
    int comm_size = comm->getSize();
    int offset = (0==comm_rank) ? dim + 1 : 0;

    // Create some coordinates.
    Teuchos::Array<double> src_coords(num_src_coords);
    Teuchos::Array<std::size_t> src_gids(num_src_points+offset);
    for ( int i = 0; i < num_src_points; ++i )
    {
	src_gids[i+offset] = comm_rank*num_src_points + i;
	src_coords[dim*i] = double(std::rand())/RAND_MAX;
	src_coords[dim*i+1] = double(std::rand())/RAND_MAX;
	src_coords[dim*i+2] = double(std::rand())/RAND_MAX;
    }
    for ( int i = 0; i < offset; ++i )
    {
	src_gids[i] = num_src_points*comm_size + i;
    }

    // Create a map.
    Teuchos::RCP<const Tpetra::Map<int,std::size_t> > map =
	Tpetra::createNonContigMap<int,std::size_t>( src_gids(), comm );

    // Create a pairing.
    DataTransferKit::SplineInterpolationPairing<dim> pairing( 
	src_coords(), src_coords(), radius );

    // Create the coefficient matrices.
    DataTransferKit::WuBasis<4> basis( radius );
    DataTransferKit::SplineCoefficientMatrix<DataTransferKit::WuBasis<4>,dim>
	coeff_mtx( map, src_coords(), src_gids(offset,num_src_points), 
		   src_coords(), src_gids(offset,num_src_points),
		   pairing, basis );

    // Create an abstract wrapper for P.
    Teuchos::RCP<const Thyra::VectorSpaceBase<double> > thyra_range_vector_space_P =
    	Thyra::createVectorSpace<double>( coeff_mtx.getP()->getRangeMap() );
    Teuchos::RCP<const Thyra::VectorSpaceBase<double> > thyra_domain_vector_space_P =
    	Thyra::createVectorSpace<double>( coeff_mtx.getP()->getDomainMap() );
    Teuchos::RCP<const Thyra::TpetraLinearOp<double,int,std::size_t> > thyra_P =
    	Teuchos::rcp( new Thyra::TpetraLinearOp<double,int,std::size_t>() );
    Teuchos::rcp_const_cast<Thyra::TpetraLinearOp<double,int,std::size_t> >(
	thyra_P)->constInitialize( 
	    thyra_range_vector_space_P, thyra_domain_vector_space_P, coeff_mtx.getP() );

    // Create an abstract wrapper for M.
    Teuchos::RCP<const Thyra::VectorSpaceBase<double> > thyra_range_vector_space_M =
    	Thyra::createVectorSpace<double>( coeff_mtx.getM()->getRangeMap() );
    Teuchos::RCP<const Thyra::VectorSpaceBase<double> > thyra_domain_vector_space_M =
    	Thyra::createVectorSpace<double>( coeff_mtx.getM()->getDomainMap() );
    Teuchos::RCP<const Thyra::TpetraLinearOp<double,int,std::size_t> > thyra_M =
    	Teuchos::rcp( new Thyra::TpetraLinearOp<double,int,std::size_t>() );
    Teuchos::rcp_const_cast<Thyra::TpetraLinearOp<double,int,std::size_t> >(
	thyra_M)->constInitialize( 
	    thyra_range_vector_space_M, thyra_domain_vector_space_M, coeff_mtx.getM() );

    // Create a transpose of P.
    Teuchos::RCP<const Thyra::LinearOpBase<double> > thyra_P_T =
	Thyra::transpose<double>( thyra_P );

    // Create a composite operator C = (P + M + P^T)
    Teuchos::RCP<const Thyra::LinearOpBase<double> > thyra_PpM =
	Thyra::add<double>( thyra_P, thyra_M );
    Teuchos::RCP<const Thyra::LinearOpBase<double> > thyra_C =
	Thyra::add<double>( thyra_PpM, thyra_P_T );

    // Create some vectors.
    int num_vec = 1;
    Teuchos::RCP<Tpetra::MultiVector<double,int,std::size_t> > x =
	Tpetra::createMultiVector<double,int,std::size_t>( map, num_vec );
    Teuchos::RCP<Tpetra::MultiVector<double,int,std::size_t> > y =
	Tpetra::createMultiVector<double,int,std::size_t>( map, num_vec );
    Teuchos::RCP<Thyra::MultiVectorBase<double> > thyra_x =
	Thyra::createMultiVector( x );
    Teuchos::RCP<Thyra::MultiVectorBase<double> > thyra_y =
	Thyra::createMultiVector( y );

    // Apply C to the vectors.
    x->putScalar(1.0);
    thyra_P->apply( Thyra::NOTRANS, *thyra_x, thyra_y.ptr(), 1.0, 0.0 );
    thyra_M->apply( Thyra::NOTRANS, *thyra_x, thyra_y.ptr(), 1.0, 0.0 );
    thyra_P_T->apply( Thyra::NOTRANS, *thyra_x, thyra_y.ptr(), 1.0, 0.0 );
    thyra_PpM->apply( Thyra::NOTRANS, *thyra_x, thyra_y.ptr(), 1.0, 0.0 );
    thyra_C->apply( Thyra::NOTRANS, *thyra_x, thyra_y.ptr(), 1.0, 0.0 );

    // Setup a belos solver.
    Teuchos::RCP<Teuchos::ParameterList> parameters = Teuchos::parameterList();
    parameters->set("Verbosity",
		    Belos::Errors|Belos::Warnings|
		    Belos::IterationDetails|Belos::OrthoDetails|
		    Belos::FinalSummary|Belos::TimingDetails|
		    Belos::StatusTestDetails|Belos::Debug);
    parameters->set("Output Frequency",1);
    typedef Thyra::MultiVectorBase<double> MV;
    typedef Thyra::LinearOpBase<double> OP;
    x->randomize();
    y->randomize();
    Belos::LinearProblem<double,MV,OP> problem( thyra_C, thyra_x, thyra_y );
    problem.setProblem();
    Belos::PseudoBlockGmresSolMgr<double,MV,OP> solver( 
	Teuchos::rcpFromRef(problem), parameters );
    solver.getCurrentParameters()->print();

    // Solve the problem with belos.
    solver.solve();
}

//---------------------------------------------------------------------------//
// end tstSplineCoefficientMatrix.cpp
//---------------------------------------------------------------------------//
