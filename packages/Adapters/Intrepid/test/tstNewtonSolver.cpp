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
 * \file   tstNewtonSolver.cpp
 * \author Stuart R. Slattery
 * \brief  Newton solver unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <cstdlib>

#include <DTK_NewtonSolver.hpp>
#include <DTK_NonlinearProblemTraits.hpp>

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Ptr.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"

#include "Intrepid_FieldContainer.hpp"

//---------------------------------------------------------------------------//
// HELPER FUNCTIONS
//---------------------------------------------------------------------------//

// Get the default communicator.
Teuchos::RCP<const Teuchos::Comm<int> > getDefaultComm()
{
#ifdef HAVE_MPI
    return Teuchos::DefaultComm<int>::getComm();
#else
    return Teuchos::rcp(new Teuchos::SerialComm<int>() );
#endif
}

//---------------------------------------------------------------------------//
// Scalar nonlinear problem. u^2 + u = a
//---------------------------------------------------------------------------//
class ScalarNonlinearProblem
{
  public:

    typedef Intrepid::FieldContainer<double> md_array_type;
    typedef md_array_type::scalar_type       scalar_type;

    ScalarNonlinearProblem( const double a )
	: d_a( a )
    { /* ... */ }

    ~ScalarNonlinearProblem() { /* ... */ }

    void evaluateResidual( const md_array_type& u, md_array_type& F ) const
    {
	F(0,0) = u(0,0)*u(0,0) + u(0,0) - d_a;
    }

    void evaluateJacobian( const md_array_type& u, md_array_type& J ) const
    {
	J(0,0,0) = 2.0*u(0,0) + 1;
    }

  private:
    
    double d_a;
};

// Traits implementation.
namespace DataTransferKit
{
template<>
class NonlinearProblemTraits<ScalarNonlinearProblem>
{
  public:

    typedef ScalarNonlinearProblem nonlinear_problem_type;
    typedef typename nonlinear_problem_type::md_array_type MDArray;
    typedef typename nonlinear_problem_type::scalar_type Scalar;

    static inline void updateState( 
	ScalarNonlinearProblem& problem, const MDArray& u )
    { /* ... */ }

    static inline void evaluateResidual( const ScalarNonlinearProblem& problem,
					 const MDArray& u, 
					 MDArray& F )
    {
	problem.evaluateResidual( u, F );
    }

    static inline void evaluateJacobian( const ScalarNonlinearProblem& problem,
					 const MDArray& u, 
					 MDArray& J )
    {
	problem.evaluateJacobian( u, J );
    }
};
}

//---------------------------------------------------------------------------//
// Vector nonlinear problem. u_1^2 + u_2 = a, u_1*u_2^2 - u_1 = b
//---------------------------------------------------------------------------//
class VectorNonlinearProblem
{
  public:

    typedef Intrepid::FieldContainer<double> md_array_type;
    typedef md_array_type::scalar_type       scalar_type;

    VectorNonlinearProblem( const double a, const double b )
	: d_a( a )
	, d_b( b )
    { /* ... */ }

    ~VectorNonlinearProblem() { /* ... */ }

    void evaluateResidual( const md_array_type& u, md_array_type& F ) const
    {
	F(0,0) = u(0,0)*u(0,0) + u(0,1) - d_a;
	F(0,1) = u(0,0)*(u(0,1)*u(0,1) - 1.0) - d_b;
    }

    void evaluateJacobian( const md_array_type& u, md_array_type& J ) const
    {
	J(0,0,0) = 2.0*u(0,0);
	J(0,0,1) = 1.0;
	J(0,1,0) = u(0,1)*u(0,1) - 1.0;
	J(0,1,1) = 2.0*u(0,0)*u(0,1);
    }

  private:
    
    double d_a;
    double d_b;
};

// Traits implementation.
namespace DataTransferKit
{
template<>
class NonlinearProblemTraits<VectorNonlinearProblem>
{
  public:

    typedef VectorNonlinearProblem nonlinear_problem_type;
    typedef typename nonlinear_problem_type::md_array_type MDArray;
    typedef typename nonlinear_problem_type::scalar_type Scalar;

    static inline void updateState( 
	VectorNonlinearProblem& problem, const MDArray& u )
    { /* ... */ }

    static inline void evaluateResidual( const VectorNonlinearProblem& problem,
					 const MDArray& u, 
					 MDArray& F )
    {
	problem.evaluateResidual( u, F );
    }

    static inline void evaluateJacobian( const VectorNonlinearProblem& problem,
					 const MDArray& u, 
					 MDArray& J )
    {
	problem.evaluateJacobian( u, J );
    }
};
}

//---------------------------------------------------------------------------//
// Tests.
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( NewtonSolver, scalar_test )
{
    DataTransferKit::NewtonSolver<ScalarNonlinearProblem> solver;
    double tolerance = 1.0e-12;
    double max_iters = 10000;

    double a_1 = 1.0;
    ScalarNonlinearProblem problem_1( a_1 );
    Intrepid::FieldContainer<double> u_1( 1, 1 );
    solver.solve( u_1, problem_1, tolerance, max_iters );
    TEST_FLOATING_EQUALITY( u_1(0,0)*u_1(0,0) + u_1(0,0), a_1, tolerance );

    double a_2 = 1830.31;
    ScalarNonlinearProblem problem_2( a_2 );
    Intrepid::FieldContainer<double> u_2( 1, 1 );
    solver.solve( u_2, problem_2, tolerance, max_iters );
    TEST_FLOATING_EQUALITY( u_2(0,0)*u_2(0,0) + u_2(0,0), a_2, tolerance );

    double a_3 = 0.0000231;
    ScalarNonlinearProblem problem_3( a_3 );
    Intrepid::FieldContainer<double> u_3( 1, 1 );
    solver.solve( u_3, problem_3, tolerance, max_iters );
    TEST_FLOATING_EQUALITY( u_3(0,0)*u_3(0,0) + u_3(0,0), a_3, tolerance );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( NewtonSolver, vector_test )
{
    DataTransferKit::NewtonSolver<VectorNonlinearProblem> solver;
    double tolerance = 1.0e-12;
    double max_iters = 10000;

    double a_1 = 1.0;
    double b_1 = 1.0;
    VectorNonlinearProblem problem_1( a_1, b_1 );
    Intrepid::FieldContainer<double> u_1( 1, 2 );
    solver.solve( u_1, problem_1, tolerance, max_iters );
    TEST_FLOATING_EQUALITY( u_1(0,0)*u_1(0,0) + u_1(0,1), a_1, tolerance );
    TEST_FLOATING_EQUALITY(
	u_1(0,0)*u_1(0,1)*u_1(0,1) - u_1(0,0), b_1, tolerance );

    double a_2 = 483.20;
    double b_2 = 0.32;
    VectorNonlinearProblem problem_2( a_2, b_2 );
    Intrepid::FieldContainer<double> u_2( 1, 2 );
    solver.solve( u_2, problem_2, tolerance, max_iters );
    TEST_FLOATING_EQUALITY( u_2(0,0)*u_2(0,0) + u_2(0,1), a_2, tolerance );
    TEST_FLOATING_EQUALITY(
	u_2(0,0)*u_2(0,1)*u_2(0,1) - u_2(0,0), b_2, tolerance );

    double a_3 = 0.00322;
    double b_3 = 1.987;
    VectorNonlinearProblem problem_3( a_3, b_3 );
    Intrepid::FieldContainer<double> u_3( 1, 2 );
    solver.solve( u_3, problem_3, tolerance, max_iters );
    TEST_FLOATING_EQUALITY( u_3(0,0)*u_3(0,0) + u_3(0,1), a_3, tolerance );
    TEST_FLOATING_EQUALITY(
	u_3(0,0)*u_3(0,1)*u_3(0,1) - u_3(0,0), b_3, tolerance );
}

//---------------------------------------------------------------------------//
// end tstNewtonSolver.cpp
//---------------------------------------------------------------------------//
