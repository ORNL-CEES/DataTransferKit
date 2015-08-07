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
 * \file   DTK_ProjectionPrimitiveNonlinearProblems_impl.hpp
 * \author Stuart Slattery
 * \brief  Nonlinear problem definitions for projection primitives.
 */
//---------------------------------------------------------------------------//

#include <cmath>

#include "DTK_DBC.hpp"
#include "DTK_ProjectionPrimitiveNonlinearProblems.hpp"

#include <Teuchos_RCP.hpp>

#include <Shards_CellTopology.hpp>

#include <Intrepid_Types.hpp>
#include <Intrepid_Basis.hpp>
#include <Intrepid_RealSpaceTools.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// ProjectPointToFaceNonlinearProblem Implementation.
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
ProjectPointToFaceNonlinearProblem::
ProjectPointToFaceNonlinearProblem( 
    const Teuchos::RCP<
    Intrepid::Basis<Scalar,Intrepid::FieldContainer<double> > >& face_basis,
    const Intrepid::FieldContainer<double>& point,
    const Intrepid::FieldContainer<double>& face_nodes,
    const Intrepid::FieldContainer<double>& face_node_normals )
    : d_face_basis( face_basis )
    , d_point( point )
    , d_face_nodes( face_nodes )
    , d_face_node_normals( face_node_normals )
{
    d_space_dim = d_point.dimension(0);
    d_topo_dim = d_space_dim - 1;
    d_cardinality = d_face_basis->getCardinality();
    int num_points = 1;
    d_basis_evals = Intrepid::FieldContainer<double>( 
	d_cardinality, num_points );
    d_grad_evals = Intrepid::FieldContainer<double>( 
	d_cardinality, num_points, d_topo_dim );
    d_eval_points = Intrepid::FieldContainer<double>( num_points, d_topo_dim );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Update the state of the problem given the new solution vector.
 */
void ProjectPointToFaceNonlinearProblem::updateState( 
    const Intrepid::FieldContainer<double>& u )
{
    // Extract the current natural coordinates from the solution vector.
    for ( int i = 0; i < d_topo_dim; ++i )
    {
	d_eval_points(0,i) = u(0,i);
    }

    // Evaluate the basis at the current natural coordinates.
    d_face_basis->getValues( 
	d_basis_evals, d_eval_points, Intrepid::OPERATOR_VALUE );

    // Evaluate the basis gradient at the current natural coordinates.
    d_face_basis->getValues( 
	d_grad_evals, d_eval_points, Intrepid::OPERATOR_GRAD );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Evaluate the nonlinear residual.
 */
void ProjectPointToFaceNonlinearProblem::evaluateResidual( 
    const Intrepid::FieldContainer<double>& u, 
    Intrepid::FieldContainer<double>& F ) const
{
    DTK_REQUIRE( 2 == u.rank() );
    DTK_REQUIRE( 2 == F.rank() );
    DTK_REQUIRE( F.dimension(0) == u.dimension(0) );
    DTK_REQUIRE( F.dimension(1) == u.dimension(1) );

    // Build the nonlinear residual.
    for ( int i = 0; i < d_space_dim; ++i )
    {
	F(0,i) = -d_point(i);
	for ( int n = 0; n < d_cardinality; ++n )
	{
	    F(0,i) += d_basis_evals(n,0) * ( 
		d_face_nodes(n,i)
		- u(0,d_topo_dim)*d_face_node_normals(n,i) 
		);
	}
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Evaluate the jacobian.
 */
void ProjectPointToFaceNonlinearProblem::evaluateJacobian( 
    const Intrepid::FieldContainer<double>& u, 
    Intrepid::FieldContainer<double>& J ) const
{
    DTK_REQUIRE( 2 == u.rank() );
    DTK_REQUIRE( 3 == J.rank() );
    DTK_REQUIRE( J.dimension(0) == u.dimension(0) );
    DTK_REQUIRE( J.dimension(1) == u.dimension(1) );
    DTK_REQUIRE( J.dimension(2) == u.dimension(1) );

    // Build the Jacobian. 
    for ( int i = 0; i < d_space_dim; ++i )
    {
	// Start with the basis gradient contributions.
	for ( int j = 0; j < d_topo_dim; ++j )
	{
	    J(0,i,j) = 0.0;
	    for ( int n = 0; n < d_cardinality; ++n )
	    {
		J(0,i,j) += d_grad_evals(n,0,j) * (
		    d_face_nodes(n,i) -
		    u(0,d_topo_dim)*d_face_node_normals(n,i)
		    );
	    }
	}

	// Then add the point normal contribution.
	J(0,i,d_space_dim-1) = 0.0;
	for ( int n = 0; n < d_cardinality; ++n )
	{
	    J(0,i,d_space_dim-1) -= 
		d_basis_evals(n,0) * d_face_node_normals(n,i); 
	}
    }
}

//---------------------------------------------------------------------------//
// PointInFaceVolumeOfInfluenceNonlinearProblem Implementation.
//---------------------------------------------------------------------------//
/*! 
 * \brief Constructor.
 */
PointInFaceVolumeOfInfluenceNonlinearProblem::
PointInFaceVolumeOfInfluenceNonlinearProblem( 
    const Intrepid::FieldContainer<double>& point,
    const Intrepid::FieldContainer<double>& face_edge_nodes,
    const Intrepid::FieldContainer<double>& face_edge_node_normals,
    const double c )
    : d_point( point )
    , d_face_edge_nodes( face_edge_nodes )
    , d_face_edge_node_normals( face_edge_node_normals )
    , d_c( c )
{
    DTK_CHECK( d_point.dimension(0) == 
		 d_face_edge_nodes.dimension(1) );
    DTK_CHECK( d_point.dimension(0) == 
		 d_face_edge_node_normals.dimension(1) );
    DTK_CHECK( d_face_edge_nodes.dimension(0) ==
		 d_face_edge_node_normals.dimension(0) );
    d_space_dim = point.dimension(0);
    d_topo_dim = d_space_dim - 1;
    d_cardinality = 4;
    int num_points = 1;
    d_basis_evals = Intrepid::FieldContainer<double>( 
	d_cardinality, num_points );
    d_grad_evals = Intrepid::FieldContainer<double>( 
	d_cardinality, num_points, d_topo_dim );
    d_eval_points = Intrepid::FieldContainer<double>( num_points, d_topo_dim );
    d_proj_normal = Intrepid::FieldContainer<double>( d_space_dim );

    // Compute the two extra psuedo-nodes along the node normal
    // directions.
    d_face_normal_nodes = Intrepid::FieldContainer<double>( 2, d_space_dim );
    for ( int j = 0; j < d_space_dim; ++j )
    {
	d_face_normal_nodes(0,j) = 
	    face_edge_nodes(0,j) + d_c*d_face_edge_node_normals(0,j);
	d_face_normal_nodes(1,j) = 
	    face_edge_nodes(1,j) + d_c*d_face_edge_node_normals(1,j);
    }
}


//---------------------------------------------------------------------------//
/*!
 * \brief Update the state of the problem given the new solution vector.
 */
void PointInFaceVolumeOfInfluenceNonlinearProblem::updateState( 
    const Intrepid::FieldContainer<double>& u )
{
    // Evaluate the basis at the current natural coordinates.
    d_basis_evals(0,0) = (1.0 - u(0,0))*(1.0 - u(0,1));
    d_basis_evals(1,0) = u(0,0)*(1.0 - u(0,1));
    d_basis_evals(2,0) = u(0,0)*u(0,1);
    d_basis_evals(3,0) = (1.0 - u(0,0))*u(0,1);

    // Evaluate the basis gradient at the current natural coordinates.
    d_grad_evals(0,0,0) = u(0,1) - 1.0;
    d_grad_evals(1,0,0) = 1.0 - u(0,1);
    d_grad_evals(2,0,0) = u(0,1);
    d_grad_evals(3,0,0) = -u(0,1);
    d_grad_evals(0,0,1) = u(0,0) - 1.0;
    d_grad_evals(1,0,1) = -u(0,0);
    d_grad_evals(2,0,1) = u(0,0);
    d_grad_evals(3,0,1) = 1.0 - u(0,0);

    // Compute the current projection normal direction.
    Intrepid::FieldContainer<double> v1( d_space_dim );
    Intrepid::FieldContainer<double> v2( d_space_dim );
    for ( int i = 0; i < d_space_dim; ++i )
    {
	v1(i) = d_face_edge_nodes(1,i) - d_face_edge_nodes(0,i) +
		d_c*u(0,1)*( d_face_edge_node_normals(1,i) - 
			     d_face_edge_node_normals(0,i) );
	v2(i) = u(0,0)*d_face_edge_node_normals(1,i) +
		(1-u(0,0))*d_face_edge_node_normals(0,i);
    }
    Intrepid::RealSpaceTools<Scalar>::vecprod(
	d_proj_normal, v1, v2 );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Evaluate the nonlinear residual.
 */
void PointInFaceVolumeOfInfluenceNonlinearProblem::evaluateResidual( 
    const Intrepid::FieldContainer<double>& u, 
    Intrepid::FieldContainer<double>& F ) const
{
    DTK_REQUIRE( 2 == u.rank() );
    DTK_REQUIRE( 2 == F.rank() );
    DTK_REQUIRE( F.dimension(0) == u.dimension(0) );
    DTK_REQUIRE( F.dimension(1) == u.dimension(1) );

    // Build the nonlinear residual.
    for ( int i = 0; i < d_space_dim; ++i )
    {
	F(0,i) = d_basis_evals(0,0)*d_face_edge_nodes(0,i) +
		 d_basis_evals(1,0)*d_face_edge_nodes(1,i) +
		 d_basis_evals(2,0)*d_face_normal_nodes(1,i) +
		 d_basis_evals(3,0)*d_face_normal_nodes(0,i) +
		 u(0,2)*d_proj_normal(i) -
		 d_point(i);
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Evaluate the jacobian.
 */
void PointInFaceVolumeOfInfluenceNonlinearProblem::evaluateJacobian( 
    const Intrepid::FieldContainer<double>& u, 
    Intrepid::FieldContainer<double>& J ) const
{
    DTK_REQUIRE( 2 == u.rank() );
    DTK_REQUIRE( 3 == J.rank() );
    DTK_REQUIRE( J.dimension(0) == u.dimension(0) );
    DTK_REQUIRE( J.dimension(1) == u.dimension(1) );
    DTK_REQUIRE( J.dimension(2) == u.dimension(1) );

    // Build the Jacobian. 
    for ( int i = 0; i < d_space_dim; ++i )
    {
	// Start with the basis gradient contributions.
	for ( int j = 0; j < d_topo_dim; ++j )
	{
	    J(0,i,j) = d_grad_evals(0,0,j) * d_face_edge_nodes(0,i) +
		       d_grad_evals(1,0,j) * d_face_edge_nodes(1,i) +
		       d_grad_evals(2,0,j) * d_face_normal_nodes(1,i) +
		       d_grad_evals(3,0,j) * d_face_normal_nodes(0,i);
	}

	// Then add the point normal contribution.
	J(0,i,d_space_dim-1) = d_proj_normal(i);
    }
}

//---------------------------------------------------------------------------//
// ProjectPointFeatureToEdgeFeatureNonlinearProblem Implementation.
//---------------------------------------------------------------------------//
/*! 
 * \brief Constructor.
 */
ProjectPointFeatureToEdgeFeatureNonlinearProblem::
ProjectPointFeatureToEdgeFeatureNonlinearProblem( 
	const Intrepid::FieldContainer<double>& point,
	const Intrepid::FieldContainer<double>& edge_nodes,
	const Intrepid::FieldContainer<double>& edge_node_normals,
	const Intrepid::FieldContainer<double>& edge_node_binormals )
    : d_point( point )
    , d_edge_nodes( edge_nodes )
    , d_edge_node_normals( edge_node_normals )
    , d_edge_node_binormals( edge_node_binormals )
{
    DTK_CHECK( 1 == d_point.rank() );
    DTK_CHECK( 2 == d_edge_nodes.rank() );
    DTK_CHECK( 2 == d_edge_node_normals.rank() );
    DTK_CHECK( 2 == d_edge_node_binormals.rank() );
    DTK_CHECK( d_point.dimension(0) == 
		 d_edge_nodes.dimension(1) );
    DTK_CHECK( d_point.dimension(0) == 
		 d_edge_node_normals.dimension(1) );
    DTK_CHECK( d_edge_nodes.dimension(0) ==
		 d_edge_node_normals.dimension(0) );
    DTK_CHECK( d_point.dimension(0) == 
		 d_edge_node_binormals.dimension(1) );
    DTK_CHECK( d_edge_nodes.dimension(0) ==
		 d_edge_node_binormals.dimension(0) );
    d_space_dim = d_point.dimension(0);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Update the state of the problem given the new solution vector.
 */
void 
ProjectPointFeatureToEdgeFeatureNonlinearProblem::updateState( 
    const Intrepid::FieldContainer<double>& u )
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Evaluate the nonlinear residual.
 */
void 
ProjectPointFeatureToEdgeFeatureNonlinearProblem::evaluateResidual(
    const Intrepid::FieldContainer<double>& u, 
    Intrepid::FieldContainer<double>& F ) const
{
    DTK_REQUIRE( 2 == u.rank() );
    DTK_REQUIRE( 2 == F.rank() );
    DTK_REQUIRE( F.dimension(0) == u.dimension(0) );
    DTK_REQUIRE( F.dimension(1) == u.dimension(1) );

    // Build the nonlinear residual.
    for ( int i = 0; i < d_space_dim; ++i )
    {
	F(0,i) = d_edge_nodes(0,i) +
		 u(0,0)*(d_edge_nodes(1,i)-d_edge_nodes(0,i)) +
		 u(0,1)*( d_edge_node_normals(0,i) +
			  u(0,0)*(d_edge_node_normals(1,i)-
				  d_edge_node_normals(0,i)) ) +
		 u(0,2)*( d_edge_node_binormals(0,i) +
			  u(0,0)*(d_edge_node_binormals(1,i)-
				  d_edge_node_binormals(0,i)) ) -
		 d_point(i);
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Evaluate the jacobian.
 */
void 
ProjectPointFeatureToEdgeFeatureNonlinearProblem::evaluateJacobian(
    const Intrepid::FieldContainer<double>& u, 
    Intrepid::FieldContainer<double>& J ) const
{
    DTK_REQUIRE( 2 == u.rank() );
    DTK_REQUIRE( 3 == J.rank() );
    DTK_REQUIRE( J.dimension(0) == u.dimension(0) );
    DTK_REQUIRE( J.dimension(1) == u.dimension(1) );
    DTK_REQUIRE( J.dimension(2) == u.dimension(1) );

    // Build the Jacobian. 
    for ( int i = 0; i < d_space_dim; ++i )
    {
	J(0,i,0) = (d_edge_nodes(1,i)-d_edge_nodes(0,i)) +
		   u(0,1)*( d_edge_node_normals(1,i)-
			    d_edge_node_normals(0,i) ) +
		   u(0,2)*( d_edge_node_binormals(1,i)-
			    d_edge_node_binormals(0,i) );

	J(0,i,1) = d_edge_node_normals(0,i) +
		   u(0,0)*(d_edge_node_normals(1,i)-
			   d_edge_node_normals(0,i));

	J(0,i,2) = d_edge_node_binormals(0,i) +
		   u(0,0)*(d_edge_node_binormals(1,i)-
			   d_edge_node_binormals(0,i));
    }
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_ProjectionPrimitiveNonlinearProblems.cpp
//---------------------------------------------------------------------------//

