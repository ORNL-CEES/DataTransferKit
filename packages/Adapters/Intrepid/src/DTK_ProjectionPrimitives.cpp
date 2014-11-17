//---------------------------------------------------------------------------//
/*
  Copyright (c) 2013, Stuart R. Slattery
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
 * \file   DTK_ProjectionPrimitives.cpp
 * \author Stuart Slattery
 * \brief  A stateless class of projection primitive operations.
 */
//---------------------------------------------------------------------------//

#include <limits>
#include <cmath>
#include <algorithm>

#include "DTK_DBC.hpp"
#include "DTK_ProjectionPrimitives.hpp"
#include "DTK_IntrepidBasisFactory.hpp"
#include "DTK_NewtonSolver.hpp"
#include "DTK_ProjectionPrimitiveNonlinearProblems.hpp"

#include <Teuchos_as.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_RCP.hpp>

#include <Shards_BasicTopologies.hpp>

#include <Intrepid_RealSpaceTools.hpp>

namespace DataTransferKit
{

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Get the center of the reference cell of the given topology.
 */
void ProjectionPrimitives::referenceCellCenter( 
    const shards::CellTopology& cell_topo, 
    Intrepid::FieldContainer<double>& cell_center )
{
    DTK_REQUIRE( 2 == cell_center.rank() );
    DTK_REQUIRE( Teuchos::as<unsigned>(cell_center.dimension(1)) == 
		   cell_topo.getDimension() );

    int num_cells = cell_center.dimension(0);    

    switch( cell_topo.getKey() ){
	case shards::Line<2>::key:
	case shards::Line<3>::key:
	    for ( int n = 0; n < num_cells; ++n )
	    {
		cell_center(n,0) = Teuchos::ScalarTraits<double>::zero();
	    }
	    break;
      
	case shards::Triangle<3>::key:
	case shards::Triangle<4>::key:
	case shards::Triangle<6>::key:    
	    for ( int n = 0; n < num_cells; ++n )
	    {

		cell_center(n,0) = 1.0/3.0;
		cell_center(n,1) = 1.0/3.0;  
	    }
	    break;
      
	case shards::Quadrilateral<4>::key:
	case shards::Quadrilateral<8>::key:
	case shards::Quadrilateral<9>::key:
	    for ( int n = 0; n < num_cells; ++n )
	    {
		cell_center(n,0) = Teuchos::ScalarTraits<double>::zero();      
		cell_center(n,1) = Teuchos::ScalarTraits<double>::zero();    
	    }
	    break;
      
	case shards::Tetrahedron<4>::key:
	case shards::Tetrahedron<10>::key:
	case shards::Tetrahedron<11>::key:
	    for ( int n = 0; n < num_cells; ++n )
	    {
		cell_center(n,0) = 1.0/6.0;    
		cell_center(n,1) = 1.0/6.0;    
		cell_center(n,2) = 1.0/6.0; 
	    }
	    break;
      
	case shards::Hexahedron<8>::key:
	case shards::Hexahedron<20>::key:
	case shards::Hexahedron<27>::key:
	    for ( int n = 0; n < num_cells; ++n )
	    {
		cell_center(n,0) = Teuchos::ScalarTraits<double>::zero();
		cell_center(n,1) = Teuchos::ScalarTraits<double>::zero();
		cell_center(n,2) = Teuchos::ScalarTraits<double>::zero();
	    }
	    break;
      
	case shards::Wedge<6>::key:
	case shards::Wedge<15>::key:
	case shards::Wedge<18>::key:
	    for ( int n = 0; n < num_cells; ++n )
	    {
		cell_center(n,0) = 1.0/3.0;
		cell_center(n,1) = 1.0/3.0;
		cell_center(n,2) = Teuchos::ScalarTraits<double>::zero();
	    }
	    break;

	case shards::Pyramid<5>::key:
	case shards::Pyramid<13>::key:
	case shards::Pyramid<14>::key:
	    for ( int n = 0; n < num_cells; ++n )
	    {
		cell_center(n,0) = Teuchos::ScalarTraits<double>::zero();
		cell_center(n,1) = Teuchos::ScalarTraits<double>::zero();
		cell_center(n,2) = 1.0/4.0;
	    }
	    break;

	default:
	    bool cell_topo_supported = false;
	    DTK_INSIST( cell_topo_supported );
	    break;
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Determine if a point is within the volume of influence of a face.
 *
 * \param parameters Projection parameters. 
 *
 * \param point Point coordinates (Dim).
 *
 * \param face_nodes Face node coordinates (Node,Dim).
 *
 * \param face_node_normals Normal vectors for each face node (Node,Dim).
 *
 * \param face_topology Cell topology of the face.
 *
 * \param max_newton_iters Maximum number of Newton iterations allowed.
 */
bool ProjectionPrimitives::pointInFaceVolumeOfInfluence(
    const Teuchos::ParameterList& parameters,
    const Intrepid::FieldContainer<double>& point,
    const Intrepid::FieldContainer<double>& face_nodes,
    const Intrepid::FieldContainer<double>& face_node_normals,
    const shards::CellTopology& face_topology )
{
    DTK_REQUIRE( 1 == point.rank() );
    DTK_REQUIRE( 2 == face_nodes.rank() );
    DTK_REQUIRE( 2 == face_node_normals.rank() );
    DTK_REQUIRE( point.dimension(0) == face_nodes.dimension(1) );
    DTK_REQUIRE( point.dimension(0) == face_node_normals.dimension(1) );
    DTK_REQUIRE( face_nodes.dimension(0) == face_node_normals.dimension(0) );
    DTK_REQUIRE( Teuchos::as<unsigned>(face_nodes.dimension(0)) == 
		   face_topology.getNodeCount() );

    // Get the spatial dimension
    int space_dim = point.dimension(0);
    DTK_CHECK( 3 == space_dim );

    // Get the geometric tolerance.
    double geometric_tolerance = parameters.get<double>("Geometric Tolerance");

    // Project the point onto each bilinear surface formed by the face vertex
    // normals and the face edges. If the point is on the correct face of each
    // surface then it is within the volume of influence of the
    // face. Counter-clockwise ordering of the nodes on the face about the
    // outward facing normal is required.
    Intrepid::FieldContainer<double> face_edge_nodes( 2, space_dim );
    Intrepid::FieldContainer<double> face_edge_node_normals( 2, space_dim );
    double distance_to_surface = 0.0;
    int face_num_edges = face_nodes.dimension(0);

    // First edges.
    for ( int e = 0; e < face_num_edges-1; ++e )
    {
	for ( int i = 0; i < space_dim; ++i )
	{
	    face_edge_nodes(0,i) = face_nodes(e,i);
	    face_edge_nodes(1,i) = face_nodes(e+1,i);
	    face_edge_node_normals(0,i) = face_node_normals(e,i);
	    face_edge_node_normals(1,i) = face_node_normals(e+1,i);
	}
	distance_to_surface = distanceToFaceBilinearSurface(
	    parameters, point, face_edge_nodes, face_edge_node_normals );
	if ( distance_to_surface > geometric_tolerance )
	{
	    return false;
	}
    }

    // Last edge.
    for ( int i = 0; i < space_dim; ++i )
    {
	face_edge_nodes(0,i) = face_nodes(face_num_edges-1,i);
	face_edge_nodes(1,i) = face_nodes(0,i);
	face_edge_node_normals(0,i) = face_node_normals(face_num_edges-1,i);
	face_edge_node_normals(1,i) = face_node_normals(0,i);
    }
    distance_to_surface = distanceToFaceBilinearSurface(
	parameters, point, face_edge_nodes, face_edge_node_normals );
    if ( distance_to_surface > geometric_tolerance )
    {
	return false;
    }

    // If we got here then the point is in the volume of influence of the
    // face as it was to the left of all bilinear surfaces formed by the
    // edges and their node normals.
    return true;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Project a point onto a face and return the physical
 * coordinates of the projected point on that face. This requires
 * the solution of a nonlinear parameterized projection problem.
 *
 * \param parameters Projection parameters. 
 *
 * \param point point coordinates (Dim)
 *
 * \param face_nodes face node coordinates (Node,Dim)
 *
 * \param face_node_normals face node normal vectors (Node,Dim)
 *
 * \param face_topology Cell topology of the face.
 *
 * \param projected_point Projected point in physical coordinates (Dim)
 *
 * \param face_edge_id The local id of the face edge onto which the node
 * projected. Will be -1 if the node did not project onto any of the face
 * edges.
 *
 * \param face_node_id The local id of the face node onto which the node
 * projected. Will be -1 if the node did not project onto any of the face
 * nodes.
 *
 * \return True if the point projected onto the face.
 */
void ProjectionPrimitives::projectPointToFace( 
    const Teuchos::ParameterList& parameters,
    const Intrepid::FieldContainer<double>& point,
    const Intrepid::FieldContainer<double>& face_nodes,
    const Intrepid::FieldContainer<double>& face_node_normals,
    const shards::CellTopology& face_topology,
    Intrepid::FieldContainer<double>& parametric_point,
    Intrepid::FieldContainer<double>& physical_point,
    int& face_edge_id,
    int& face_node_id )
{
    DTK_REQUIRE( 1 == point.rank() );
    DTK_REQUIRE( 2 == face_nodes.rank() );
    DTK_REQUIRE( 2 == face_node_normals.rank() );
    DTK_REQUIRE( 1 == physical_point.rank() );
    DTK_REQUIRE( point.dimension(0) == face_nodes.dimension(1) );
    DTK_REQUIRE( point.dimension(0) == face_node_normals.dimension(1) );
    DTK_REQUIRE( point.dimension(0) == physical_point.dimension(0) );
    DTK_REQUIRE( face_nodes.dimension(0) == 
		   face_node_normals.dimension(0) );
    DTK_REQUIRE( Teuchos::as<unsigned>(face_nodes.dimension(0)) == 
		   face_topology.getNodeCount() );
    DTK_REQUIRE( 2 == parametric_point.rank() );
    DTK_REQUIRE( 1 == parametric_point.dimension(0) );
    DTK_REQUIRE( point.dimension(0) == parametric_point.dimension(1) );


    // Get dimensions.
    int space_dim = point.dimension(0);
    int topo_dim = face_topology.getDimension();

    // Get the geometric tolerance.
    double geometric_tolerance = parameters.get<double>("Geometric Tolerance");

    // Get the newton tolerance.
    double newton_tolerance = parameters.get<double>("Newton Tolerance");

    // Get the max newton iterations.
    double max_newton_iters = parameters.get<int>("Max Newton Iterations");

    // Get the basis functions for the face cell topology.
    Teuchos::RCP<Intrepid::Basis<double,Intrepid::FieldContainer<double> > >
	face_basis = IntrepidBasisFactory::create( face_topology );

    // Initializeunknown vector.
    parametric_point.initialize( 0.0 );
   
    // Set the initial solution guess to the center of the cell for the
    // parametric coordinates and 0 for distance.
    int num_points = 1;
    Intrepid::FieldContainer<double> face_center( num_points, topo_dim );
    referenceCellCenter( face_topology, face_center );
    for ( int n = 0; n < topo_dim; ++n )
    {
	parametric_point(0,n) = face_center(0,n);
    }
    parametric_point(0,topo_dim) = Teuchos::ScalarTraits<double>::zero();

    // Build the nonlinear problem data.
    ProjectPointToFaceNonlinearProblem nonlinear_problem( 
	face_basis, point, face_nodes, face_node_normals );

    // Solve the nonlinear problem.
    NewtonSolver<ProjectPointToFaceNonlinearProblem>::solve(
	parametric_point, nonlinear_problem, newton_tolerance, max_newton_iters );

    // Apply tolerancing. If the point projected near a face edge or
    // node within the tolerance, move it to that point.
    if ( std::abs(parametric_point(0,0)) < geometric_tolerance )
    {
	parametric_point(0,0) = 0.0;
    }
    else if ( std::abs(1.0-parametric_point(0,0)) < geometric_tolerance )
    {
	parametric_point(0,0) = 1.0;
    }
    else if ( std::abs(1.0+parametric_point(0,0)) < geometric_tolerance )
    {
	parametric_point(0,0) = -1.0;
    }
    if ( std::abs(parametric_point(0,1)) < geometric_tolerance )
    {
	parametric_point(0,1) = 0.0;
    }
    else if ( std::abs(1.0-parametric_point(0,1)) < geometric_tolerance )
    {
	parametric_point(0,1) = 1.0;
    }
    else if ( std::abs(1.0+parametric_point(0,1)) < geometric_tolerance )
    {
	parametric_point(0,1) = -1.0;
    }
    // This case is for the diagonal of the triangle reference geometry.
    if ( std::abs(1.0-parametric_point(0,0)-parametric_point(0,1)) < 
	 geometric_tolerance )
    {
	parametric_point(0,0) = 1.0 - parametric_point(0,1);
    }

    // Compute the physical coordinates from the projection in the reference
    // space.
    nonlinear_problem.updateState( parametric_point );
    for ( int d = 0; d < space_dim; ++d )
    {
	physical_point(d) = 0.0;
	for ( int n = 0; n < nonlinear_problem.d_cardinality; ++n )
	{
	    physical_point(d) += 
		nonlinear_problem.d_basis_evals(n,0) * face_nodes(n,d);
	}
    }

    // Determine if the point projected onto any of the face nodes or edges.
    face_edge_id = -1;
    face_node_id = -1;

    // Triangle case.
    if ( (shards::Triangle<3>::key == face_topology.getKey()) ||
	 (shards::Triangle<4>::key == face_topology.getKey()) ||
	 (shards::Triangle<6>::key == face_topology.getKey()) )
    {
	if ( (0.0 == parametric_point(0,0)) && (0.0 == parametric_point(0,1)) )
	{
	    face_node_id = 0;
	}
	else if ( (1.0 == parametric_point(0,0)) && (0.0 == parametric_point(0,1)) )
	{
	    face_node_id = 1;
	}
	else if ( (0.0 == parametric_point(0,0)) && (1.0 == parametric_point(0,1)) )
	{
	    face_node_id = 2;
	}
	else if ( 0.0 == parametric_point(0,1) )
	{
	    face_edge_id = 0;
	}
	else if ( 1.0 == parametric_point(0,0) + parametric_point(0,1) )
	{
	    face_edge_id = 1;
	}
	else if ( 0.0 == parametric_point(0,0) )
	{
	    face_edge_id = 2;
	}
    }
    // Quadrilateral case.
    else if ( (shards::Quadrilateral<4>::key == face_topology.getKey()) ||
	      (shards::Quadrilateral<8>::key == face_topology.getKey()) ||
	      (shards::Quadrilateral<9>::key == face_topology.getKey()) )
    {
	if ( (-1.0 == parametric_point(0,0)) && (-1.0 == parametric_point(0,1)) )
	{
	    face_node_id = 0;
	}
	else if ( (1.0 == parametric_point(0,0)) && (-1.0 == parametric_point(0,1)) )
	{
	    face_node_id = 1;
	}
	else if ( (1.0 == parametric_point(0,0)) && (1.0 == parametric_point(0,1)) )
	{
	    face_node_id = 2;
	}
	else if ( (-1.0 == parametric_point(0,0)) && (1.0 == parametric_point(0,1)) )
	{
	    face_node_id = 3;
	}
	else if( -1.0 == parametric_point(0,1) )
	{
	    face_edge_id = 0;
	}
	else if( 1.0 == parametric_point(0,0) )
	{
	    face_edge_id = 1;
	}
	else if( 1.0 == parametric_point(0,1) )
	{
	    face_edge_id = 2;
	}
	else if( -1.0 == parametric_point(0,0) )
	{
	    face_edge_id = 3;
	}
    }
    // Unsupported case.
    else
    {
	DTK_INSIST(
	    (shards::Triangle<3>::key == face_topology.getKey()) ||
	    (shards::Triangle<4>::key == face_topology.getKey()) ||
	    (shards::Triangle<6>::key == face_topology.getKey()) ||
	    (shards::Quadrilateral<4>::key == face_topology.getKey()) ||
	    (shards::Quadrilateral<8>::key == face_topology.getKey()) ||
	    (shards::Quadrilateral<9>::key == face_topology.getKey()) );
    }
}

//---------------------------------------------------------------------------//
/*
 * \brief Project a feature point to a feature edge. Return false if
 * the point did not project onto the edge.
 *
 * \param parameters Projection parameters. 
 *
 * \param point Point coordinates (Dim)
 *
 * \param edge_nodes Edge node coordinates (Node,Dim)
 *
 * \param edge_node_normals Edge node normal vectors (Node,Dim)
 *
 * \param projected_point Projected point in physical coordinates (Dim)
 *
 * \param edge_node_id The local id of the edge node onto which the
 * feature node projected. Will be -1 if the feature node did not
 * project onto any of the edge nodes.
 *
 * \return Return true of the feature node projected onto the feature edge.
 */
bool ProjectionPrimitives::projectPointFeatureToEdgeFeature(
    const Teuchos::ParameterList& parameters,
    const Intrepid::FieldContainer<double>& point,
    const Intrepid::FieldContainer<double>& point_normal,
    const Intrepid::FieldContainer<double>& edge_nodes,
    const Intrepid::FieldContainer<double>& edge_node_normals,
    Intrepid::FieldContainer<double>& projected_point,
    int& edge_node_id )
{
    DTK_REQUIRE( 1 == point.rank() );
    DTK_REQUIRE( 1 == point_normal.rank() );
    DTK_REQUIRE( 2 == edge_nodes.rank() );
    DTK_REQUIRE( 2 == edge_node_normals.rank() );
    DTK_REQUIRE( 1 == projected_point.rank() );
    DTK_REQUIRE( point.dimension(0) == edge_nodes.dimension(1) );
    DTK_REQUIRE( point.dimension(0) == edge_node_normals.dimension(1) );
    DTK_REQUIRE( point.dimension(0) == projected_point.dimension(0) );
    DTK_REQUIRE( edge_nodes.dimension(0) == 
		   edge_node_normals.dimension(0) );
    DTK_REQUIRE( edge_nodes.dimension(0) == 2 );

    // Get dimensions.
    int space_dim = point.dimension(0);

    // Get the geometric tolerance.
    double geometric_tolerance = parameters.get<double>("Geometric Tolerance");

    // Get the normal tolerance.
    double normal_tolerance = parameters.get<double>("Normal Tolerance");

    // Get the newton tolerance.
    double newton_tolerance = parameters.get<double>("Newton Tolerance");

    // Get the max newton iterations.
    double max_newton_iters = parameters.get<int>("Max Newton Iterations");

    // Check to make sure that the projection direction is opposite the node
    // normal directions. 
    double dot_1 = 0.0;
    double dot_2 = 0.0;
    for ( int d = 0; d < space_dim; ++d )
    {
	dot_1 += point_normal(d)*edge_node_normals(0,d);
	dot_2 += point_normal(d)*edge_node_normals(1,d);
    }
    if ( dot_1 > -0.5*normal_tolerance || dot_2 > -0.5*normal_tolerance )
    {
	return false;
    }

    // Build the binormals.
    Intrepid::FieldContainer<double> gvec( space_dim );
    for ( int d = 0; d < space_dim; ++d )
    {
	gvec(d) = edge_nodes(1,d) - edge_nodes(0,d);
    }
    Intrepid::FieldContainer<double> normal( space_dim );
    Intrepid::FieldContainer<double> l( space_dim );
    Intrepid::FieldContainer<double> edge_node_binormals( 
	edge_node_normals.dimension(0),
	edge_node_normals.dimension(1) );
    for ( int n = 0; n < 2; ++n )
    {
	for ( int d = 0; d < space_dim; ++d )
	{
	    normal(d) = edge_node_normals(n,d);
	}
	Intrepid::RealSpaceTools<double>::vecprod( l, normal, gvec );
	for ( int d = 0; d < space_dim; ++d )
	{
	    edge_node_binormals(n,d) = l(d);
	}
    }

    // Allocate unknown vector.
    int num_points = 1;
    Intrepid::FieldContainer<double> u( num_points, space_dim );
   
    // Set the initial solution guess to 0.5.
    for ( int n = 0; n < space_dim; ++n )
    {
	u(0,n) = 0.5;
    }

    // Build the nonlinear problem data.
    ProjectPointFeatureToEdgeFeatureNonlinearProblem
	nonlinear_problem( point, edge_nodes, 
			   edge_node_normals, edge_node_binormals );

    // Solve the nonlinear problem.
    NewtonSolver<ProjectPointFeatureToEdgeFeatureNonlinearProblem>::solve(
	u, nonlinear_problem, newton_tolerance, max_newton_iters );

    // Apply tolerancing. 
    if ( std::abs(u(0,0)) < geometric_tolerance )
    {
	u(0,0) = 0.0;
    }
    else if ( std::abs(1.0-u(0,0)) < geometric_tolerance )
    {
	u(0,0) = 1.0;
    }

    // See if the solution was outside the parameteric bounds of the edge.
    if ( 0.0 > u(0,0) || 1.0 < u(0,0) )
    {
	return false;
    }

    // Compute the physical coordinates from the projection in the reference
    // space.
    for ( int d = 0; d < space_dim; ++d )
    {
	projected_point(d) = edge_nodes(0,d) + u(0,0)*gvec(d);
    }
    
    // Check to see if the point projected onto one of the vertices of the
    // edge.
    edge_node_id = -1;
    if ( u(0,0) == 0.0 )
    {
	edge_node_id = 0;
    }
    else if ( u(0,0) == 1.0 )
    {
	edge_node_id = 1;
    }

    // There was an interesection so return true.
    return true;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Intersect two edges in 3 dimensions and return their intersection
 * point realized on both edges. Return the intersection type.
 *
 * \param parameters Intersection parameters. 
 *
 * \param edge_1 edge 1 node coordinates (node,dim).
 *
 * \param edge_2 edge 2 node coordinates (node,dim).
 *
 * \param edge_2_node_normals edge 2 node normal vectors (node,dim).
 *
 * \param edge_1_intersection The coordinates of the intersection realized on
 * the edge 1 (dim).
 *
 * \param edge_2_intersection The coordinates of the intersection realized on
 * the edge 2 (dim).
 *
 * \param edge_1_node_id The local id of the edge 1 node onto which the edge 2
 * nodes projected. Will be -1 if the edge 2 did not intersect with any of
 * the edge 1 nodes.
 *
 * \param edge_2_node_id The local id of the edge 2 node onto which the edge 1
 * nodes projected. Will be -1 if the edge 1 did not intersect with any of the
 * edge 2 nodes.
 *
 * \return True if there was an intersection.
*/
bool ProjectionPrimitives::edgeEdgeIntersection( 
    const Teuchos::ParameterList& parameters,
    const Intrepid::FieldContainer<double>& edge_1,
    const Intrepid::FieldContainer<double>& edge_2,
    const Intrepid::FieldContainer<double>& edge_2_node_normals,
    Intrepid::FieldContainer<double>& edge_1_intersection,
    Intrepid::FieldContainer<double>& edge_2_intersection,
    int& edge_1_node_id,
    int& edge_2_node_id )
{
    DTK_REQUIRE( 2 == edge_1.rank() );
    DTK_REQUIRE( 2 == edge_2.rank() );
    DTK_REQUIRE( 2 == edge_2_node_normals.rank() );
    DTK_REQUIRE( 1 == edge_1_intersection.rank() );
    DTK_REQUIRE( 1 == edge_2_intersection.rank() );
    DTK_REQUIRE( 2 == edge_1.dimension(0) );
    DTK_REQUIRE( 2 == edge_2.dimension(0) );
    DTK_REQUIRE( 2 == edge_2_node_normals.dimension(0) );
    DTK_REQUIRE( edge_1.dimension(1) == edge_2.dimension(1) );
    DTK_REQUIRE( edge_1.dimension(1) == edge_2_node_normals.dimension(1) );
    DTK_REQUIRE( edge_1.dimension(1) == edge_1_intersection.dimension(0) );
    DTK_REQUIRE( edge_1.dimension(1) == edge_2_intersection.dimension(0) );

    // The condition number tolerance.
    double cond_tol = 1.0 / std::sqrt( std::numeric_limits<double>::epsilon() );

    // Get the spatial dimension.
    int space_dim = edge_1.dimension(1);

    // Get the geometric tolerance.
    double geometric_tolerance = parameters.get<double>("Geometric Tolerance");

    // Get the merge tolerance.
    double merge_tolerance = parameters.get<double>("Merge Tolerance");

    // Build intermediate vectors.
    Intrepid::FieldContainer<double> bvec( space_dim );
    Intrepid::FieldContainer<double> gvec( space_dim );
    Intrepid::FieldContainer<double> dvec( space_dim );
    Intrepid::FieldContainer<double> d0( space_dim );
    Intrepid::FieldContainer<double> g0minusb0( space_dim );
    for ( int i = 0; i < space_dim; ++i )
    {
	bvec(i) = edge_1(1,i) - edge_1(0,i);
	gvec(i) = edge_2(1,i) - edge_2(0,i);
	dvec(i) = edge_2_node_normals(1,i) - edge_2_node_normals(0,i);
	d0(i) = edge_2_node_normals(0,i);
	g0minusb0(i) = edge_2(0,i) - edge_1(0,i);
    }

    // Intersection parametric coordinates.
    double alpha = 0.0;
    double beta = 0.0;

    // Compute the parameterized realization of the intersection on edge
    // 2. Start by building the quadratic equation coefficients.
    Intrepid::FieldContainer<double> work_vec( space_dim );
    Intrepid::RealSpaceTools<double>::vecprod( work_vec, bvec, dvec );
    double a = Intrepid::RealSpaceTools<double>::dot( work_vec, gvec );
    double b = Intrepid::RealSpaceTools<double>::dot( work_vec, g0minusb0 );
    Intrepid::RealSpaceTools<double>::vecprod( work_vec, bvec, d0 );
    b += Intrepid::RealSpaceTools<double>::dot( work_vec, gvec );
    double c = Intrepid::RealSpaceTools<double>::dot( work_vec, g0minusb0 );

    // Solve the quadratic equation if a != 0 to get the parameterized
    // realization of the intersection on the edge 2.
    if ( std::abs(a) > geometric_tolerance )
    {
	double ax2 = 2*a;
	double sqrt_arg = b*b - 4*a*c;

	if ( std::abs(sqrt_arg) < geometric_tolerance )
	{
	    sqrt_arg = 0.0;
	}

	// If we get a negative argument for the square root then beta
	// will be in non-physical (complex) coordinates.
	if ( sqrt_arg < 0.0 )
	{
	    return false;
	}

	// Compute beta.
	double sqrt_val = std::sqrt(sqrt_arg);
	double beta_1 = (-b + sqrt_val) / (ax2);
	double beta_2 = (-b - sqrt_val) / (ax2);

	// Resolve inconsistencies by perturbing small beta onto vertices.
	// These small values will result from the inexact primitive
	// solutions.
	if ( std::abs(beta_1) < merge_tolerance )
	{
	    beta_1 = 0.0;
	}
	else if ( std::abs(1.0 - beta_1) < merge_tolerance )
	{
	    beta_1 = 1.0;
	}
	if ( std::abs(beta_2) < merge_tolerance )
	{
	    beta_2 = 0.0;
	}
	else if ( std::abs(1.0 - beta_2) < merge_tolerance )
	{
	    beta_2 = 1.0;
	}

	// If there are two valid solutions then return no intersection as
	// these edges are nearly parallel.
	if ( (beta_1 >= 0.0 && beta_1 <= 1.0) &&
	     (beta_2 >= 0.0 && beta_2 <= 1.0) )
	{
	    return false;
	}

	// Choose the solution that is in the accepted range.
	if ( beta_1 >= 0.0 && beta_1 <= 1.0 )
	{
	    beta = beta_1;
	}
	else if ( beta_2 >= 0.0 && beta_2 <= 1.0 )
	{
	    beta = beta_2;
	}
	else
	{
	    // There is no intersection if neither solution to the quadratic
	    // equation is between 0 and 1.
	    return false;
	}
    }
    // If a = 0 and b != 0 then the equation is linear.
    else if ( std::abs(b) > geometric_tolerance )
    {
	beta = -c / b;
	
	// Resolve inconsistencies by perturbing small beta onto vertices.
	// These small values will result from the inexact primitive
	// solutions.
	if ( std::abs(beta) < merge_tolerance )
	{
	    beta = 0.0;
	}
	else if ( std::abs(1.0 - beta) < merge_tolerance )
	{
	    beta = 1.0;
	}
    }
    // If a = 0 and b = 0 then beta is undefined.
    else
    {
	return false;
    }

    // Check the condition number of beta.
    Intrepid::FieldContainer<double> nvec( space_dim );
    double n_length = 0.0;
    double g_length = 0.0;
    for ( int i = 0; i < space_dim; ++i )
    {
	nvec(i) = g0minusb0(i) + beta*gvec(i);
	n_length += nvec(i)*nvec(i);
	g_length += gvec(i)*gvec(i);
    }
    Intrepid::FieldContainer<double> ncrossg( space_dim );
    Intrepid::RealSpaceTools<double>::vecprod( ncrossg, nvec, gvec );
    double ncrossg_length = 0.0;
    for ( int i = 0; i < space_dim; ++i )
    {
	ncrossg_length += ncrossg(i)*ncrossg(i);
    }

    // Only check the condition number if it is defined.
    if ( ncrossg_length > geometric_tolerance &&
	 n_length > geometric_tolerance )
    {
	double beta_cond = std::sqrt(n_length) * std::sqrt(g_length) / 
			   std::sqrt(ncrossg_length);
	if ( beta_cond > cond_tol )
	{
	    return false;
	}
    }

    // Compute the parameterized realization of the intersection on the edge
    // 1.
    Intrepid::FieldContainer<double> l( space_dim );
    Intrepid::RealSpaceTools<double>::scale( dvec, beta );
    Intrepid::RealSpaceTools<double>::add( work_vec, d0, dvec );
    Intrepid::RealSpaceTools<double>::vecprod( l, gvec, work_vec );
    double dot1 = Intrepid::RealSpaceTools<double>::dot( l, g0minusb0 );
    double dot2 = Intrepid::RealSpaceTools<double>::dot( l, bvec );

    if ( std::abs(dot2) < geometric_tolerance )
    {
	return false;
    }
    alpha = dot1 / dot2;

    // Resolve inconsistencies by pertubing small alpha onto vertices. These
    // small values will result from the inexact primitive solutions.
    if ( std::abs(alpha) < merge_tolerance )
    {
	alpha = 0.0;
    }
    else if ( std::abs(1.0 - alpha) < merge_tolerance )
    {
	alpha = 1.0;
    }

    // Check the condition number of alpha.
    Intrepid::FieldContainer<double> lvec( space_dim );
    double l_length = 0.0;
    double b_length = 0.0;
    for ( int i = 0; i < space_dim; ++i )
    {
	lvec(i) = alpha*bvec(i) - g0minusb0(i);
	l_length += lvec(i)*lvec(i);
	b_length += bvec(i)*bvec(i);
    }
    Intrepid::FieldContainer<double> lcrossb( space_dim );
    Intrepid::RealSpaceTools<double>::vecprod( lcrossb, lvec, bvec );
    double lcrossb_length = 0.0;
    for ( int i = 0; i < space_dim; ++i )
    {
	lcrossb_length += lcrossb(i)*lcrossb(i);
    }

    // Only check the condition number if it is defined.
    if ( lcrossb_length > geometric_tolerance && l_length > geometric_tolerance )
    {
	double alpha_cond = std::sqrt(l_length) * std::sqrt(b_length) / 
			    std::sqrt(lcrossb_length);
	if ( alpha_cond > cond_tol )
	{
	    return false;
	}
    }

    // Check for an intersection on the edge 1.
    if ( alpha < 0.0 || alpha > 1.0 )
    {
	return false;
    }

    // Check for an intersection on the edge 2.
    if ( beta < 0.0 || beta > 1.0 )
    {
	return false;
    }

    // Build the physical intersection location on the edge 2.
    for ( int i = 0; i < space_dim; ++i )
    {
	edge_2_intersection(i) = edge_2(0,i) + beta*gvec(i);
    }

    // Build the physical intersection location on the edge 1.
    for ( int i = 0; i < space_dim; ++i )
    {
	edge_1_intersection(i) = edge_1(0,i) + alpha*bvec(i);
    }

    // Check the error in the solution.
    double d_length = 0.0;
    Intrepid::FieldContainer<double> v_vec( space_dim );
    for ( int i = 0; i < space_dim; ++i )
    {
	d_length += work_vec(i)*work_vec(i);
	v_vec(i) = edge_1_intersection(i) - edge_2_intersection(i);
    }
    d_length *= -d_length;
    d_length += 1.0;
    for ( int i = 0; i < space_dim; ++i )
    {
	v_vec(i) *= d_length;
    }
    double v_length = 0.0;
    for ( int i = 0; i < space_dim; ++i )
    {
	v_length += v_vec(i)*v_vec(i);
    }
    double error = std::sqrt(v_length) / 
		   std::min( std::sqrt(b_length), std::sqrt(g_length) );
    if ( error > cond_tol || 
	 b_length < geometric_tolerance ||
	 g_length < geometric_tolerance )
    {
	return false;
    }

    // Determine if any of the edge 1 nodes were intersected.
    if ( 0.0 == alpha )
    {
	edge_1_node_id = 0;
    }
    else if ( 1.0 == alpha )
    {
	edge_1_node_id = 1;
    }
    else
    {
	edge_1_node_id = -1;
    }

    // Determine if any of the edge 2 nodes were intersected.
    if ( 0.0 == beta )
    {
	edge_2_node_id = 0;
    }
    else if ( 1.0 == beta )
    {
	edge_2_node_id = 1;
    }
    else
    {
	edge_2_node_id = -1;
    }

    // Both alpha and beta are between 0 and 1 so there is an
    // intersection.
    return true;
}

//---------------------------------------------------------------------------//
// PRIVATE IMPLEMENTATION
//---------------------------------------------------------------------------//
/*!
 * \brief Compute the distance of a projected point onto a bilinear surface
 * formed by a face edge and its normals.
 *
 * \param parameters Intersection parameters. 
 *
 * \param point Point coordinates (Dim).
 *
 * \param face_edge_nodes Face node edge coordinates (Node,Dim).
 *
 * \param face_edge_node_normals Normal vectors for each face edge node
 * (Node,Dim).
 */
typename Intrepid::FieldContainer<double>::scalar_type 
ProjectionPrimitives::distanceToFaceBilinearSurface(
    const Teuchos::ParameterList& parameters,
    const Intrepid::FieldContainer<double>& point,
    const Intrepid::FieldContainer<double>& face_edge_nodes,
    const Intrepid::FieldContainer<double>& face_edge_node_normals )
{
    // Get the geometric tolerance.
    double geometric_tolerance = parameters.get<double>("Geometric Tolerance");

    // Get the newton tolerance.
    double newton_tolerance = parameters.get<double>("Newton Tolerance");

    // Get the max newton iterations.
    double max_newton_iters = parameters.get<int>("Max Newton Iterations");

    // Compute the scale factor. Choose half the length of the face for
    // simplicity.
    double xdist = face_edge_nodes(1,0) - face_edge_nodes(0,0);
    double ydist = face_edge_nodes(1,1) - face_edge_nodes(0,1);
    double zdist = face_edge_nodes(1,2) - face_edge_nodes(0,2);
    double c = 0.5*std::sqrt( xdist*xdist + ydist*ydist + zdist*zdist );

    // Get the distance to the bilinear surface.
    Intrepid::FieldContainer<double> u( 1, point.dimension(0) );
    u(0,0) = 0.5;
    u(0,1) = 0.5;
    u(0,2) = Teuchos::ScalarTraits<double>::zero();
    PointInFaceVolumeOfInfluenceNonlinearProblem nonlinear_problem(
	point, face_edge_nodes, face_edge_node_normals, c );
    NewtonSolver<PointInFaceVolumeOfInfluenceNonlinearProblem>::solve(
	u, nonlinear_problem, newton_tolerance, max_newton_iters );

    // Check for degeneracy.
    nonlinear_problem.updateState( u );
    Intrepid::FieldContainer<double> v1( nonlinear_problem.d_space_dim );
    Intrepid::FieldContainer<double> v2( nonlinear_problem.d_space_dim );
    for ( int i = 0; i < nonlinear_problem.d_space_dim; ++i )
    {
	v1(i) = nonlinear_problem.d_face_edge_nodes(1,i) - 
		nonlinear_problem.d_face_edge_nodes(0,i);
	v2(i) = v1(i) + nonlinear_problem.d_c*u(0,1)*( 
	    nonlinear_problem.d_face_edge_node_normals(1,i) - 
	    nonlinear_problem.d_face_edge_node_normals(0,i) );
    }

    // Return a positive number if degenerate.
    if ( Intrepid::RealSpaceTools<double>::dot(v1,v2) < geometric_tolerance )
    {
	return 1.0;
    }

    // Return the solution if not degenerate.
    return u(0,2);
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_ProjectionPrimitives.cpp
//---------------------------------------------------------------------------//
