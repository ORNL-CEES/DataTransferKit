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
 * \file   DTK_ProjectionPrimitiveNonlinearProblems.hpp
 * \author Stuart Slattery
 * \brief  Nonlinear problem definitions for projection primitives.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_PROJECTIONPRIMITIVENONLINEARPROBLEMS_HPP
#define DTK_PROJECTIONPRIMITIVENONLINEARPROBLEMS_HPP

#include <cmath>

#include "DTK_NonlinearProblemTraits.hpp"

#include <Teuchos_RCP.hpp>

#include <Shards_CellTopology.hpp>

#include <Intrepid_Basis.hpp>
#include <Intrepid_FieldContainer.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \class ProjectPointToFaceNonlinearProblem
 * \brief Nonlinear problem for projecting a point into the reference
 * frame of a face.
 */
//---------------------------------------------------------------------------//
class ProjectPointToFaceNonlinearProblem
{
  public:

    //! Multidimensional array typedefs.
    typedef Intrepid::FieldContainer<double>                       md_array_type;
    typedef typename Intrepid::FieldContainer<double>::scalar_type Scalar;

    // Constructor.
    ProjectPointToFaceNonlinearProblem( 
	const Teuchos::RCP<Intrepid::Basis<Scalar,Intrepid::FieldContainer<double> > >& face_basis,
	const Intrepid::FieldContainer<double>& point,
	const Intrepid::FieldContainer<double>& face_nodes,
	const Intrepid::FieldContainer<double>& face_node_normals );

    //! Destructor.
    ~ProjectPointToFaceNonlinearProblem() { /* ... */ }

    //! Update the state of the problem given the new solution vector.
    void updateState( const Intrepid::FieldContainer<double>& u );

    // Evaluate the nonlinear residual.
    void evaluateResidual( const Intrepid::FieldContainer<double>& u, 
			   Intrepid::FieldContainer<double>& F ) const;

    // Evaluate the jacobian.
    void evaluateJacobian( const Intrepid::FieldContainer<double>& u, 
			   Intrepid::FieldContainer<double>& J ) const;

  public:

    // Basis of the face we are projecting to. 
    Teuchos::RCP<Intrepid::Basis<Scalar,Intrepid::FieldContainer<double> > > d_face_basis;

    // Point coordinates to be projected.
    Intrepid::FieldContainer<double> d_point;

    // Physical coordinates of the face nodes.
    Intrepid::FieldContainer<double> d_face_nodes;

    // Normal vectors of the face nodes.
    Intrepid::FieldContainer<double> d_face_node_normals;

    // Spatial dimension.
    int d_space_dim;
    
    // Face topological dimension.
    int d_topo_dim;
    
    // Face basis cardinality.
    int d_cardinality;

    // Face basis evaluations.
    Intrepid::FieldContainer<double> d_basis_evals;

    // Face basis gradient evalations.
    Intrepid::FieldContainer<double> d_grad_evals;

    // Evaluation points in natural coordinates.
    Intrepid::FieldContainer<double> d_eval_points;

}; // end ProjectToFaceNonlinearProblem

//---------------------------------------------------------------------------//
// ProjectPointToFaceNonlinearProblem traits implementation.
//---------------------------------------------------------------------------//
template<>
class NonlinearProblemTraits<ProjectPointToFaceNonlinearProblem>
{
  public:

    typedef ProjectPointToFaceNonlinearProblem 
    nonlinear_problem_type;

    typedef typename nonlinear_problem_type::md_array_type MDArray;

    typedef typename nonlinear_problem_type::Scalar Scalar;

    static inline void updateState( 
	nonlinear_problem_type& problem, const Intrepid::FieldContainer<double>& u )
    {
	problem.updateState( u );
    }

    static inline void evaluateResidual( const nonlinear_problem_type& problem,
					 const Intrepid::FieldContainer<double>& u, 
					 Intrepid::FieldContainer<double>& F )
    {
	problem.evaluateResidual( u, F );
    }

    static inline void evaluateJacobian( const nonlinear_problem_type& problem,
					 const Intrepid::FieldContainer<double>& u, 
					 Intrepid::FieldContainer<double>& J )
    {
	problem.evaluateJacobian( u, J );
    }
};

//---------------------------------------------------------------------------//
/*!
 * \class PointInFaceVolumeOfInfluenceNonlinearProblem.
 * \brief Nonlinear problem struct for pointInFaceVolumeOfInfluence. This
 * problem projects a point onto the surface defined by a vertex in the face
 * we are checking and the normal vector of that vertex. 
 */
//---------------------------------------------------------------------------//
class PointInFaceVolumeOfInfluenceNonlinearProblem
{
  public:

    //! Multidimensional array typedefs.
    typedef Intrepid::FieldContainer<double>                       md_array_type;
    typedef typename Intrepid::FieldContainer<double>::scalar_type Scalar;

    // Constructor.
    PointInFaceVolumeOfInfluenceNonlinearProblem( 
	const Intrepid::FieldContainer<double>& point,
	const Intrepid::FieldContainer<double>& face_edge_nodes,
	const Intrepid::FieldContainer<double>& face_edge_node_normals,
	const double c );

    //! Destructor.
    ~PointInFaceVolumeOfInfluenceNonlinearProblem() { /* ... */ }

    // Update the state of the problem given the new solution vector.
    void updateState( const Intrepid::FieldContainer<double>& u );

    // Evaluate the nonlinear residual.
    void evaluateResidual( const Intrepid::FieldContainer<double>& u, 
			   Intrepid::FieldContainer<double>& F ) const;

    // Evaluate the jacobian.
    void evaluateJacobian( const Intrepid::FieldContainer<double>& u, 
			   Intrepid::FieldContainer<double>& J ) const;

  public:

    // Point to check for inclusion in the volume of influence.
    Intrepid::FieldContainer<double> d_point;

    // Physical coordinates of the nodes constructing the edge of the face
    // we are currently checking.
    Intrepid::FieldContainer<double> d_face_edge_nodes;

    // Normal vectors of the nodes constructing the edge of the face that
    // we are currently checking.
    Intrepid::FieldContainer<double> d_face_edge_node_normals;

    // Scale factor.
    double d_c;

    // Spatial dimension.
    int d_space_dim;

    // Face topological dimension.
    int d_topo_dim;

    // Pseudo-face cardinality.
    int d_cardinality;

    // Psuedo-face basis evaluations.
    Intrepid::FieldContainer<double> d_basis_evals;

    // Psuedo-face basis gradient evaluations.
    Intrepid::FieldContainer<double> d_grad_evals;

    // Psuedo-face basis projection natural coordinates.
    Intrepid::FieldContainer<double> d_eval_points;

    // Current projection normal direction.
    Intrepid::FieldContainer<double> d_proj_normal;

    // Extra nodes constructing the pseudo-face.
    Intrepid::FieldContainer<double> d_face_normal_nodes;

}; // end PointInFaceVolumeOfInfluenceNonlinearProblem

//---------------------------------------------------------------------------//
// PointInFaceVolumeOfInfluenceNonlinearProblem traits implementation.
//---------------------------------------------------------------------------//
template<>
class NonlinearProblemTraits<PointInFaceVolumeOfInfluenceNonlinearProblem>
{
  public:

    typedef PointInFaceVolumeOfInfluenceNonlinearProblem 
    nonlinear_problem_type;

    typedef typename nonlinear_problem_type::md_array_type MDArray;

    typedef typename nonlinear_problem_type::Scalar Scalar;

    static inline void updateState( 
	nonlinear_problem_type& problem, const Intrepid::FieldContainer<double>& u )
    {
	problem.updateState( u );
    }

    static inline void evaluateResidual( const nonlinear_problem_type& problem,
					 const Intrepid::FieldContainer<double>& u, 
					 Intrepid::FieldContainer<double>& F )
    {
	problem.evaluateResidual( u, F );
    }

    static inline void evaluateJacobian( const nonlinear_problem_type& problem,
					 const Intrepid::FieldContainer<double>& u, 
					 Intrepid::FieldContainer<double>& J )
    {
	problem.evaluateJacobian( u, J );
    }
};

//---------------------------------------------------------------------------//
/*!
 * \class ProjectPointFeatureToEdgeFeatureNonlinearProblem.
 * \brief Nonlinear problem struct for ProjectBlueFeatureToGreenFeature. This
 * problem projects a feature point onto the a feature edge.
 */
//---------------------------------------------------------------------------//
class ProjectPointFeatureToEdgeFeatureNonlinearProblem
{
  public:

    //! Multidimensional array typedefs.
    typedef Intrepid::FieldContainer<double>                       md_array_type;
    typedef typename Intrepid::FieldContainer<double>::scalar_type Scalar;

    // Constructor.
    ProjectPointFeatureToEdgeFeatureNonlinearProblem( 
	const Intrepid::FieldContainer<double>& point,
	const Intrepid::FieldContainer<double>& edge_nodes,
	const Intrepid::FieldContainer<double>& edge_node_normals,
	const Intrepid::FieldContainer<double>& edge_node_binormals );

    //! Destructor.
    ~ProjectPointFeatureToEdgeFeatureNonlinearProblem() { /* ... */ }

    // Update the state of the problem given the new solution vector.
    void updateState( const Intrepid::FieldContainer<double>& u );

    // Evaluate the nonlinear residual.
    void evaluateResidual( const Intrepid::FieldContainer<double>& u, 
			   Intrepid::FieldContainer<double>& F ) const;

    // Evaluate the jacobian.
    void evaluateJacobian( const Intrepid::FieldContainer<double>& u, 
			   Intrepid::FieldContainer<double>& J ) const;

  public:

    // point to project.
    Intrepid::FieldContainer<double> d_point;

    // Coordinates of the edge nodes.
    Intrepid::FieldContainer<double> d_edge_nodes;

    // Normal vectors of the edge nodes.
    Intrepid::FieldContainer<double> d_edge_node_normals;

    // Binormal vectors of the edge nodes.
    Intrepid::FieldContainer<double> d_edge_node_binormals;

    // Spatial dimension.
    int d_space_dim;

}; // end ProjectPointFeatureToEdgeFeatureNonlinearProblem

//---------------------------------------------------------------------------//
// ProjectPointFeatureToEdgeFeatureNonlinearProblem traits implementation.
//---------------------------------------------------------------------------//
template<>
class NonlinearProblemTraits<ProjectPointFeatureToEdgeFeatureNonlinearProblem>
{
  public:

    typedef ProjectPointFeatureToEdgeFeatureNonlinearProblem 
    nonlinear_problem_type;

    typedef typename nonlinear_problem_type::md_array_type MDArray;

    typedef typename nonlinear_problem_type::Scalar Scalar;

    static inline void updateState( 
	nonlinear_problem_type& problem, const Intrepid::FieldContainer<double>& u )
    {
	problem.updateState( u );
    }

    static inline void evaluateResidual( const nonlinear_problem_type& problem,
					 const Intrepid::FieldContainer<double>& u, 
					 Intrepid::FieldContainer<double>& F )
    {
	problem.evaluateResidual( u, F );
    }

    static inline void evaluateJacobian( const nonlinear_problem_type& problem,
					 const Intrepid::FieldContainer<double>& u, 
					 Intrepid::FieldContainer<double>& J )
    {
	problem.evaluateJacobian( u, J );
    }
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_PROJECTIONPRIMITIVENONLINEARPROBLEMS_HPP

//---------------------------------------------------------------------------//
// end DTK_ProjectionPrimitiveNonlinearProblems.hpp
//---------------------------------------------------------------------------//

