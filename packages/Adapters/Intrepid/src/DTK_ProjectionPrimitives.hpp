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
 * \file   DTK_ProjectionPrimitives.hpp
 * \author Stuart Slattery
 * \brief  A stateless class of projection primitive operations.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_PROJECTIONPRIMITIVES_HPP
#define DTK_PROJECTIONPRIMITIVES_HPP

#include <cmath>

#include "DTK_DBC.hpp"

#include <Teuchos_ParameterList.hpp>

#include <Shards_CellTopology.hpp>

#include <Intrepid_Basis.hpp>
#include <Intrepid_FieldContainer.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \class ProjectionPrimitives
 * \brief A stateless class of projection primitive operations.
 */
//---------------------------------------------------------------------------//
class ProjectionPrimitives
{
  public:

    // Constructor.
    ProjectionPrimitives() { /* ... */ }

    // Destructor.
    ~ProjectionPrimitives() { /* ... */ }

    // Get the center of the reference cell of the given topology.
    static void referenceCellCenter( const shards::CellTopology& cell_topo,
				     Intrepid::FieldContainer<double>& center );

    // Determine if a point is within the volume of influence of a face.
    static bool pointInFaceVolumeOfInfluence(
	const Teuchos::ParameterList& parameters,
	const Intrepid::FieldContainer<double>& point,
	const Intrepid::FieldContainer<double>& face_nodes,
	const Intrepid::FieldContainer<double>& face_node_normals,
	const shards::CellTopology& face_topology );

    // Project a point onto a face and return the physical and parametric
    // coordinates of the projected point on that face. This
    // requires the solution of a nonlinear parameterized projection problem.
    static void projectPointToFace(
	const Teuchos::ParameterList& parameters,
	const Intrepid::FieldContainer<double>& point,
	const Intrepid::FieldContainer<double>& face_nodes,
	const Intrepid::FieldContainer<double>& face_node_normals,
	const shards::CellTopology& face_topology,
	Intrepid::FieldContainer<double>& parametric_point,
	Intrepid::FieldContainer<double>& physical_point,
	int& face_edge_id,
	int& face_node_id );

    // Project a feature point to a feature edge.
    static bool projectPointFeatureToEdgeFeature(
	const Teuchos::ParameterList& parameters,
	const Intrepid::FieldContainer<double>& point,
	const Intrepid::FieldContainer<double>& point_normal,
	const Intrepid::FieldContainer<double>& edge_nodes,
	const Intrepid::FieldContainer<double>& edge_node_normals,
	Intrepid::FieldContainer<double>& projected_point,
	int& edge_node_id );

    // Intersect two edges in 3 dimensions and return their intersection point
    // realized on both edges.
    static bool edgeEdgeIntersection( 
	const Teuchos::ParameterList& parameters,
	const Intrepid::FieldContainer<double>& edge_1,
	const Intrepid::FieldContainer<double>& edge_2,
	const Intrepid::FieldContainer<double>& edge_2_node_normals,
	Intrepid::FieldContainer<double>& edge_1_intersection,
	Intrepid::FieldContainer<double>& edge_2_intersection,
	int& edge_1_node_id,
	int& edge_2_node_id );

  private:

    // Compute the distance of a projected point onto a bilinear surface
    // formed by a face edge and its normals.
    static double distanceToFaceBilinearSurface(
	const Teuchos::ParameterList& parameters,
	const Intrepid::FieldContainer<double>& point,
	const Intrepid::FieldContainer<double>& face_edge_nodes,
	const Intrepid::FieldContainer<double>& face_edge_node_normals );
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_PROJECTIONPRIMITIVES_HPP

//---------------------------------------------------------------------------//
// end DTK_ProjectionPrimitives.hpp
//---------------------------------------------------------------------------//

