//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   example/Physics_B.cc
 * \author Stuart Slattery
 * \date   Wed Oct 05 09:38:40 2011
 * \brief  Physics_B member definitions.
 */
//---------------------------------------------------------------------------//
// $Id: template.cc,v 1.3 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#include "Physics_B.hh"
#include "Post_Proc.hh"
#include <cassert>

namespace physics_B
{
//---------------------------------------------------------------------------//
// Constructor.
Physics_B::Physics_B(double x_min, double x_max,
		     double y_min, double y_max,
		     double x_nodes, double y_nodes)
{
    // Make sure max > min.
    assert(x_max > x_min);
    assert(y_max > y_min);

    // Generate the mesh (regular cartesian grid).
    x_edges.resize(x_nodes);
    y_edges.resize(y_nodes);

    // x_edges
    double x_width = (x_max - x_min) / (x_nodes-1);
    for (int i = 0; i < x_nodes; ++i)
    {
	x_edges[i] = x_min + i*x_width;
    }

    // y_edges
    double y_width = (y_max - y_min) / (y_nodes-1);
    for (int j = 0; j < y_nodes; ++j)
    {
	y_edges[j] = y_min + j*y_width;
    }

    // Resize the source vector to correspond to nodal values and fill
    // with 0.
    b.resize(x_nodes*y_nodes);
    Vector_Dbl::iterator it;
    for (it = b.begin(); it != b.end(); ++it)
    {
	*it = 0.0;
    }
}

//---------------------------------------------------------------------------//
// Destructor.
Physics_B::~Physics_B()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Set an element in the source vector.
void Physics_B::set_source(int index, double source)
{
    assert( index < b.size() );
    b[index] = source;
}

//---------------------------------------------------------------------------//
// Plot the source term.
void Physics_B::plot_source()
{
    // Make a new post-processing object.
    utils::Post_Proc post_proc(x_edges, y_edges);

    // Tag the state vector onto the mesh vertices.
    post_proc.add_vertex_tag(b, "b");

    // Write the database to file.
    post_proc.write("Physics_B.vtk");
}

} // end namespace physics_B

//---------------------------------------------------------------------------//
//                 end of Physics_B.cc
//---------------------------------------------------------------------------//
