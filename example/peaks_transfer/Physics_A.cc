//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   peaks_example/Physics_A.cc
 * \author Stuart Slattery
 * \date   Wed Oct 05 09:38:34 2011
 * \brief  Physics_A member definitions.
 */
//---------------------------------------------------------------------------//
// $Id: template.cc,v 1.3 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#include "Physics_A.hh"
#include "Post_Proc.hh"
#include <cassert>
#include <cmath>
#include <algorithm>
#include <iostream>

namespace physics_A
{
//---------------------------------------------------------------------------//
// Constructor.
Physics_A::Physics_A(double x_min, double x_max,
		     double y_min, double y_max,
		     double x_cells, double y_cells)
{
    // Make sure max > min.
    assert(x_max > x_min);
    assert(y_max > y_min);

    // Generate the mesh (regular cartesian grid).
    x_edges.resize(x_cells+1);
    y_edges.resize(y_cells+1);

    // x_edges
    double x_width = (x_max - x_min) / x_cells;
    for (int i = 0; i < x_cells+1; ++i)
    {
	x_edges[i] = x_min + i*x_width;
    }

    // y_edges
    double y_width = (y_max - y_min) / y_cells;
    for (int j = 0; j < y_cells+1; ++j)
    {
	y_edges[j] = y_min + j*y_width;
    }

    // Resize the state vector to correspond to cell-centered values and fill
    // with 0.
    X.resize(x_cells*y_cells);
    Vector_Dbl::iterator it;
    for (it = X.begin(); it != X.end(); ++it)
    {
	*it = 0.0;
    }
}

//---------------------------------------------------------------------------//
// Destructor.
Physics_A::~Physics_A()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Do standalone solve.
void Physics_A::solve()
{
    // Populate the state vector with the peaks function. Here, the function
    // is applied to cell centered values.
    int index;
    double x_center;
    double y_center;
    for (int j = 0; j < (y_edges.size() - 1); ++j)
    {
	for (int i = 0; i < (x_edges.size() - 1); ++i)
	{
	    index = i + j*(x_edges.size() - 1);
	    assert(index < X.size());
	    
	    x_center = (x_edges[i+1] + x_edges[i]) / 2;
	    y_center = (y_edges[j+1] + y_edges[j]) / 2;

	    X[index] = 
		(
		       3
		       *(1 - x_center)*(1 - x_center)
		       *exp(-(x_center*x_center) - (y_center+1)*(y_center+1))
		)
		       -
		(
		    10
		    *(x_center/5 
		      - x_center*x_center*x_center
		      - y_center*y_center*y_center*y_center*y_center)
		    *exp(-(x_center*x_center) - y_center*y_center)
		)
		    -
		(
		    (1/3)
		    *exp(-(x_center+1)*(x_center+1) - y_center*y_center)
		);
	}
    }
}

//---------------------------------------------------------------------------//
// Given a (x,y) coordinates, return true if in the local domain, false if
// not. If in the local domain, populate a handle argument with an indicator
// to the cell that we found the point in.
bool Physics_A::get_xy_info(double x, 
			    double y, 
			    int handle)
{
    // Search the x domain.
    int i;
    i = (std::lower_bound(x_edges.begin(), x_edges.end(), x))
	- x_edges.begin();
    if (i > x_edges.size() - 2) i = x_edges.size() - 2;

    // Search the y domain.
    int j;
    j = (std::lower_bound(y_edges.begin(), y_edges.end(), y))
	- y_edges.begin();
    if (j > y_edges.size() - 2) j = y_edges.size() - 2;

    // Get the index into the state vector and associate this with the given
    // handle. 
    int index = i + j*(x_edges.size() - 1);
    handle_map.insert( std::pair<int,int>(handle, index);

    return true;
}

//---------------------------------------------------------------------------//
// Given a handle, get that part of the state vector.
void Physics_A::get_state(int handle, double data)
{
    data = X[ handle_map[handle] ];
}

//---------------------------------------------------------------------------//
// Plot the state vector.
void Physics_A::plot_state()
{
    // Make a new post-processing object.
    utils::Post_Proc post_proc(x_edges, y_edges);

    // Tag the state vector onto the mesh elements.
    post_proc.add_quad_tag(X, "X");

    // Write the database to file.
    post_proc.write("Physics_A.vtk");
}

//---------------------------------------------------------------------------//

} // end namespace physics_A

//---------------------------------------------------------------------------//
//                 end of Physics_A.c
//---------------------------------------------------------------------------//
