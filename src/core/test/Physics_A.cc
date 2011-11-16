//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   implementation/Physics_A.cc
 * \author Stuart Slattery
 * \date   Wed Oct 05 09:38:34 2011
 * \brief  Physics_A member definitions.
 */
//---------------------------------------------------------------------------//
// $Id: template.cc,v 1.3 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#include <cassert>
#include <cmath>
#include <algorithm>
#include <iostream>

#include "Physics_A.hh"

namespace physics_A
{
//---------------------------------------------------------------------------//
// Constructor.
Physics_A::Physics_A(Communicator comm,
		     double x_min, double x_max,
		     double y_min, double y_max,
		     int x_cells, int y_cells)
    : d_comm(comm)
{
    nemesis::set_internal_comm(d_comm);

    // Make sure max > min.
    assert(x_max > x_min);
    assert(y_max > y_min);

    // Generate the mesh (regular cartesian grid).
    double x_width = (x_max - x_min) / x_cells;
    double y_width = (y_max - y_min) / y_cells;

    // Local mesh bounds.
    double x_min_local;
    double x_max_local;
    double y_min_local;
    double y_max_local;

    // 1 process partitioning
    if ( nemesis::nodes() == 1 )
    {
	x_edges.resize(x_cells+1);
	y_edges.resize(y_cells+1);

	// x_edges
	for (int i = 0; i < x_cells+1; ++i)
	{
	    x_edges[i] = x_min + i*x_width;
	}

	// y_edges
	for (int j = 0; j < y_cells+1; ++j)
	{
	    y_edges[j] = y_min + j*y_width;
	}

	// Resize the data vector to correspond to cell-centered values and
	// fill with 0.
	X.resize(x_cells*y_cells);
	Vector_Dbl::iterator it;
	for (it = X.begin(); it != X.end(); ++it)
	{
	    *it = 0.0;
	}
    }

    // 2 process partitioning
    else if ( nemesis::nodes() == 2 )
    {

	// process 0
	if ( nemesis::node() == 0 )
	{
	    x_edges.resize( (int) (x_cells+1) / 2);
	    y_edges.resize(y_cells+1);

	    x_min_local = x_min;
	    x_max_local = x_min_local + (x_edges.size() - 1)*x_width;
	    y_min_local = y_min;
	    y_max_local = y_max;

	    // x_edges
	    for (int i = 0; i < x_cells+1; ++i)
	    {
		x_edges[i] = x_min_local + i*x_width;
	    }

	    // y_edges
	    for (int j = 0; j < y_cells+1; ++j)
	    {
		y_edges[j] = y_min_local + j*y_width;
	    }

	    // Resize the data vector to correspond to cell-centered values
	    // and fill witth 0.
	    X.resize(x_cells*y_cells);
	    Vector_Dbl::iterator it;
	    for (it = X.begin(); it != X.end(); ++it)
	    {
		*it = 0.0;
	    }
	}

	// process 1
	if ( nemesis::node() == 1 )
	{
	    x_edges.resize( (int) x_cells+1 - (x_cells+1) / 2 );
	    y_edges.resize(y_cells+1);

	    x_min_local = x_min_local + (x_edges.size() - 1)*x_width;
	    x_max_local = x_max;
	    y_min_local = y_min;
	    y_max_local = y_max;

	    // x_edges
	    for (int i = 0; i < x_cells+1; ++i)
	    {
		x_edges[i] = x_min_local + i*x_width;
	    }

	    // y_edges
	    for (int j = 0; j < y_cells+1; ++j)
	    {
		y_edges[j] = y_min_local + j*y_width;
	    }

	    // Resize the data vector to correspond to cell-centered values
	    // and fill with 0.
	    X.resize(x_cells*y_cells);
	    Vector_Dbl::iterator it;
	    for (it = X.begin(); it != X.end(); ++it)
	    {
		*it = 0.0;
	    }
	}
    }

    // 3 process partitioning
    else if ( nemesis::nodes() == 3 )
    {
	// process 0
	if ( nemesis::node() == 0 )
	{
	    x_edges.resize( (int) (x_cells+1) / 2);
	    y_edges.resize( (int) (y_cells+1) / 2);

	    x_min_local = x_min;
	    x_max_local = x_min_local + (x_edges.size() - 1)*x_width;
	    y_min_local = y_min;
	    y_max_local = y_min_local + (y_edges.size() - 1)*y_width;

	    // x_edges
	    for (int i = 0; i < x_cells+1; ++i)
	    {
		x_edges[i] = x_min_local + i*x_width;
	    }

	    // y_edges
	    for (int j = 0; j < y_cells+1; ++j)
	    {
		y_edges[j] = y_min_local + j*y_width;
	    }

	    // Resize the data vector to correspond to cell-centered values
	    // and fill with 0.
	    X.resize(x_cells*y_cells);
	    Vector_Dbl::iterator it;
	    for (it = X.begin(); it != X.end(); ++it)
	    {
		*it = 0.0;
	    }
	}

	// process 1
	if ( nemesis::node() == 1 )
	{
	    x_edges.resize( (int) x_cells+1 - (x_cells+1) / 2 );
	    y_edges.resize( (int) (y_cells+1) / 2 );

	    x_min_local = x_min_local + (x_edges.size() - 1)*x_width;
	    x_max_local = x_max;
	    y_min_local = y_min;
	    y_max_local = y_min_local + (y_edges.size() - 1)*y_width;

	    // x_edges
	    for (int i = 0; i < x_cells+1; ++i)
	    {
		x_edges[i] = x_min_local + i*x_width;
	    }

	    // y_edges
	    for (int j = 0; j < y_cells+1; ++j)
	    {
		y_edges[j] = y_min_local + j*y_width;
	    }

	    // Resize the data vector to correspond to cell-centered values
	    // and fill with 0.
	    X.resize(x_cells*y_cells);
	    Vector_Dbl::iterator it;
	    for (it = X.begin(); it != X.end(); ++it)
	    {
		*it = 0.0;
	    }
	}

	// process 2
	if ( nemesis::node() == 2 )
	{
	    x_edges.resize( x_cells+1 );
	    y_edges.resize( (int) y_cells+1 - (y_cells+1) / 2 );

	    x_min_local = x_min;
	    x_max_local = x_max;
	    y_min_local = y_min + (y_edges.size() - 1)*y_width;
	    y_max_local = y_max;

	    // x_edges
	    for (int i = 0; i < x_cells+1; ++i)
	    {
		x_edges[i] = x_min_local + i*x_width;
	    }

	    // y_edges
	    for (int j = 0; j < y_cells+1; ++j)
	    {
		y_edges[j] = y_min_local + j*y_width;
	    }

	    // Resize the data vector to correspond to cell-centered values
	    // and fill with 0.
	    X.resize(x_cells*y_cells);
	    Vector_Dbl::iterator it;
	    for (it = X.begin(); it != X.end(); ++it)
	    {
		*it = 0.0;
	    }
	}
    }

    // 4 process partitioning
    else if ( nemesis::nodes() == 4 )
    {
	// process 0
	if ( nemesis::node() == 0 )
	{
	    x_edges.resize( x_cells+1 );
	    y_edges.resize( (int) y_cells+1 - (y_cells+1) / 2);

	    x_min_local = x_min;
	    x_max_local = x_min_local + (x_edges.size() - 1)*x_width;
	    y_min_local = y_min;
	    y_max_local = y_min_local + (y_edges.size() - 1)*y_width;

	    // x_edges
	    for (int i = 0; i < x_cells+1; ++i)
	    {
		x_edges[i] = x_min_local + i*x_width;
	    }

	    // y_edges
	    for (int j = 0; j < y_cells+1; ++j)
	    {
		y_edges[j] = y_min_local + j*y_width;
	    }

	    // Resize the data vector to correspond to cell-centered values
	    // and fill with 0.
	    X.resize(x_cells*y_cells);
	    Vector_Dbl::iterator it;
	    for (it = X.begin(); it != X.end(); ++it)
	    {
		*it = 0.0;
	    }
	}

	// process 1
	if ( nemesis::node() == 1 )
	{
	    x_edges.resize( (int) x_cells - (x_cells+1) / 2 );
	    y_edges.resize( (int) y_cells - (y_cells+1) / 2 );

	    x_min_local = x_min_local + (x_edges.size() - 1)*x_width;
	    x_max_local = x_max;
	    y_min_local = y_min;
	    y_max_local = y_min_local + (y_edges.size() - 1)*y_width;

	    // x_edges
	    for (int i = 0; i < x_cells+1; ++i)
	    {
		x_edges[i] = x_min_local + i*x_width;
	    }

	    // y_edges
	    for (int j = 0; j < y_cells+1; ++j)
	    {
		y_edges[j] = y_min_local + j*y_width;
	    }

	    // Resize the data vector to correspond to cell-centered values
	    // and fill with 0.
	    X.resize(x_cells*y_cells);
	    Vector_Dbl::iterator it;
	    for (it = X.begin(); it != X.end(); ++it)
	    {
		*it = 0.0;
	    }
	}

	// process 2
	if ( nemesis::node() == 2 )
	{
	    x_edges.resize( x_cells+1 );
	    y_edges.resize( (int) y_cells - (y_cells+1) / 2 );

	    x_min_local = x_min;
	    x_max_local = x_min_local + (x_edges.size() - 1)*x_width;
	    y_min_local = y_min + (y_edges.size() - 1)*y_width;
	    y_max_local = y_max;

	    // x_edges
	    for (int i = 0; i < x_cells+1; ++i)
	    {
		x_edges[i] = x_min_local + i*x_width;
	    }

	    // y_edges
	    for (int j = 0; j < y_cells+1; ++j)
	    {
		y_edges[j] = y_min_local + j*y_width;
	    }

	    // Resize the data vector to correspond to cell-centered values
	    // and fill with 0.
	    X.resize(x_cells*y_cells);
	    Vector_Dbl::iterator it;
	    for (it = X.begin(); it != X.end(); ++it)
	    {
		*it = 0.0;
	    }
	}
	
	// process 3
	if (nemesis::node() == 3)
	{
	    x_edges.resize( x_cells+1 );
	    y_edges.resize( (int) y_cells - (y_cells+1) / 2 );

	    x_min_local = x_min_local + (x_edges.size() - 1)*x_width;
	    x_max_local = x_max;
	    y_min_local = y_min + (y_edges.size() - 1)*y_width;
	    y_max_local = y_max;

	    // x_edges
	    for (int i = 0; i < x_cells+1; ++i)
	    {
		x_edges[i] = x_min_local + i*x_width;
	    }

	    // y_edges
	    for (int j = 0; j < y_cells+1; ++j)
	    {
		y_edges[j] = y_min_local + j*y_width;
	    }

	    // Resize the data vector to correspond to cell-centered values
	    // and fill with 0.
	    X.resize(x_cells*y_cells);
	    Vector_Dbl::iterator it;
	    for (it = X.begin(); it != X.end(); ++it)
	    {
		*it = 0.0;
	    }
	}
    }

    nemesis::reset_internal_comm();
}

//---------------------------------------------------------------------------//
// Destructor.
Physics_A::~Physics_A()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Get the communicator.
const Communicator& Physics_A::comm()
{
    return d_comm;
}

//---------------------------------------------------------------------------//
// Do standalone solve.
void Physics_A::solve()
{
    // Populate the data vector with the its process id. Here, the data
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

	    X[index] = nemesis::node();
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

    // Get the index into the data vector and associate this with the given
    // handle. 
    int index = i + j*(x_edges.size() - 1);
    handle_map.insert( std::pair<int,int>(handle, index) );

    return true;
}

//---------------------------------------------------------------------------//
// Given a handle, get that part of the data vector.
void Physics_A::get_data(const std::vector<int> &handles,
			 std::vector<double> &data)
{
    std::vector<double>::const_iterator handle_it;
    for (handle_it = handles.begin(); handle_it != handles.end(); ++handle_it) 
    {
	data.push_back( X[ handle_map[handle] ] );
    }
}

//---------------------------------------------------------------------------//

} // end namespace physics_A

//---------------------------------------------------------------------------//
//                 end of Physics_A.c
//---------------------------------------------------------------------------//
