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
Physics_B::Physics_B(Communicator comm,
		     double x_min, double x_max,
		     double y_min, double y_max,
		     double x_nodes, double y_nodes)
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
	    x_max_local = x_max;
	    y_min_local = y_min;
	    y_max_local = y_min + (y_edges.size() - 1)*y_width;

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

    Nemesis::reset_internal_comm();
}

//---------------------------------------------------------------------------//
// Destructor.
Physics_B::~Physics_B()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Get the communicator.
const Communicator& Physics_B::comm()
{
    return d_comm;
}

//---------------------------------------------------------------------------//
// Set an element in the source vector.
void Physics_B::set_data(const std::vector<int> &handles, 
			 const std::vector<double> &data)
{
    std::vector<int>::const_iterator handle_it;
    std::vector<double>::const_iterator data_it = data.begin();
    for (handle_it = handles.begin(); handle_it != handles.end(); ++handle_it)
    {
	b[*handle_it] = *data_it;
	++data_it;
    }
}

} // end namespace physics_B

//---------------------------------------------------------------------------//
//                 end of Physics_B.cc
//---------------------------------------------------------------------------//
