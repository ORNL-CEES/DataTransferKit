//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   example/Post_Proc.cc
 * \author Stuart Slattery
 * \date   Sun Oct 02 16:02:51 2011
 * \brief  Post processing utilities.
 */
//---------------------------------------------------------------------------//
// $Id: template.cc,v 1.3 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#include "Post_Proc.hh"
#include <vector>
#include <cassert>
#include <iostream>

namespace utils
{
//---------------------------------------------------------------------------//
// Constructor.
Post_Proc::Post_Proc(Edge_Vector x_edges, Edge_Vector y_edges)
{
    // Create the MOAB instance.
    MBI = new moab::Core();

    // Create the vertices.
    std::vector<double> coords;
    int num_verts = x_edges.size()*y_edges.size();
    coords.resize(num_verts*3);

    int ival, jval, idx;
    for (jval = 0; jval < y_edges.size(); ++jval) {
	for (ival = 0; ival < x_edges.size(); ++ival) {
	    idx = ival + x_edges.size()*jval;
	    assert(idx < coords.size());
	    coords[3*idx] = x_edges[ival];
	    coords[3*idx + 1] = y_edges[jval];
	    coords[3*idx + 2] = 0.0;
	}
    }

    rval = MBI->create_vertices(&coords[0], num_verts, vtx_range);
    assert(moab::MB_SUCCESS == rval);

    // Create the quads.
    moab::EntityHandle conn[4];
    int iq, jq;
    for (jq = 0; jq != y_edges.size()-1; jq++) {
	for (iq = 0; iq != x_edges.size()-1; iq++) {
	    idx = iq + (x_edges.size())*jq;
	    assert(idx < vtx_range.size());
	    conn[0] = vtx_range[idx];
	    conn[1] = vtx_range[idx + 1];
	    conn[2] = vtx_range[idx + x_edges.size() + 1];
	    conn[3] = vtx_range[idx + x_edges.size()];

	    moab::EntityHandle this_quad;
	    rval = MBI->create_element(moab::MBQUAD, conn, 4, this_quad);
	    assert(moab::MB_SUCCESS == rval);
	    quad_range.insert(this_quad);
	}
    }
}
 
//---------------------------------------------------------------------------//
// Destructor.
Post_Proc::~Post_Proc()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Tag a field vector onto the mesh quads.
void Post_Proc::add_quad_tag(std::vector<double> field, std::string name)
{
    // Create the tag on the mesh.
    moab::Tag this_tag;
    rval = MBI->tag_get_handle(&name[0], 
			       1,
			       moab::MB_TYPE_DOUBLE,
			       this_tag,
			       moab::MB_TAG_DENSE|moab::MB_TAG_CREAT);
    assert(moab::MB_SUCCESS == rval);

    // Write the field vector data into the tags.
    assert(quad_range.size() == field.size());
    rval = MBI->tag_set_data(this_tag, quad_range, &field[0]);
    assert(moab::MB_SUCCESS == rval);
}

//---------------------------------------------------------------------------//
// Tag a field vector onto the mesh vertices.
void Post_Proc::add_vertex_tag(std::vector<double> field, std::string name)
{
    // Create the tag on the mesh.
    moab::Tag this_tag;
    rval = MBI->tag_get_handle(&name[0], 
			       1,
			       moab::MB_TYPE_DOUBLE,
			       this_tag,
			       moab::MB_TAG_DENSE|moab::MB_TAG_CREAT);
    assert(moab::MB_SUCCESS == rval);

    // Write the field vector data into the tags.
    assert(vtx_range.size() == field.size());
    rval = MBI->tag_set_data(this_tag, vtx_range, &field[0]);
    assert(moab::MB_SUCCESS == rval);
}

//---------------------------------------------------------------------------//
// Write the mesh database to a file.
void Post_Proc::write(std::string filename)
{
    rval = MBI->write_mesh(&filename[0]);
    assert(moab::MB_SUCCESS == rval);
}

//---------------------------------------------------------------------------//

} // end namespace utils

//---------------------------------------------------------------------------//
//                 end of Post_Proc.cc
//---------------------------------------------------------------------------//
