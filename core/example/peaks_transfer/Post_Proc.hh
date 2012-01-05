//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   example/Post_Proc.hh
 * \author Stuart Slattery
 * \date   Sun Oct 02 16:02:51 2011
 * \brief  Post-Processing Utilities
 */
//---------------------------------------------------------------------------//
// $Id: template.hh,v 1.4 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#ifndef example_Post_Proc_hh
#define example_Post_Proc_hh

#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include "moab/Range.hpp"
#include <string>
#include <vector>

namespace utils
{

//===========================================================================//
/*!
 * \class Post_Proc
 * \brief Post-Processing Utilities.
 */
//===========================================================================//

class Post_Proc 
{
  public:

    //@{
    //! Useful Typedefs.
    typedef std::vector<double>                 Edge_Vector;

  private:

    // MOAB instance.
    moab::Interface *MBI;

    // MOAB errorcode;
    moab::ErrorCode rval;

    // Vertex range.
    moab::Range vtx_range;

    // Quad range.
    moab::Range quad_range;

 public:

    // Constructor.
    Post_Proc(Edge_Vector x_edges, Edge_Vector y_edges);

    // Destructor.
    ~Post_Proc();

    // Tag a field vector onto the mesh quads.
    void add_quad_tag(std::vector<double> field, std::string name);

    // Tag a field vector onto the mesh vertices.
    void add_vertex_tag(std::vector<double> field, std::string name);    

    // Write the mesh database to a file
    void write(std::string filename);
};

} // end namespace utils

#endif // example_Post_Proc_hh

//---------------------------------------------------------------------------//
//              end of example/Post_Proc.hh
//---------------------------------------------------------------------------//
