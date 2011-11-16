//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   implemenation/Physics_A.hh
 * \author Stuart Slattery
 * \date   Wed Oct 05 09:38:34 2011
 * \brief  Physics_A class definition.
 */
//---------------------------------------------------------------------------//
// $Id: template.hh,v 1.4 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#ifndef implementation_Physics_A_hh
#define implementation_Physics_A_hh

#include <vector>
#include <map>

#include "comm/global.hh"

namespace physics_A
{

//===========================================================================//
/*!
 * \class Physics_A
 * \brief This class is an implementation for unit testing data transfer with
 * the coupler. Physics_A populates a data vector on a uniform grid at the
 * center of the grid cells. 
 */
//===========================================================================//

class Physics_A 
{
  public:
    //@{
    //! Useful typedefs.
    typedef nemesis::Communicator_t               Communicator;
    //@}

  private:

    // MPI Communicator
    Communicator d_comm;

    // Mesh boundaries in x direction.
    std::vector<double> d_x_edges;

    // Mesh boundaries in y direction.
    std::vector<double> d_y_edges;

    // Data vector.
    std::vector<double> d_X;

    // Handle map
    std::map<int,int> handle_map;


  public:

    // Constructor.
    Physics_A(Communicator comm,
	      double x_min, double x_max,
	      double y_min, double y_max,
	      double x_cells, double y_cells);

    // Destructor.
    ~Physics_A();

    // Get the communicator.
    const Communicator& comm();

    // Standalone solve.
    void solve();

    // Given (x,y) coordinates, return true if it exists in the local
    // process rank domain, false if not.
    bool get_xy_info(double x, 
		     double y, 
		     int handle);

    // Given a handle, get that part of the data vector.
    void get_data(const std::vector<int> &handles, 
		  std::vector<double> data);

    //@{ 
    //! Get the mesh edges for testing purposes.
    void get_x_edges( const std::vector<double> &x)
    { x = d_x_edges; }

    void get_y_edges( const std::vector<double> &y)
    { y = d_y_edges; }

    void get_z_edges( const std::vector<double> &z)
    { z = d_z_edges; }
    //@}
};

} // end namespace physics_A

#endif // implementation_Physics_A_hh

//---------------------------------------------------------------------------//
//              end of implementation/Physics_A.hh
//---------------------------------------------------------------------------//
