//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   example/Physics_A.hh
 * \author Stuart Slattery
 * \date   Wed Oct 05 09:38:34 2011
 * \brief  Physics_A class definition.
 */
//---------------------------------------------------------------------------//
// $Id: template.hh,v 1.4 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#ifndef example_Physics_A_hh
#define example_Physics_A_hh

#include <vector>

namespace physics_A
{

//===========================================================================//
/*!
 * \class Physics_A
 * \brief This class is part of a super simple example for data
 * transfer. Physics_A populates a solution vector with a function on a
 * uniform grid. 
 */
//===========================================================================//

class Physics_A 
{
  public:
    //@{
    // Useful typedefs.
    typedef std::vector<double>                   Vector_Dbl;
    //@}

  private:

    // Mesh boundaries in x direction.
    Vector_Dbl x_edges;

    // Mesh boundaries in y direction.
    Vector_Dbl y_edges;

    // State vector.
    Vector_Dbl X;

    // Source Vector.
    Vector_Dbl b;

  public:

    // Constructor.
    Physics_A(double x_min, double x_max,
	      double y_min, double y_max,
	      double x_cells, double y_cells);

    // Destructor.
    ~Physics_A();

    // Standalone solve.
    void solve();

    // Return a const reference to the x_edges vector.
    const Vector_Dbl x_domain() { return x_edges; }

    // Return a const reference to the y_edges vector.
    const Vector_Dbl y_domain() { return y_edges; }

    // Given a (x,y) coordinates, return the local process rank in
    // which that point exists and the index into the local state vector that
    // will be applied at that point. Return true if point is in the local
    // domain, false if not.
    bool get_xy_info(double x, 
		     double y, 
		     int handle);

    // Given a handle, get that part of the state vector.
    void get_state(int handle, double data);

    // Plot the state vector.
    void plot_state();
};

} // end namespace physics_A

#endif // example_Physics_A_hh

//---------------------------------------------------------------------------//
//              end of example/Physics_A.hh
//---------------------------------------------------------------------------//
