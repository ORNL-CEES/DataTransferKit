//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   example/Physics_B.hh
 * \author Stuart Slattery
 * \date   Wed Oct 05 09:38:40 2011
 * \brief  Physics_B class definition.
 */
//---------------------------------------------------------------------------//
// $Id: template.hh,v 1.4 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#ifndef example_Physics_B_hh
#define example_Physics_B_hh

#include <vector>

namespace physics_B
{

//===========================================================================//
/*!
 * \class Physics_B
 * \brief Physics_B is part of the simple data transfer example. It acquires
 * its source term from the state vector of Physics_A.
 */
//===========================================================================//

class Physics_B 
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

    // Source vector.
    Vector_Dbl b;

  public:

    // Constructor.
    Physics_B(double x_min, double x_max,
	      double y_min, double y_max,
	      double x_nodes, double y_nodes);

    // Destructor.
    ~Physics_B();

    // Return a const pointer state vector.
    const Vector_Dbl* state() { return &X; }

    // Return a pointer to the source vector.
    Vector_Dbl* source() { return &b; }

    // Return a const reference to the x_edges vector.
    const Vector_Dbl x_domain() { return x_edges; }

    // Return a const reference to the y_edges vector.
    const Vector_Dbl y_domain() { return y_edges; }

    // Set an element in the source vector.
    void set_source(int index, double source);
    
    // Plot the source term.
    void plot_source();
};

} // end namespace physics_B

#endif // example_Physics_B_hh

//---------------------------------------------------------------------------//
//              end of example/Physics_B.hh
//---------------------------------------------------------------------------//
