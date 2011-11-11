//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   implementation/Physics_B.hh
 * \author Stuart Slattery
 * \date   Wed Oct 05 09:38:40 2011
 * \brief  Physics_B class definition.
 */
//---------------------------------------------------------------------------//
// $Id: template.hh,v 1.4 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#ifndef implementation_Physics_B_hh
#define implementation_Physics_B_hh

#include <vector>

#include "comm/global.hh"

namespace physics_B
{

//===========================================================================//
/*!
 * \class Physics_B
 * \brief This class is an implementation for unit testing data transfer with
 * the coupler. It acquires its data vector from the data vector of
 * Physics_A. The Physics_B data vector is associated with the mesh vertices.
 */
//===========================================================================//

class Physics_B 
{
  public:
    //@{
    // Useful typedefs.
    typedef std::vector<double>                   Vector_Dbl;
    typedef nemesis::Communicator_t               Communicator;
    //@}

  private:

    // MPI communicator.
    Communicator d_comm;

    // Mesh boundaries in x direction.
    Vector_Dbl d_x_edges;

    // Mesh boundaries in y direction.
    Vector_Dbl d_y_edges;

    // Data vector.
    Vector_Dbl d_X;

  public:

    // Constructor.
    Physics_B(Communicator comm,
	      double x_min, double x_max,
	      double y_min, double y_max,
	      double x_nodes, double y_nodes);

    // Destructor.
    ~Physics_B();

    // Return the communicator.
    const Communicator& comm();

    // Set an element in the source vector.
    void set_data(int handle, double data);

    // Return a const reference data vector.
    const Vector_Dbl& data() { return X; }
};

} // end namespace physics_B

#endif // implementation_Physics_B_hh

//---------------------------------------------------------------------------//
//              end of implementation/Physics_B.hh
//---------------------------------------------------------------------------//
