//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   implementation/TE_Physics_A.hh
 * \author Stuart Slattery
 * \date   Wed Oct 05 10:57:30 2011
 * \brief  Transfer Evaluator for Physics_A.
 */
//---------------------------------------------------------------------------//
// $Id: template.hh,v 1.4 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#ifndef implementation_TE_Physics_A_hh
#define implementation_TE_Physics_A_hh

#include "Physics_A.hh"
#include "core/Transfer_Evaluator.hh"
#include "comm/global.hh"

namespace coupler
{

//===========================================================================//
/*!
 * \class TE_Physics_A
 * \brief Transfer Evaluator for Physics_A. Derives from Transfer_Evaluator.
 */
//===========================================================================//

class TE_Physics_A : public Transfer_Evaluator
{
  public:

    typedef nemesis::Communicator_t            Communicator;

  private:

    // Pointer to physics_A object
    physics_A::Physics_A* a;

  public:

    // Constructor.
    TE_Physics_A(physics_A::Physics_A* a_);

    // Destructor
    ~TE_Physics_A();

    //! Register the communicator.
    void register_comm(Communicator &comm);

    //! Register a field associated with the entities.
    bool register_field(std::string field_name);

    // Given a (x,y,z) coordinates, return true if point is in the local
    // domain, false if not. Provide a handle to the entity it was located in.
    bool find_xyz(double x, 
		  double y, 
		  double z,
		  Handle &handle);

    // Given an entity handle, get the field data associated with that
    // handle. 
    void pull_data(std::string field_name,
		   Handle handle,
		   DataType &data);
};

} // end namespace coupler

#endif // implementation_TE_Physics_A_hh

//---------------------------------------------------------------------------//
//              end of implementation/TE_Physics_A.hh
//---------------------------------------------------------------------------//
