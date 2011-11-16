//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   implemenation/TE_Physics_B.hh
 * \author Stuart Slattery
 * \date   Wed Oct 05 10:57:33 2011
 * \brief  Transfer Evaluator for Physics_B.
 */
//---------------------------------------------------------------------------//
// $Id: template.hh,v 1.4 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#ifndef implementation_TE_Physics_B_hh
#define implementation_TE_Physics_B_hh

#include <vector>

#include "Physics_B.hh"
#include "core/Transfer_Evaluator.hh"
#include "comm/global.hh"

namespace coupler
{

//===========================================================================//
/*!
 * \class TE_Physics_B
 * \brief Transfer Evaluator for Physics_B.
 */
//===========================================================================//

class TE_Physics_B : public Transfer_Evaluator
{

  public:

    typedef nemesis::Communicator_t          Communicator;

  private:

    // Pointer to physics_A object
    physics_B::Physics_B* b;

    // Point coordinate vector.
    std::vector<CoordinateType> points;

    // Handle vector.
    std::vector<HandleType> handles;

  public:

    // Constructor.
    TE_Physics_B(physics_B::Physics_B* b_);

    // Destructor
    ~TE_Physics_B();

    //! Register the communicator.
    void register_comm(const Communicator &comm);

    //! Register a field associated with the entities.
    bool register_field(const std::string &field_name);

    // Register a vector of interleaved point coordinates.
    void register_points(const std::string &field_name,
			 std::vector<HandleType> &handles,
			 std::vector<CoordinateType> &coordinates);

    // Push data onto a field.
    void push_data(const std::string &field_name,
		   const std::vector<HandleType> &handles,
		   const std::vector<DataType> &data);
};

} // end namespace coupler

#endif // implmentation_TE_Physics_B_hh

//---------------------------------------------------------------------------//
//              end of implemenation/TE_Physics_B.hh
//---------------------------------------------------------------------------//
