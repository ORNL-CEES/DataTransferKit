//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   example/TE_Physics_A.hh
 * \author Stuart Slattery
 * \date   Wed Oct 05 10:57:30 2011
 * \brief  Transfer Evaluator for Physics_A.
 */
//---------------------------------------------------------------------------//
// $Id: template.hh,v 1.4 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#ifndef example_TE_Physics_A_hh
#define example_TE_Physics_A_hh

#include "Physics_A.hh"
#include "../../src/core/Transfer_Evaluator.hh"

namespace dtransfer
{

//===========================================================================//
/*!
 * \class TE_Physics_A
 * \brief Transfer Evaluator for Physics_A. Derives from Transfer_Evaluator.
 */
//===========================================================================//

class TE_Physics_A : public Transfer_Evaluator
{
  private:

    // Pointer to physics_A object
    physics_A::Physics_A* a;

  public:

    // Constructor.
    TE_Physics_A(physics_A::Physics_A* a_);

    // Destructor
    ~TE_Physics_A();

    // Register a vector of interleaved point coordinates.
    void register_xyz(std::vector<double> &points);

    //! Register a field associated with the entities.
    bool register_field(std::string field_name);

    // Given a (x,y,z) coordinates, return the local process rank in
    // which that point exists and the index into the local state vector that
    // will be applied at that point. Return true if point is in the local
    // domain, false if not.
    bool find_xyz(double x, 
		  double y, 
		  double z, 
		  int &rank,
                  Const_Iterator &value_iterator)
};

} // end namespace dtransfer

#endif // example_TE_Physics_A_hh

//---------------------------------------------------------------------------//
//              end of example/TE_Physics_A.hh
//---------------------------------------------------------------------------//
