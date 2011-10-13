//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   example/TE_Physics_B.hh
 * \author Stuart Slattery
 * \date   Wed Oct 05 10:57:33 2011
 * \brief  Transfer Evaluator for Physics_B.
 */
//---------------------------------------------------------------------------//
// $Id: template.hh,v 1.4 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#ifndef example_TE_Physics_B_hh
#define example_TE_Physics_B_hh

#include "Physics_B.hh"
#include "../../src/core/Transfer_Evaluator.hh"

namespace dtransfer
{

//===========================================================================//
/*!
 * \class TE_Physics_B
 * \brief Transfer Evaluator for Physics_B.
 */
//===========================================================================//

class TE_Physics_B : public Transfer_Evaluator
{
  private:

    // Pointer to physics_A object
    physics_B::Physics_B* b;

  public:

    // Constructor.
    TE_Physics_B(physics_B::Physics_B* b_);

    // Destructor
    ~TE_Physics_B();

    // Register a vector of interleaved point coordinates.
    void register_xyz(std::vector<double> &points);

    //! Register a field associated with the entities.
    void register_field(std::string field_name);

    // Given a (x,y,z) coordinates, return the local process rank in
    // which that point exists and the index into the local state vector that
    // will be applied at that point. Return true if point is in the local
    // domain, false if not.
    void find_xyz(double x, 
		  double y, 
		  double z, 
		  int &rank,
		  int &index,
		  bool domain);

    // Pull from a field.
    void pull_data(int index, double &data);

    // Push data onto a field.
    void push_data(int index, double data);
};

} // end namespace dtransfer

#endif // example_TE_Physics_B_hh

//---------------------------------------------------------------------------//
//              end of example/TE_Physics_B.hh
//---------------------------------------------------------------------------//
