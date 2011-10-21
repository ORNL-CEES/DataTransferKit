//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   peaks_transfer/TE_Physics_B.hh
 * \author Stuart Slattery
 * \date   Wed Oct 05 10:57:33 2011
 * \brief  Transfer Evaluator for Physics_B.
 */
//---------------------------------------------------------------------------//
// $Id: template.hh,v 1.4 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#ifndef peaks_transfer_TE_Physics_B_hh
#define peaks_transfer_TE_Physics_B_hh

#include "Physics_B.hh"
#include "../../src/core/Transfer_Evaluator.hh"

namespace peaks_transfer
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

    //! Register a field associated with the entities.
    bool register_field(std::string field_name);

    // Register a vector of interleaved point coordinates.
    void register_xyz(std::string field_name,
		      std::vector<double> &points,
		      Vec_Handle &handles);

    // Push data onto a field.
    void push_data(std::string field_name,
		   Handle handle,
		   double data);
};

} // end namespace peaks_transfer

#endif // peaks_transfer_TE_Physics_B_hh

//---------------------------------------------------------------------------//
//              end of peaks_transfer/TE_Physics_B.hh
//---------------------------------------------------------------------------//
