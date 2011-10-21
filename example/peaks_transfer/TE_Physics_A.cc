//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   example/TE_Physics_A.cc
 * \author Stuart Slattery 
 * \date   Wed Oct 05 10:57:30 2011
 * \brief  Transfer Evaluator for Physics_A.
 */
//---------------------------------------------------------------------------//
// $Id: template.cc,v 1.3 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#include "TE_Physics_A.hh"

namespace dtransfer
{
//---------------------------------------------------------------------------//
// Constructor.
TE_Physics_A::TE_Physics_A(physics_A::Physics_A* a_)
    : a(a_)
{ /* ... */ }

//---------------------------------------------------------------------------//
// Destructor.
TE_Physics_A::~TE_Physics_A()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Register a vector of interleaved point coordinates.
void TE_Physics_A::register_xyz(std::vector<double> &points)
{
}

//---------------------------------------------------------------------------//
//! Register a field associated with the entities.
bool TE_Physics_A::register_field(std::string field_name)
{
    if (field_name == "PEAKS") return true;
    else return false;
}

//---------------------------------------------------------------------------//
// Given a (x,y,z) coordinates, return the local process rank in
// which that point exists and the index into the local state vector that
// will be applied at that point. Return true if point is in the local
// domain, false if not.
bool TE_Physics_A::find_xyz(double x, 
			    double y, 
			    double z,
			    Handle &handle)
{
    return a->get_xy_info(x, y, handle);
}

//---------------------------------------------------------------------------//
// Given an entity handle, get the field data associated with that handle.
void TE_Physics_A::pull_data(std::string field_name,
			     Handle handle,
			     double &data)
{
    if (field_name == "PEAKS")
    {
	a->get_state(handle, data);
    }
}

//---------------------------------------------------------------------------//

} // end namespace dtransfer

//---------------------------------------------------------------------------//
//                 end of TE_Physics_A.cc
//---------------------------------------------------------------------------//
