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
    return true;
}

//---------------------------------------------------------------------------//
// Given a (x,y,z) coordinates, return the local process rank in
// which that point exists and the index into the local state vector that
// will be applied at that point. Return true if point is in the local
// domain, false if not.
bool TE_Physics_A::find_xyz(double x, 
			    double y, 
			    double z, 
			    int &rank,
			    Const_Iterator &value_iterator)
{
    domain = a->get_xy_info(x, y, rank, index);
}

//---------------------------------------------------------------------------//
// Pull data from a field.
void TE_Physics_A::pull_data(int index, double &data)
{
    a->get_state(index, data);
}

//---------------------------------------------------------------------------//
// Push data onto a field.
void TE_Physics_A::push_data(int index, double data)
{

}

//---------------------------------------------------------------------------//

} // end namespace dtransfer

//---------------------------------------------------------------------------//
//                 end of TE_Physics_A.cc
//---------------------------------------------------------------------------//
