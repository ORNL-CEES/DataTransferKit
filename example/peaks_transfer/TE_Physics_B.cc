//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   example/TE_Physics_B.cc
 * \author Stuart Slattery
 * \date   Wed Oct 05 10:57:33 2011
 * \brief  Transfer Evaluator for Physics_B.
 */
//---------------------------------------------------------------------------//
// $Id: template.cc,v 1.3 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#include "TE_Physics_B.hh"

namespace dtransfer
{

//---------------------------------------------------------------------------//
// Constructor.
TE_Physics_B::TE_Physics_B(physics_B::Physics_B* b_)
    : b(b_)
{ /* ... */ }

//---------------------------------------------------------------------------//
// Destructor.
TE_Physics_B::~TE_Physics_B()
{ /* ... */ }

//---------------------------------------------------------------------------//
//! Register a field associated with the entities.
bool TE_Physics_B::register_field(std::string field_name)
{
    if (field_name == "PEAKS") return true;
    else return false;
}

//---------------------------------------------------------------------------//
// Register a vector of interleaved point coordinates.
void TE_Physics_B::register_xyz(std::string field_name,
				std::vector<double> &points,
				std::vector<Handle> &handles)
{
    if (field_name == "PEAKS")
    {
	// Coordinate vector.
	points.resize( 3*(b->x_domain().size())*(b->y_domain().size()) );

	// Set iterators.
	physics_B::Physics_B::Vector_Dbl::const_iterator x_it;
	physics_B::Physics_B::Vector_Dbl::const_iterator y_it;
	std::vector<double>::iterator coord_it = points.begin();
	Vec_Handle::iterator handle_it = handles.begin();
	int handle_counter = 0;
	// Populate the coordinate vector.
	for (y_it = b->y_domain().begin(); y_it != b->y_domain().end(); ++y_it)
	{
	    for (x_it = b->x_domain().begin(); x_it != b->x_domain().end(); ++x_it)
	    {
		*coord_it = *x_it;
		++coord_it;

		*coord_it = *y_it;
		++coord_it;

		*coord_it = 0.0;
		++coord_it;

		*handle_it = handle_counter;
		++handle_it;
		++handle_counter;
	    }
	}
    }
}

//---------------------------------------------------------------------------//
// Push data onto a field.
void TE_Physics_B::push_data(std::string field_name,
			     Handle handle, 
			     double data)
{
    b->set_source(handle, data);
}

//---------------------------------------------------------------------------//

} // end namespace dtransfer

//---------------------------------------------------------------------------//
//                 end of TE_Physics_B.cc
//---------------------------------------------------------------------------//
