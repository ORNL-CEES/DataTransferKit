//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   peaks_example/TE_Physics_B.cc
 * \author Stuart Slattery
 * \date   Wed Oct 05 10:57:33 2011
 * \brief  Transfer Evaluator for Physics_B.
 */
//---------------------------------------------------------------------------//
// $Id: template.cc,v 1.3 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#include "TE_Physics_B.hh"

namespace peaks_example
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
// Register cartesian coordinates with a field. The coordinate vector
// should be interleaved. The handle vector should consist of globally
// unique handles. These iterators imply contiguous memory storage.
void TE_Physics_B::register_xyz(std::string field_name,
				Coord_Iterator &points_begin,
				Coord_Iterator &points_end,
				Handle_Iterator &handles_begin,
				Handle_Iterator &handles_end)
{
    if (field_name == "PEAKS")
    {
	// Vector setup.
	points.resize( 3*(b->x_domain().size())*(b->y_domain().size()) );
	handles.resize( points.size() / 3 );

	// Set iterators.
	Coord_Iterator x_it;
	Coord_Iterator y_it;
	Coord_Iterator coord_it = points.begin();
	Handle_Iterator handle_it = handles.begin();

	// Populate the coordinate and handle vectors.
	int handle_counter = 0;
	for (y_it = b->y_domain().begin(); 
	     y_it != b->y_domain().end(); 
	     ++y_it)
	{
	    for (x_it = b->x_domain().begin(); 
		 x_it != b->x_domain().end(); 
		 ++x_it)
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

	// Return the iterators.
	points_begin = points.begin();
	points_end = points.end();
	handles_begin = handles.begin();
	handles_end = handles.end();
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

} // end namespace peaks_example

//---------------------------------------------------------------------------//
//                 end of TE_Physics_B.cc
//---------------------------------------------------------------------------//
