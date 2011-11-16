//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   implemenation/TE_Physics_B.cc
 * \author Stuart Slattery
 * \date   Wed Oct 05 10:57:33 2011
 * \brief  Transfer Evaluator for Physics_B.
 */
//---------------------------------------------------------------------------//
// $Id: template.cc,v 1.3 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#include "TE_Physics_B.hh"

namespace coupler
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
//! Register the communicator.
void TE_Physics_B::register_comm(const Communicator &comm)
{
    comm = b->comm();
}

//---------------------------------------------------------------------------//
//! Register a field associated with the entities.
bool TE_Physics_B::register_field(const std::string &field_name)
{
    if (field_name == "ORDINATE") return true;
    else return false;
}

//---------------------------------------------------------------------------//
// Register cartesian coordinates with a field. The coordinate vector
// should be interleaved. The handle vector should consist of globally
// unique handles. These iterators imply contiguous memory storage.
void TE_Physics_B::register_points(const std::string &field_name,
				   std::vector<HandleType> &handles,
				   std::vector<CoordinateType> &coordinates)
{
    if (field_name == "ORDINATE")
    {
	// Vector setup.
	points.resize( 3*(b->x_domain().size())*(b->y_domain().size()) );
	handles.resize( points.size() / 3 );

	// Set iterators.
	std::vector<double>::const_iterator x_it;
	std::vector<double>::const_iterator y_it;
	std::vector<CoordinateType>::const_iterator coord_it = points.begin();
	std::vector<HandleType>::const_iterator handle_it = handles.begin();

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
    }
}

//---------------------------------------------------------------------------//
// Push data onto a field.
void TE_Physics_B::push_data(const std::string &field_name,
			     const std::vector<HandleType> &handles,
			     const std::vector<DataType> &data)
{
    b->set_source(handles, data);
}

//---------------------------------------------------------------------------//

} // end namespace coupler

//---------------------------------------------------------------------------//
//                 end of TE_Physics_B.cc
//---------------------------------------------------------------------------//
