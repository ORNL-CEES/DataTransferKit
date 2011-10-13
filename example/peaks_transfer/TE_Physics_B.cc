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
// Register a vector of interleaved point coordinates.
void TE_Physics_B::register_xyz(std::vector<double> &points)
{

    // Coordinate vector.
    points.resize( 3*(b->x_domain().size())*(b->y_domain().size()) );

    // Set iterators.
    physics_B::Physics_B::Vector_Dbl::const_iterator x_it;
    physics_B::Physics_B::Vector_Dbl::const_iterator y_it;
    std::vector<double>::iterator coord_it = points.begin();
    
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
	}
    }
}

//---------------------------------------------------------------------------//
//! Register a field associated with the entities.
void TE_Physics_B::register_field(std::string field_name)
{

}

//---------------------------------------------------------------------------//
// Given a (x,y,z) coordinates, return the local process rank in
// which that point exists and the index into the local state vector that
// will be applied at that point. Return true if point is in the local
// domain, false if not.
void TE_Physics_B::find_xyz(double x, 
			    double y, 
			    double z, 
			    int &rank,
			    int &index,
			    bool domain)
{
    rank = 0;
    index = 0;
    domain = true;
}

//---------------------------------------------------------------------------//
// Pull data from a field.
void TE_Physics_B::pull_data(int index, double &data)
{

}

//---------------------------------------------------------------------------//
// Push data onto a field.
void TE_Physics_B::push_data(int index, double data)
{
    b->set_source(index, data);
}

//---------------------------------------------------------------------------//

} // end namespace dtransfer

//---------------------------------------------------------------------------//
//                 end of TE_Physics_B.cc
//---------------------------------------------------------------------------//
