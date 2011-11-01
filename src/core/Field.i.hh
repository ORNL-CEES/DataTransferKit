//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   field/Field.i.hh
 * \author Stuart Slattery
 * \date   Thu Oct 13 09:42:01 2011
 * \brief  Member definitions of class Field.
 */
//---------------------------------------------------------------------------//
// $Id: template.i.hh,v 1.4 2008/01/04 22:50:12 9te Exp $
//---------------------------------------------------------------------------//

#ifndef field_Field_i_hh
#define field_Field_i_hh

namespace field
{

//---------------------------------------------------------------------------//
// Constructor.
Field::Field()
{
    // Initialize the maps.
    d_map_a2b = new Transfer_Map<Handle,Ordinate>();
    d_map_b2a = new Transfer_Map<Handle,Ordinate>();
}

//---------------------------------------------------------------------------//
// Destructor.
Field::~Field()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Get the map for A to B.
Field::SP_Transfer_Map map_a2b()
{
    return d_map_a2b;
}

//---------------------------------------------------------------------------//
// Get the map for B to A.
Field::SP_Transfer_Map map_b2a()
{
    return d_map_b2a;
}

//---------------------------------------------------------------------------//

} // end namespace field

#endif // field_Field_i_hh

//---------------------------------------------------------------------------//
//              end of field/Field.i.hh
//---------------------------------------------------------------------------//
